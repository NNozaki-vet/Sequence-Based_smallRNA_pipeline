# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: initializedcheck=False
# cython: cdivision=True

from cpython.bytes cimport PyBytes_AS_STRING, PyBytes_GET_SIZE
from libc.string cimport memcmp
cimport cython

cdef int MAX_MM_SUPPORTED = 3

cdef inline int _resolve_group_rank(
    int mm_count,
    int rank0,
    int rank1,
    int rank2,
    int rank3,
) noexcept:
    if mm_count == 0:
        return rank0
    elif mm_count == 1:
        return rank1
    elif mm_count == 2:
        return rank2
    elif mm_count == 3:
        return rank3
    return -1


cdef inline tuple _build_hit_tuple(
    int payload_id,
    int del5,
    int mismatch_count,
    int p1,
    int rb1,
    int qb1,
    int p2,
    int rb2,
    int qb2,
    int p3,
    int rb3,
    int qb3,
):
    return (
        payload_id,
        del5,
        mismatch_count,
        p1, rb1, qb1,
        p2, rb2, qb2,
        p3, rb3, qb3,
    )


cdef inline int _count_mismatches_ptr(
    const unsigned char* ref_ptr,
    Py_ssize_t start,
    const unsigned char* core_ptr,
    Py_ssize_t core_len,
    int max_mm,
    int prefix_check_len,
    int* p1,
    int* rb1,
    int* qb1,
    int* p2,
    int* rb2,
    int* qb2,
    int* p3,
    int* rb3,
    int* qb3,
) noexcept nogil:
    cdef Py_ssize_t i
    cdef Py_ssize_t prefix_len
    cdef int mm_count = 0
    cdef unsigned char rb, qb

    p1[0] = 0; rb1[0] = 0; qb1[0] = 0
    p2[0] = 0; rb2[0] = 0; qb2[0] = 0
    p3[0] = 0; rb3[0] = 0; qb3[0] = 0

    if core_len > 0 and memcmp(<const void*>(ref_ptr + start), <const void*>core_ptr, core_len) == 0:
        return 0

    if prefix_check_len > 0:
        prefix_len = prefix_check_len
        if prefix_len > core_len:
            prefix_len = core_len

        for i in range(prefix_len):
            rb = ref_ptr[start + i]
            qb = core_ptr[i]
            if rb != qb:
                mm_count += 1
                if mm_count == 1:
                    p1[0] = <int>(i + 1); rb1[0] = <int>rb; qb1[0] = <int>qb
                elif mm_count == 2:
                    p2[0] = <int>(i + 1); rb2[0] = <int>rb; qb2[0] = <int>qb
                elif mm_count == 3:
                    p3[0] = <int>(i + 1); rb3[0] = <int>rb; qb3[0] = <int>qb
                if mm_count > max_mm:
                    return -1

        for i in range(prefix_len, core_len):
            rb = ref_ptr[start + i]
            qb = core_ptr[i]
            if rb != qb:
                mm_count += 1
                if mm_count == 1:
                    p1[0] = <int>(i + 1); rb1[0] = <int>rb; qb1[0] = <int>qb
                elif mm_count == 2:
                    p2[0] = <int>(i + 1); rb2[0] = <int>rb; qb2[0] = <int>qb
                elif mm_count == 3:
                    p3[0] = <int>(i + 1); rb3[0] = <int>rb; qb3[0] = <int>qb
                if mm_count > max_mm:
                    return -1
    else:
        for i in range(core_len):
            rb = ref_ptr[start + i]
            qb = core_ptr[i]
            if rb != qb:
                mm_count += 1
                if mm_count == 1:
                    p1[0] = <int>(i + 1); rb1[0] = <int>rb; qb1[0] = <int>qb
                elif mm_count == 2:
                    p2[0] = <int>(i + 1); rb2[0] = <int>rb; qb2[0] = <int>qb
                elif mm_count == 3:
                    p3[0] = <int>(i + 1); rb3[0] = <int>rb; qb3[0] = <int>qb
                if mm_count > max_mm:
                    return -1

    return mm_count


cdef tuple _scan_packed_group_at_fixed_start(
    const unsigned char* ref_ptr,
    Py_ssize_t ref_len,
    object group_obj,
    Py_ssize_t fixed_start,
    int max_mm,
    int current_best_group_rank,
    int prefix_check_len,
):
    """
    Compares a packed group at the fixed reference position fixed_start.

    packed group layout:
      group_obj[0] = group_best_rank
      group_obj[1] = core_len
      group_obj[2] = n_payloads
      group_obj[3] = payload_ids      array('I')
      group_obj[4] = best_ranks       array('I')
      group_obj[5] = rank0            array('I')
      group_obj[6] = rank1            array('I')
      group_obj[7] = rank2            array('I')
      group_obj[8] = rank3            array('I')
      group_obj[9] = core_blob_bytes  bytes
    """
    cdef list best_hits = []
    cdef int local_best_rank = current_best_group_rank

    cdef int group_best_rank
    cdef Py_ssize_t core_len
    cdef Py_ssize_t n_payloads

    cdef unsigned int[::1] payload_ids
    cdef unsigned int[::1] best_ranks
    cdef unsigned int[::1] rank0
    cdef unsigned int[::1] rank1
    cdef unsigned int[::1] rank2
    cdef unsigned int[::1] rank3

    cdef bytes core_blob_bytes
    cdef const unsigned char* core_blob_ptr

    cdef Py_ssize_t i
    cdef Py_ssize_t core_off
    cdef const unsigned char* core_ptr

    cdef int mm_count
    cdef int group_rank
    cdef int p1, rb1, qb1
    cdef int p2, rb2, qb2
    cdef int p3, rb3, qb3

    if not group_obj:
        return (-1, best_hits)

    group_best_rank = <int>group_obj[0]
    core_len = <Py_ssize_t>group_obj[1]
    n_payloads = <Py_ssize_t>group_obj[2]

    if n_payloads <= 0:
        return (-1, best_hits)

    if local_best_rank >= 0 and group_best_rank > local_best_rank:
        return (-1, best_hits)

    if fixed_start < 0:
        return (-1, best_hits)

    if core_len <= 0 or core_len > ref_len:
        return (-1, best_hits)

    if fixed_start + core_len > ref_len:
        return (-1, best_hits)

    payload_ids = group_obj[3]
    best_ranks = group_obj[4]
    rank0 = group_obj[5]
    rank1 = group_obj[6]
    rank2 = group_obj[7]
    rank3 = group_obj[8]

    core_blob_bytes = <bytes>group_obj[9]
    core_blob_ptr = <const unsigned char*>PyBytes_AS_STRING(core_blob_bytes)

    for i in range(n_payloads):
        if local_best_rank >= 0 and <int>best_ranks[i] > local_best_rank:
            continue

        core_off = i * core_len
        core_ptr = core_blob_ptr + core_off

        mm_count = _count_mismatches_ptr(
            ref_ptr,
            fixed_start,
            core_ptr,
            core_len,
            max_mm,
            prefix_check_len,
            &p1, &rb1, &qb1,
            &p2, &rb2, &qb2,
            &p3, &rb3, &qb3,
        )

        if mm_count < 0:
            continue

        group_rank = _resolve_group_rank(
            mm_count,
            <int>rank0[i],
            <int>rank1[i],
            <int>rank2[i],
            <int>rank3[i],
        )

        if group_rank < 0:
            continue

        if local_best_rank >= 0 and group_rank > local_best_rank:
            continue

        if local_best_rank < 0 or group_rank < local_best_rank:
            local_best_rank = group_rank
            best_hits = []

        if group_rank == local_best_rank:
            best_hits.append(
                _build_hit_tuple(
                    <int>payload_ids[i],
                    <int>fixed_start,
                    mm_count,
                    p1, rb1, qb1,
                    p2, rb2, qb2,
                    p3, rb3, qb3,
                )
            )

    return (local_best_rank, best_hits)


cdef tuple _scan_both_exactlen_packed(
    bytes ref_bytes,
    object packed_group,
    int max_mm,
    int current_best_group_rank,
    int prefix_check_len,
):
    cdef Py_ssize_t ref_len = PyBytes_GET_SIZE(ref_bytes)
    cdef const unsigned char* ref_ptr = <const unsigned char*>PyBytes_AS_STRING(ref_bytes)
    cdef Py_ssize_t core_len

    if not packed_group:
        return (-1, [])

    core_len = <Py_ssize_t>packed_group[1]
    if ref_len != core_len:
        return (-1, [])

    return _scan_packed_group_at_fixed_start(
        ref_ptr,
        ref_len,
        packed_group,
        0,
        max_mm,
        current_best_group_rank,
        prefix_check_len,
    )

cdef tuple _scan_left_groups_packed(
    bytes ref_bytes,
    object packed_groups,
    int max_mm,
    int current_best_group_rank,
    int prefix_check_len,
):
    """
    For LEFT_ONLY.
    Payloads with ext5 > 0 are fixed to the left end of the reference and are therefore compared at start = 0.
    """
    cdef Py_ssize_t ref_len = PyBytes_GET_SIZE(ref_bytes)
    cdef const unsigned char* ref_ptr = <const unsigned char*>PyBytes_AS_STRING(ref_bytes)

    cdef list best_hits = []
    cdef int local_best_rank = current_best_group_rank

    cdef object group_obj
    cdef tuple subres
    cdef int sub_best
    cdef list sub_hits
    cdef Py_ssize_t core_len

    if not packed_groups:
        return (-1, best_hits)

    for group_obj in packed_groups:
        if not group_obj:
            continue

        core_len = <Py_ssize_t>group_obj[1]
        if core_len <= 0 or core_len > ref_len:
            continue

        subres = _scan_packed_group_at_fixed_start(
            ref_ptr,
            ref_len,
            group_obj,
            0,
            max_mm,
            local_best_rank,
            prefix_check_len,
        )

        sub_best = <int>subres[0]
        sub_hits = <list>subres[1]

        if sub_best >= 0:
            if local_best_rank < 0 or sub_best < local_best_rank:
                local_best_rank = sub_best
                best_hits = sub_hits
            elif sub_best == local_best_rank:
                best_hits.extend(sub_hits)

    return (local_best_rank, best_hits)

cdef tuple _scan_right_groups_packed(
    bytes ref_bytes,
    object packed_groups,
    int max_mm,
    int current_best_group_rank,
    int prefix_check_len,
):
    """
    For RIGHT_ONLY.
    Payloads with ext3 > 0 are fixed to the right end of the reference,
    and are therefore compared at start = ref_len - core_len.
    """
    cdef Py_ssize_t ref_len = PyBytes_GET_SIZE(ref_bytes)
    cdef const unsigned char* ref_ptr = <const unsigned char*>PyBytes_AS_STRING(ref_bytes)

    cdef list best_hits = []
    cdef int local_best_rank = current_best_group_rank

    cdef object group_obj
    cdef tuple subres
    cdef int sub_best
    cdef list sub_hits
    cdef Py_ssize_t core_len
    cdef Py_ssize_t fixed_start

    if not packed_groups:
        return (-1, best_hits)

    for group_obj in packed_groups:
        if not group_obj:
            continue

        core_len = <Py_ssize_t>group_obj[1]
        if core_len <= 0 or core_len > ref_len:
            continue

        fixed_start = ref_len - core_len

        subres = _scan_packed_group_at_fixed_start(
            ref_ptr,
            ref_len,
            group_obj,
            fixed_start,
            max_mm,
            local_best_rank,
            prefix_check_len,
        )

        sub_best = <int>subres[0]
        sub_hits = <list>subres[1]

        if sub_best >= 0:
            if local_best_rank < 0 or sub_best < local_best_rank:
                local_best_rank = sub_best
                best_hits = sub_hits
            elif sub_best == local_best_rank:
                best_hits.extend(sub_hits)

    return (local_best_rank, best_hits)




cdef tuple _scan_anywhere_packed(
    bytes ref_bytes,
    object packed_groups,
    int max_mm,
    int current_best_group_rank,
    int prefix_check_len,
):
    """
    For ANYWHERE.
    Payloads with ext5 = ext3 = 0 are searched at every start position within the reference.
    """
    cdef Py_ssize_t ref_len = PyBytes_GET_SIZE(ref_bytes)
    cdef const unsigned char* ref_ptr = <const unsigned char*>PyBytes_AS_STRING(ref_bytes)

    cdef list best_hits = []
    cdef int local_best_rank = current_best_group_rank

    cdef object group_obj
    cdef int group_best_rank
    cdef Py_ssize_t core_len
    cdef Py_ssize_t n_payloads

    cdef unsigned int[::1] payload_ids
    cdef unsigned int[::1] best_ranks
    cdef unsigned int[::1] rank0
    cdef unsigned int[::1] rank1
    cdef unsigned int[::1] rank2
    cdef unsigned int[::1] rank3

    cdef bytes core_blob_bytes
    cdef const unsigned char* core_blob_ptr

    cdef Py_ssize_t max_start
    cdef Py_ssize_t start
    cdef Py_ssize_t i
    cdef Py_ssize_t core_off
    cdef const unsigned char* core_ptr

    cdef int mm_count
    cdef int group_rank
    cdef int p1, rb1, qb1
    cdef int p2, rb2, qb2
    cdef int p3, rb3, qb3

    if not packed_groups:
        return (-1, best_hits)

    for group_obj in packed_groups:
        if not group_obj:
            continue

        group_best_rank = <int>group_obj[0]
        if local_best_rank >= 0 and group_best_rank > local_best_rank:
            continue

        core_len = <Py_ssize_t>group_obj[1]
        n_payloads = <Py_ssize_t>group_obj[2]

        if n_payloads <= 0:
            continue

        if core_len <= 0 or core_len > ref_len:
            continue

        payload_ids = group_obj[3]
        best_ranks = group_obj[4]
        rank0 = group_obj[5]
        rank1 = group_obj[6]
        rank2 = group_obj[7]
        rank3 = group_obj[8]

        core_blob_bytes = <bytes>group_obj[9]
        core_blob_ptr = <const unsigned char*>PyBytes_AS_STRING(core_blob_bytes)

        max_start = ref_len - core_len

        for i in range(n_payloads):
            if local_best_rank >= 0 and <int>best_ranks[i] > local_best_rank:
                continue

            core_off = i * core_len
            core_ptr = core_blob_ptr + core_off

            for start in range(max_start + 1):
                mm_count = _count_mismatches_ptr(
                    ref_ptr,
                    start,
                    core_ptr,
                    core_len,
                    max_mm,
                    prefix_check_len,
                    &p1, &rb1, &qb1,
                    &p2, &rb2, &qb2,
                    &p3, &rb3, &qb3,
                )

                if mm_count < 0:
                    continue

                group_rank = _resolve_group_rank(
                    mm_count,
                    <int>rank0[i],
                    <int>rank1[i],
                    <int>rank2[i],
                    <int>rank3[i],
                )

                if group_rank < 0:
                    continue

                if local_best_rank >= 0 and group_rank > local_best_rank:
                    continue

                if local_best_rank < 0 or group_rank < local_best_rank:
                    local_best_rank = group_rank
                    best_hits = []

                if group_rank == local_best_rank:
                    best_hits.append(
                        _build_hit_tuple(
                            <int>payload_ids[i],
                            <int>start,
                            mm_count,
                            p1, rb1, qb1,
                            p2, rb2, qb2,
                            p3, rb3, qb3,
                        )
                    )

    return (local_best_rank, best_hits)

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef tuple find_best_hits_multi_payloads_cy(
    bytes ref_bytes,
    object payloads_both_exactlen,
    object payloads_left_by_corelen,
    object payloads_right_by_corelen,
    object anywhere_by_corelen,
    int max_mm,
    int current_best_group_rank=-1,
    int prefix_check_len=0,
):
    """
    Public function called from Python.

    All inputs are expected to use the packed-group format:
      - payloads_both_exactlen
      - payloads_left_by_corelen
      - payloads_right_by_corelen
      - anywhere_by_corelen
    """
    cdef int local_best_rank = current_best_group_rank
    cdef list best_hits = []

    cdef tuple subres
    cdef int sub_best
    cdef list sub_hits

    if max_mm < 0:
        return (-1, best_hits)

    if max_mm > MAX_MM_SUPPORTED:
        raise ValueError("max_mm exceeds supported range in this packed Cython module")

    # BOTH_ENDS
    if payloads_both_exactlen:
        subres = _scan_both_exactlen_packed(
            ref_bytes,
            payloads_both_exactlen,
            max_mm,
            local_best_rank,
            prefix_check_len,
        )

        sub_best = <int>subres[0]
        sub_hits = <list>subres[1]

        if sub_best >= 0:
            if local_best_rank < 0 or sub_best < local_best_rank:
                local_best_rank = sub_best
                best_hits = sub_hits
            elif sub_best == local_best_rank:
                best_hits.extend(sub_hits)

    # LEFT_ONLY
    if payloads_left_by_corelen:
        subres = _scan_left_groups_packed(
            ref_bytes,
            payloads_left_by_corelen,
            max_mm,
            local_best_rank,
            prefix_check_len,
        )

        sub_best = <int>subres[0]
        sub_hits = <list>subres[1]

        if sub_best >= 0:
            if local_best_rank < 0 or sub_best < local_best_rank:
                local_best_rank = sub_best
                best_hits = sub_hits
            elif sub_best == local_best_rank:
                best_hits.extend(sub_hits)

    # RIGHT_ONLY
    if payloads_right_by_corelen:
        subres = _scan_right_groups_packed(
            ref_bytes,
            payloads_right_by_corelen,
            max_mm,
            local_best_rank,
            prefix_check_len,
        )

        sub_best = <int>subres[0]
        sub_hits = <list>subres[1]

        if sub_best >= 0:
            if local_best_rank < 0 or sub_best < local_best_rank:
                local_best_rank = sub_best
                best_hits = sub_hits
            elif sub_best == local_best_rank:
                best_hits.extend(sub_hits)

    # ANYWHERE
    if anywhere_by_corelen:
        subres = _scan_anywhere_packed(
            ref_bytes,
            anywhere_by_corelen,
            max_mm,
            local_best_rank,
            prefix_check_len,
        )

        sub_best = <int>subres[0]
        sub_hits = <list>subres[1]

        if sub_best >= 0:
            if local_best_rank < 0 or sub_best < local_best_rank:
                local_best_rank = sub_best
                best_hits = sub_hits
            elif sub_best == local_best_rank:
                best_hits.extend(sub_hits)

    return (local_best_rank, best_hits)