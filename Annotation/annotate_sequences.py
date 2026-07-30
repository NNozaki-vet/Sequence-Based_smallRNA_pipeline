#!/usr/bin/env python3
r"""
Parallel exact-sequence annotation against user-defined reference FASTA files.

The computational search, split-rank logic, dynamic batch sizing, reference
caching, multiprocessing scheduler, and minimum-score output procedure are
preserved from the research analysis implementation. Dataset-specific paths
and reference definitions are supplied at runtime.

Example:
    python annotate_sequences.py ^
        --query-fasta "D:\project\filtered_sequences.fa" ^
        --references-file "D:\project\references.csv" ^
        --output-dir "D:\project\annotation_results" ^
        --core-name "project_sequences" ^
        --workers 16
"""

import atexit
import argparse
import bisect
import csv
import multiprocessing as mp
import os
import pickle
import shutil
import sys
import tempfile
import time
import math
from array import array
from collections import defaultdict, deque
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


def set_csv_field_size_limit_to_max() -> None:
    max_size = sys.maxsize
    while True:
        try:
            csv.field_size_limit(max_size)
            break
        except OverflowError:
            max_size //= 10


set_csv_field_size_limit_to_max()

# Compiled annotation core
#
# The extension module can be provided in either of the following ways:
#   1. Install/build it from annotation_align_core.pyx.
#   2. Build it from annotation_align_core.c.
#   3. Place a compatible prebuilt annotation_align_core*.pyd beside this
#      script or install it in the active Python environment.
try:
    from annotation_align_core import find_best_hits_multi_payloads_cy
except ImportError as exc:
    raise ImportError(
        "Could not import 'annotation_align_core'. Build or install the "
        "compiled extension before running this script. A prebuilt "
        "annotation_align_core.cp313-win_amd64.pyd is compatible only with "
        "64-bit CPython 3.13 on Windows."
    ) from exc


# =========================================================
# Runtime input settings
# =========================================================

reference_list = []
sequence_fasta = ""
core_name = ""
output_folder = ""
CACHE_DIR = ""


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Annotate exact query sequences against one or more reference "
            "FASTA files using the compiled annotation_align_core extension."
        )
    )
    parser.add_argument(
        "--query-fasta",
        required=True,
        help=(
            "FASTA file containing query sequences. Duplicate query "
            "sequences are collapsed while preserving first-occurrence order."
        ),
    )
    parser.add_argument(
        "--references-file",
        required=True,
        help=(
            "CSV file with the required columns 'reference_name' and "
            "'file_path'. Relative reference paths are resolved from the "
            "directory containing this CSV file."
        ),
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory in which annotation CSV files and reference caches are written.",
    )
    parser.add_argument(
        "--core-name",
        default=None,
        help=(
            "Prefix used for output filenames. By default, the query FASTA "
            "filename is used after removing its FASTA extension."
        ),
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=16,
        help="Number of multiprocessing workers. Default: 16.",
    )
    return parser.parse_args()


def remove_fasta_suffix(filename):
    lower_name = filename.lower()
    for suffix in (
        ".fasta.gz",
        ".fna.gz",
        ".fa.gz",
        ".fasta",
        ".fna",
        ".fa",
    ):
        if lower_name.endswith(suffix):
            return filename[:-len(suffix)]
    return Path(filename).stem


def load_reference_definitions(references_file):
    references_file = Path(references_file).expanduser()

    if not references_file.is_file():
        raise FileNotFoundError(
            f"Reference definition CSV not found: {references_file}"
        )

    references = []

    with references_file.open(
        "r",
        newline="",
        encoding="utf-8-sig",
    ) as handle:
        reader = csv.DictReader(handle)

        if reader.fieldnames is None:
            raise ValueError(
                f"Reference definition CSV has no header: {references_file}"
            )

        reader.fieldnames = [
            str(name).replace("\ufeff", "").strip()
            for name in reader.fieldnames
        ]

        required_columns = {"reference_name", "file_path"}
        missing_columns = required_columns.difference(reader.fieldnames)

        if missing_columns:
            raise ValueError(
                "Reference definition CSV is missing required columns: "
                + ", ".join(sorted(missing_columns))
            )

        for row_number, row in enumerate(reader, start=2):
            reference_name = str(row["reference_name"]).strip()
            file_path_text = str(row["file_path"]).strip()

            if not reference_name and not file_path_text:
                continue

            if not reference_name:
                raise ValueError(
                    f"Empty reference_name at row {row_number}: "
                    f"{references_file}"
                )

            if not file_path_text:
                raise ValueError(
                    f"Empty file_path at row {row_number}: "
                    f"{references_file}"
                )

            file_path = Path(file_path_text).expanduser()
            if not file_path.is_absolute():
                file_path = references_file.parent / file_path

            references.append(
                {
                    "reference_name": reference_name,
                    "file_path": str(file_path.resolve()),
                }
            )

    if not references:
        raise ValueError(
            f"No reference definitions were found in: {references_file}"
        )

    reference_names = [item["reference_name"] for item in references]
    duplicated_names = sorted({
        name
        for name in reference_names
        if reference_names.count(name) > 1
    })

    if duplicated_names:
        raise ValueError(
            "Duplicate reference_name values were found: "
            + ", ".join(duplicated_names)
        )

    return references


def configure_runtime(args):
    global reference_list
    global sequence_fasta
    global core_name
    global output_folder
    global CACHE_DIR
    global N_WORKERS
    global DOMINANT_FIRST_MAX_SLOTS

    query_path = Path(args.query_fasta).expanduser()
    output_path = Path(args.output_dir).expanduser()

    reference_list = load_reference_definitions(args.references_file)
    sequence_fasta = str(query_path.resolve())
    output_folder = str(output_path.resolve())

    if args.core_name is None:
        core_name = remove_fasta_suffix(query_path.name)
    else:
        core_name = str(args.core_name).strip()

    if not core_name:
        raise ValueError("core_name must not be empty.")

    if any(character in core_name for character in '<>:"/\\|?*'):
        raise ValueError(
            "core_name contains a character that is not valid in a "
            f"Windows filename: {core_name!r}"
        )

    if args.workers < 1:
        raise ValueError("--workers must be >= 1.")

    N_WORKERS = int(args.workers)
    DOMINANT_FIRST_MAX_SLOTS = N_WORKERS

    os.makedirs(output_folder, exist_ok=True)
    CACHE_DIR = os.path.join(output_folder, "reference_cache")
    os.makedirs(CACHE_DIR, exist_ok=True)


# =========================================================
# Speed / output settings
# =========================================================
N_WORKERS = 16
BATCH_SIZE = 150                  # Fallback value used only when USE_DYNAMIC_BATCH_SIZE=False.
BATCHES_PER_TASK = 1
KEEP_TEMP_FILES = False
UNANNOTATED_TOTAL_SCORE = 10

# B. safe reject
PREFIX_CHECK_LEN = 6              # Start with 4; it may later be compared with 6.

# A. dynamic batch sizing
USE_DYNAMIC_BATCH_SIZE = True

MIN_WAVES_PER_REFERENCE = 0        # Ensure at least N_WORKERS×0 tasks for each reference.
MAX_WAVES_PER_REFERENCE = 40       # In principle, limit the number to at most N_WORKERS×40 tasks.
MAX_TASK_SEC_HARD = 5400           # Approximate upper limit when one task is too long; decrease it to shorten task time (e.g., 3600), or increase it if too many tasks are created (e.g., 7200).
PILOT_QUERY_COUNT = 16             # Select 16 queries at evenly spaced intervals from all queries.
BATCH_SAFETY_FACTOR = 0.80         # Effective target seconds = TARGET_TASK_SEC × 0.80

MIN_BATCH_SIZE = 40
MAX_BATCH_SIZE = 100000

# Auto target task seconds
USE_AUTO_TARGET_TASK_SEC = True

TARGET_TASK_SEC_FALLBACK = 1800    # Used only when USE_AUTO_TARGET_TASK_SEC=False.
TARGET_TASK_SEC_BASE = 1200        # Target seconds at the reference query count.
TARGET_TASK_SEC_QUERY_REF = 2000   # Reference query count.
TARGET_TASK_SEC_BASE_WORKERS = 16  # Reference worker count.

TARGET_TASK_SEC_MIN = 900          # Do not make tasks shorter than this.
TARGET_TASK_SEC_MAX = 5400         # Do not make target tasks longer than this.
TARGET_TASK_SEC_ROUND_TO = 60      # Round to 60-second increments.

TARGET_TASK_SEC_QUERY_EXP = 0.25   # Rate of increase with query count.
TARGET_TASK_SEC_WORKER_EXP = 0.50  # Adjustment based on worker count.

# Task ordering
#   heavy_first tends to place only reference_hsa first, which can leave other references until the end.
#   weighted_interleave allocates more hsa tasks according to the predicted total workload by reference,
#   while also scheduling mmu/bta/cfa from the early stage.
TASK_ORDER_STRATEGY = "weighted_interleave"
WEIGHTED_INTERLEAVE_WEIGHT_EXP = 1.0   # 1.0: proportional to predicted total workload. 0.5: more evenly balanced.
WEIGHTED_INTERLEAVE_MIN_WEIGHT = 1
WEIGHTED_INTERLEAVE_MAX_WEIGHT = 20
TASK_ORDER_PREVIEW_N = 30

# Dominant-reference first wave
#   If the reference with the largest predicted total workload accounts for most of the overall workload,
#   prioritize that reference in the first wave to fill the CPUs.
#   Then mix in other references using weighted interleave.
#   This shortens the critical path for a dominant reference such as hsa,
#   while also reducing the tail effect in which bta/cfa/mmu remain until the end.
DOMINANT_FIRST_WAVE = True
DOMINANT_FIRST_MIN_FRACTION = 0.50
DOMINANT_FIRST_MAX_SLOTS = N_WORKERS
# =========================================================
# Split-rank settings
# =========================================================
MAX_EXT_TOTAL = 3
MAX_MM = 3
MAX_MISMATCH = MAX_MM
RANK_MAP = None

OUTPUT_COLUMNS = [
    "query_sequence",
    "reference_entry_name",
    "total_mismatch_and_extension_count",
    "mismatch_count",
    "query_5prime_extension_length",
    "query_5prime_extension_sequence",
    "reference_5prime_unmatched_length",
    "reference_mismatch_position1_1based",
    "reference_to_query_mismatch1",
    "reference_mismatch_position2_1based",
    "reference_to_query_mismatch2",
    "reference_mismatch_position3_1based",
    "reference_to_query_mismatch3",
    "query_3prime_extension_length",
    "query_3prime_extension_sequence",
    "reference_3prime_unmatched_length",
]

# Worker globals
G_CURRENT_REFERENCE_KEY = None
G_REF_RECORDS = None
G_REF_LENGTHS = None
G_START_IDX_CACHE = None
G_OPEN_TEMP_REFERENCE_NAME = None
G_OPEN_TEMP_PATH = None
G_OPEN_TEMP_FILE = None
G_OPEN_TEMP_WRITER = None
G_TEMP_PART_INDEX_BY_REFERENCE = None
G_QUERY_BUNDLE_CACHE = None

REFERENCE_CACHE_FORMAT_VERSION = 2

# ---------------------------
# Basic IO / helpers
# ---------------------------
def clean_seq(seq: str) -> str:
    return str(seq).upper().replace("U", "T").replace(" ", "").replace("\n", "").replace("\r", "")


def read_fasta_sequences(fasta_path: str) -> List[str]:
    seqs: List[str] = []
    current: List[str] = []
    with open(fasta_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    seqs.append(clean_seq("".join(current)))
                    current = []
            else:
                current.append(line)
        if current:
            seqs.append(clean_seq("".join(current)))
    seen = set()
    uniq = []
    for s in seqs:
        if s not in seen:
            seen.add(s)
            uniq.append(s)
    return uniq

def validate_input_paths() -> None:
    if not os.path.isfile(sequence_fasta):
        raise FileNotFoundError(f"Query FASTA not found: {sequence_fasta}")

    if not reference_list:
        raise ValueError("reference_list is empty.")

    for ref in reference_list:
        ref_name = ref.get("reference_name", "")
        ref_path = ref.get("file_path", "")

        if not ref_name:
            raise ValueError(f"reference_name is missing in reference_list entry: {ref}")

        if not ref_path:
            raise ValueError(f"file_path is missing for reference: {ref_name}")

        if not os.path.isfile(ref_path):
            raise FileNotFoundError(f"Reference FASTA not found for {ref_name}: {ref_path}")


def chunk_list(lst: Sequence[str], n: int) -> Iterable[List[str]]:
    for i in range(0, len(lst), n):
        yield list(lst[i:i + n])


def get_reference_cache_path(reference_name: str) -> str:
    return os.path.join(CACHE_DIR, f"{core_name}_{reference_name}_preprocessed.pkl")


def load_reference_fasta(reference_fasta_path: str):
    """
    Load a reference FASTA file.

    Expected format:
      >reference_entry_name
      ATCG...

    In a conventional two-line FASTA format, the header is on an odd-numbered line,
    and the sequence is on an even-numbered line.

    The implementation processes FASTA records rather than physical line numbers.
    Therefore, line-wrapped FASTA sequences are also supported.
    """
    seq_to_names: Dict[str, List[str]] = {}

    current_entry_name = None
    current_seq_parts: List[str] = []

    def register_current_record() -> None:
        if current_entry_name is None:
            return

        ref_seq = clean_seq("".join(current_seq_parts))

        if not current_entry_name:
            raise ValueError(
                f"Empty FASTA header detected in reference file: "
                f"{reference_fasta_path}"
            )

        if not ref_seq:
            raise ValueError(
                f"Empty sequence detected for FASTA entry "
                f"{current_entry_name!r}: {reference_fasta_path}"
            )

        seq_to_names.setdefault(ref_seq, []).append(current_entry_name)

    with open(reference_fasta_path, "r", encoding="utf-8-sig") as f:
        for line_no, raw_line in enumerate(f, start=1):
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if current_entry_name is not None:
                    register_current_record()

                current_entry_name = line[1:].strip()
                current_seq_parts = []

                if not current_entry_name:
                    raise ValueError(
                        f"Empty FASTA header at line {line_no}: "
                        f"{reference_fasta_path}"
                    )
            else:
                if current_entry_name is None:
                    raise ValueError(
                        f"Sequence found before the first FASTA header "
                        f"at line {line_no}: {reference_fasta_path}"
                    )

                current_seq_parts.append(line)

    if current_entry_name is not None:
        register_current_record()

    if not seq_to_names:
        raise ValueError(
            f"No valid reference sequences found in FASTA: "
            f"{reference_fasta_path}"
        )

    records = []

    for ref_seq, entry_names in seq_to_names.items():
        ref_bytes = ref_seq.encode("ascii")
        ref_len = len(ref_seq)
        records.append((ref_seq, ref_bytes, entry_names, ref_len))

    records.sort(key=lambda x: (x[3], x[0]))
    ref_lengths = [r[3] for r in records]

    return records, ref_lengths



def get_file_metadata(file_path: str):
    st = os.stat(file_path)
    return {
        "source_path": os.path.abspath(file_path),
        "source_mtime": st.st_mtime,
        "source_size": st.st_size,
    }


def is_cache_valid(cache_payload: dict, reference_fasta_path: str) -> bool:
    if not isinstance(cache_payload, dict):
        return False

    current_meta = get_file_metadata(reference_fasta_path)

    return (
        cache_payload.get("cache_format_version")
        == REFERENCE_CACHE_FORMAT_VERSION
        and cache_payload.get("source_path") == current_meta["source_path"]
        and cache_payload.get("source_mtime") == current_meta["source_mtime"]
        and cache_payload.get("source_size") == current_meta["source_size"]
        and isinstance(cache_payload.get("records"), list)
        and isinstance(cache_payload.get("ref_lengths"), list)
    )


def load_or_build_reference_cache(reference_name: str, reference_fasta_path: str):
    cache_path = get_reference_cache_path(reference_name)

    if os.path.exists(cache_path):
        try:
            with open(cache_path, "rb") as f:
                payload = pickle.load(f)

            if is_cache_valid(payload, reference_fasta_path):
                return (
                    cache_path,
                    payload["records"],
                    payload["ref_lengths"],
                    True,
                )
        except Exception:
            pass

    records, ref_lengths = load_reference_fasta(reference_fasta_path)

    payload = get_file_metadata(reference_fasta_path)
    payload["cache_format_version"] = REFERENCE_CACHE_FORMAT_VERSION
    payload["records"] = records
    payload["ref_lengths"] = ref_lengths

    tmp_cache_path = cache_path + ".tmp"

    with open(tmp_cache_path, "wb") as f:
        pickle.dump(payload, f, protocol=pickle.HIGHEST_PROTOCOL)

    os.replace(tmp_cache_path, cache_path)

    return cache_path, records, ref_lengths, False


def prepare_all_reference_caches(reference_defs):
    cache_info_by_reference = {}
    for ref in reference_defs:
        ref_name = ref["reference_name"]
        ref_path = ref["file_path"]
        cache_path, records, ref_lengths, cache_hit = load_or_build_reference_cache(ref_name, ref_path)
        cache_info_by_reference[ref_name] = {
            "cache_path": cache_path,
            "n_unique_reference_sequences": len(records),
            "cache_hit": cache_hit,
        }
        status = "cache reused" if cache_hit else "cache rebuilt"
        print(f"[{ref_name}] {status}: {cache_path} (unique reference sequences: {len(records)})", flush=True)
    return cache_info_by_reference


def load_reference_cache(cache_path: str):
    with open(cache_path, "rb") as f:
        payload = pickle.load(f)
    return payload["records"], payload["ref_lengths"]


# ---------------------------
# Output row builders
# ---------------------------
def build_output_row_list_fast(sequence, reference_entry_name, hit=None, total_override=None):
    row = [sequence, reference_entry_name, 0, 0, 0, "", 0, "", "", "", "", "", "", 0, "", 0]
    if hit is None:
        if total_override is not None:
            row[2] = total_override
        return row

    ext5, ext5_seq, del5, mismatches, ext3, ext3_seq, del3, mismatch_count_override = hit
    mismatch_count = mismatch_count_override if mismatch_count_override is not None else len(mismatches)
    total_mismatch_and_extension_count = mismatch_count + ext5 + ext3
    row[2] = total_mismatch_and_extension_count
    row[3] = mismatch_count
    row[4] = ext5
    row[5] = ext5_seq
    row[6] = del5
    for i in range(min(3, len(mismatches))):
        pos, base_change = mismatches[i]
        row[7 + i * 2] = pos
        row[8 + i * 2] = base_change
    row[13] = ext3
    row[14] = ext3_seq
    row[15] = del3
    return row


def build_output_row_from_fixed_hit(sequence: str, reference_entry_name: str, payload: dict, fixed_hit: tuple, rlen: int):
    (
        _payload_id,
        del5,
        mismatch_count,
        p1, rb1, qb1,
        p2, rb2, qb2,
        p3, rb3, qb3,
    ) = fixed_hit

    mismatches = []
    if mismatch_count >= 1 and p1 > 0:
        mismatches.append((p1, f"{chr(rb1)}>{chr(qb1)}"))
    if mismatch_count >= 2 and p2 > 0:
        mismatches.append((p2, f"{chr(rb2)}>{chr(qb2)}"))
    if mismatch_count >= 3 and p3 > 0:
        mismatches.append((p3, f"{chr(rb3)}>{chr(qb3)}"))

    del3 = rlen - (del5 + payload["core_len"])
    full_hit = (
        payload["ext5"],
        payload["ext5_seq"],
        del5,
        tuple(mismatches),
        payload["ext3"],
        payload["ext3_seq"],
        del3,
        mismatch_count,
    )
    return build_output_row_list_fast(
        sequence=sequence,
        reference_entry_name=reference_entry_name,
        hit=full_hit,
    )


# ---------------------------
# Split-rank helpers
# ---------------------------
def build_priority_splits(max_ext: int) -> List[Tuple[int, int]]:
    splits: List[Tuple[int, int]] = []
    for ext_total in range(0, max_ext + 1):
        for ext5 in range(0, ext_total + 1):
            ext3 = ext_total - ext5
            splits.append((ext5, ext3))
    return splits


def build_split_rank_map(max_ext: int, max_mm: int) -> Dict[Tuple[int, int, int], int]:
    rank_map: Dict[Tuple[int, int, int], int] = {}
    rank = 0
    for mm in range(0, max_mm + 1):
        for ext5, ext3 in build_priority_splits(max_ext=max_ext):
            rank_map[(ext5, ext3, mm)] = rank
            rank += 1
    return rank_map


def ensure_rank_map_initialized() -> Dict[Tuple[int, int, int], int]:
    global RANK_MAP
    if RANK_MAP is None:
        RANK_MAP = build_split_rank_map(max_ext=MAX_EXT_TOTAL, max_mm=MAX_MM)
    return RANK_MAP


def build_payload_spec(payload_id, ext5, ext3, core_bytes, rank_map):
    """
    Return a payload specification for packed input.
    Returns:
        (payload_id, core_bytes, core_len, best_possible_rank, rank0, rank1, rank2, rank3)
    """
    core_len = len(core_bytes)

    rank0 = rank_map[(ext5, ext3, 0)]
    rank1 = rank_map[(ext5, ext3, 1)]
    rank2 = rank_map[(ext5, ext3, 2)]
    rank3 = rank_map[(ext5, ext3, 3)]

    best_possible_rank = rank0

    return (
        payload_id,
        core_bytes,
        core_len,
        best_possible_rank,
        rank0,
        rank1,
        rank2,
        rank3,
    )

def sort_specs_for_packing(specs):
    """
    Sort into a stable order before packing.
    key:
      1) best_possible_rank
      2) payload_id
    """
    return sorted(specs, key=lambda s: (s[3], s[0]))


def pack_specs_group(specs, expected_core_len=None):
    """
    Pack specifications with the same core_len into one group.

    Packed group format:
      (
        group_best_rank,
        core_len,
        n_payloads,
        payload_ids,
        best_ranks,
        rank0_arr,
        rank1_arr,
        rank2_arr,
        rank3_arr,
        core_blob_bytes,
      )

    where
      payload_ids, best_ranks, rank0_arr... are array('I') objects
      core_blob_bytes is bytes formed by concatenating the payload cores
    """
    if not specs:
        return None

    specs = sort_specs_for_packing(specs)

    if expected_core_len is None:
        core_len = specs[0][2]
    else:
        core_len = expected_core_len

    payload_ids = array("I")
    best_ranks = array("I")
    rank0_arr = array("I")
    rank1_arr = array("I")
    rank2_arr = array("I")
    rank3_arr = array("I")
    core_blob_parts = []

    for spec in specs:
        payload_id, core_bytes, spec_core_len, best_possible_rank, rank0, rank1, rank2, rank3 = spec

        if spec_core_len != core_len:
            raise ValueError(
                f"pack_specs_group(): core_len mismatch detected "
                f"(expected {core_len}, got {spec_core_len})"
            )

        payload_ids.append(payload_id)
        best_ranks.append(best_possible_rank)
        rank0_arr.append(rank0)
        rank1_arr.append(rank1)
        rank2_arr.append(rank2)
        rank3_arr.append(rank3)
        core_blob_parts.append(core_bytes)

    group_best_rank = min(best_ranks) if len(best_ranks) > 0 else 10**9
    core_blob_bytes = b"".join(core_blob_parts)
    n_payloads = len(payload_ids)

    return (
        group_best_rank,
        core_len,
        n_payloads,
        payload_ids,
        best_ranks,
        rank0_arr,
        rank1_arr,
        rank2_arr,
        rank3_arr,
        core_blob_bytes,
    )


def pack_group_list_from_specs_by_corelen(specs_by_corelen):
    """
    {core_len: [specs]} -> [(group_best_rank, core_len, ...), ...]
    Return in ascending group_best_rank order.
    """
    packed_groups = []

    for core_len, specs in specs_by_corelen.items():
        packed = pack_specs_group(specs, expected_core_len=core_len)
        if packed is not None:
            packed_groups.append(packed)

    packed_groups.sort(key=lambda g: (g[0], g[1]))
    return packed_groups


def pack_both_by_reflen(specs_by_reflen):
    """
    For BOTH_ENDS:
      {ref_len: [specs]} -> {ref_len: packed_group}
    """
    packed = {}
    for ref_len, specs in specs_by_reflen.items():
        if not specs:
            continue
        packed[ref_len] = pack_specs_group(specs)
    return packed


def build_query_payload_bundle_splitrank(query_seq, rank_map, max_ext_total=3, max_mm=3):
    """
    Construct the payload bundle for each query in packed format.

    Strict rules:
      - ext3 > 0          -> anchored to the right end
      - ext5 > 0          -> anchored to the left end
      - ext5 > 0 and ext3 > 0 -> matches both ends (ref_len = len(query) - ext5 - ext3)
      - ext5 = ext3 = 0   -> treated as anywhere

    Rank comparison order:
      1) mm
      2) ext_total = ext5 + ext3
      3) ext5
    """
    qbytes = query_seq.encode("ascii")
    qlen = len(qbytes)

    payload_meta_by_id = {}
    next_payload_id = 0

    both_specs_by_reflen = defaultdict(list)
    left_specs_by_corelen = defaultdict(list)
    right_specs_by_corelen = defaultdict(list)
    anywhere_specs_by_corelen = defaultdict(list)

    min_ref_len_global = qlen    # Use min_ref_len_global as the current search start.
    max_ref_len_global = qlen    # Retain max_ref_len_global for possible future upper limits on reference length.
    for ext_total in range(0, max_ext_total + 1):
        for ext5 in range(0, ext_total + 1):
            ext3 = ext_total - ext5

            core_start = ext5
            core_end = qlen - ext3
            if core_start > core_end:
                continue

            core_bytes = qbytes[core_start:core_end]
            core_len = len(core_bytes)
            if core_len <= 0:
                continue

            payload_id = next_payload_id
            next_payload_id += 1

            payload_meta_by_id[payload_id] = {
                "payload_id": payload_id,
                "ext5": ext5,
                "ext3": ext3,
                "ext5_seq": query_seq[:ext5] if ext5 > 0 else "",
                "ext3_seq": query_seq[qlen - ext3:] if ext3 > 0 else "",
                "core_len": core_len,
                "core_bytes": core_bytes,
            }

            spec = build_payload_spec(
                payload_id=payload_id,
                ext5=ext5,
                ext3=ext3,
                core_bytes=core_bytes,
                rank_map=rank_map,
            )

            ref_len = qlen - ext5 - ext3
            if ref_len < min_ref_len_global:
                min_ref_len_global = ref_len
            if ref_len > max_ref_len_global:
                max_ref_len_global = ref_len

            if ext5 > 0 and ext3 > 0:
                both_specs_by_reflen[ref_len].append(spec)
            elif ext5 > 0:
                left_specs_by_corelen[core_len].append(spec)
            elif ext3 > 0:
                right_specs_by_corelen[core_len].append(spec)
            else:
                anywhere_specs_by_corelen[core_len].append(spec)

    payloads_both_by_reflen = pack_both_by_reflen(both_specs_by_reflen)
    payloads_left_by_corelen = tuple(pack_group_list_from_specs_by_corelen(left_specs_by_corelen))
    payloads_right_by_corelen = tuple(pack_group_list_from_specs_by_corelen(right_specs_by_corelen))
    anywhere_by_corelen = tuple(pack_group_list_from_specs_by_corelen(anywhere_specs_by_corelen))

    return {
        "payload_meta_by_id": payload_meta_by_id,
        "payloads_both_by_reflen": payloads_both_by_reflen,
        "payloads_left_by_corelen": payloads_left_by_corelen,
        "payloads_right_by_corelen": payloads_right_by_corelen,
        "anywhere_by_corelen": anywhere_by_corelen,
        "min_ref_len_global": min_ref_len_global,
        "max_ref_len_global": max_ref_len_global,
        "max_mm": max_mm,
    }
# ----------------------------------------------
#helper functions
# ----------------------------------------------

def clamp(x: int, lo: int, hi: int) -> int:
    if x < lo:
        return lo
    if x > hi:
        return hi
    return x


def round_to_step(x: float, step: int) -> int:
    if step <= 0:
        return int(round(x))
    return int(step * round(float(x) / step))


def choose_target_task_sec(n_query: int, n_workers: int) -> int:
    """
    Automatically determine TARGET_TASK_SEC from the number of queries and workers.

    Basic policy:
      - use longer tasks as the query count increases to prevent excessive task counts
      - use slightly shorter tasks as the worker count increases to maintain parallelism
      - apply gradual rather than abrupt changes
    """
    if not USE_AUTO_TARGET_TASK_SEC:
        return int(TARGET_TASK_SEC_FALLBACK)

    n_query = max(1, int(n_query))
    n_workers = max(1, int(n_workers))

    query_factor = (n_query / TARGET_TASK_SEC_QUERY_REF) ** TARGET_TASK_SEC_QUERY_EXP
    worker_factor = (TARGET_TASK_SEC_BASE_WORKERS / n_workers) ** TARGET_TASK_SEC_WORKER_EXP

    target_sec = TARGET_TASK_SEC_BASE * query_factor * worker_factor
    target_sec = round_to_step(target_sec, TARGET_TASK_SEC_ROUND_TO)

    return clamp(
        int(target_sec),
        TARGET_TASK_SEC_MIN,
        TARGET_TASK_SEC_MAX,
    )


def pick_pilot_queries_evenly(sequence_list: Sequence[str], k: int) -> List[str]:
    n = len(sequence_list)
    if n == 0:
        return []
    if n <= k:
        return list(sequence_list)

    idxs = sorted(set(round((n - 1) * i / (k - 1)) for i in range(k)))
    return [sequence_list[i] for i in idxs]


def estimate_sec_per_query_for_reference(
    reference_key: str,
    cache_path: str,
    pilot_queries: Sequence[str],
) -> float:
    """
    Estimate sec/query for the reference using pilot queries.
    Perform only the search without writing CSV output.
    """
    ensure_reference_loaded(reference_key=reference_key, cache_path=cache_path)

    local_bundle_cache = {}
    start_t = time.perf_counter()

    for query_seq in pilot_queries:
        bundle = local_bundle_cache.get(query_seq)
        if bundle is None:
            bundle = build_query_payload_bundle_splitrank(
                query_seq=query_seq,
                rank_map=RANK_MAP,
                max_ext_total=MAX_EXT_TOTAL,
                max_mm=MAX_MISMATCH,
            )
            local_bundle_cache[query_seq] = bundle

        start_idx = bisect.bisect_left(G_REF_LENGTHS, bundle["min_ref_len_global"])
        current_best_group_rank = -1

        for _ref_seq, ref_bytes, _entry_names, rlen in G_REF_RECORDS[start_idx:]:
            if rlen < bundle["min_ref_len_global"]:
                continue

            payloads_both_exactlen = bundle["payloads_both_by_reflen"].get(rlen, None)

            best_rank_for_ref, hits = find_best_hits_multi_payloads_cy(
                ref_bytes=ref_bytes,
                payloads_both_exactlen=payloads_both_exactlen,
                payloads_left_by_corelen=bundle["payloads_left_by_corelen"],
                payloads_right_by_corelen=bundle["payloads_right_by_corelen"],
                anywhere_by_corelen=bundle["anywhere_by_corelen"],
                max_mm=bundle["max_mm"],
                current_best_group_rank=current_best_group_rank,
                prefix_check_len=PREFIX_CHECK_LEN,
            )

            if best_rank_for_ref < 0 or not hits:
                continue

            if current_best_group_rank < 0 or best_rank_for_ref < current_best_group_rank:
                current_best_group_rank = best_rank_for_ref

    elapsed = time.perf_counter() - start_t
    return elapsed / max(1, len(pilot_queries))


# =========================================================
# Batch size settings
# =========================================================
def choose_batch_size_ref(
    n_query: int,
    n_workers: int,
    sec_per_query: float,
    target_task_sec: int,
) -> tuple[int, dict]:
    """
    Use the query count, worker count, and reference cost in sec/query to
    determine the batch size for each reference.

    Objectives:
      - prevent excessive task counts as the query count increases
      - prevent individual tasks from becoming too long
      - ensure enough tasks for N_WORKERS
    """
    if n_query <= 0:
        return 1, {
            "target_task_sec": target_task_sec,
            "pred_total_sec": 0.0,
            "effective_target_sec": 0.0,
            "tasks_by_time": 1,
            "desired_tasks": 1,
            "estimated_tasks": 1,
            "estimated_task_sec": 0.0,
            "estimated_waves": 0.0,
            "min_tasks_ref": 1,
            "max_tasks_ref": 1,
            "hard_batch_max": 1,
        }

    sec_per_query = max(sec_per_query, 1e-9)

    pred_total_sec = sec_per_query * n_query

    # Apply the safety factor to the target task seconds rather than directly to the batch size.
    effective_target_sec = max(1.0, target_task_sec * BATCH_SAFETY_FACTOR)

    # Number of tasks based on the time target.
    tasks_by_time = max(1, math.ceil(pred_total_sec / effective_target_sec))

    # Lower and upper task-count limits based on worker count.
    min_tasks_ref = max(1, int(n_workers * MIN_WAVES_PER_REFERENCE))
    max_tasks_ref = max(min_tasks_ref, int(n_workers * MAX_WAVES_PER_REFERENCE))

    # Prevent the number of tasks from exceeding the number of queries.
    min_tasks_ref = min(min_tasks_ref, n_query)
    max_tasks_ref = min(max_tasks_ref, n_query)

    # In principle, keep the task count within this range.
    desired_tasks = clamp(tasks_by_time, min_tasks_ref, max_tasks_ref)

    # Determine the batch size from the desired task count.
    batch_size_ref = math.ceil(n_query / desired_tasks)

    # Upper batch-size limit to prevent individual tasks from becoming too long.
    hard_batch_max = max(1, int(MAX_TASK_SEC_HARD / sec_per_query))

    # Use the smaller of the physical upper limit and the hard upper limit.
    upper_batch_limit = min(MAX_BATCH_SIZE, hard_batch_max)

    # Prioritize the hard limit.
    # If upper_batch_limit < MIN_BATCH_SIZE, prioritize the hard limit over MIN_BATCH_SIZE.
    if upper_batch_limit < MIN_BATCH_SIZE:
        batch_size_ref = upper_batch_limit
    else:
        batch_size_ref = clamp(batch_size_ref, MIN_BATCH_SIZE, upper_batch_limit)

    estimated_task_sec = batch_size_ref * sec_per_query
    estimated_tasks = math.ceil(n_query / batch_size_ref)
    estimated_waves = estimated_tasks / max(n_workers, 1)

    info = {
        "target_task_sec": target_task_sec,
        "pred_total_sec": pred_total_sec,
        "effective_target_sec": effective_target_sec,
        "tasks_by_time": tasks_by_time,
        "desired_tasks": desired_tasks,
        "estimated_tasks": estimated_tasks,
        "estimated_task_sec": estimated_task_sec,
        "estimated_waves": estimated_waves,
        "min_tasks_ref": min_tasks_ref,
        "max_tasks_ref": max_tasks_ref,
        "hard_batch_max": hard_batch_max,
    }

    return batch_size_ref, info

def build_batch_size_and_cost_by_reference(sequence_list, cache_info_by_reference):
    """
    Run a pilot for each reference and
    determine the batch size from the query count, worker count, and sec/query.

    Returns:
      batch_size_by_reference
      sec_per_query_by_reference
    """
    pilot_queries = pick_pilot_queries_evenly(sequence_list, PILOT_QUERY_COUNT)
    n_query = len(sequence_list)

    target_task_sec = choose_target_task_sec(
        n_query=n_query,
        n_workers=N_WORKERS,
    )

    print(
        f"[AUTO TARGET] n_query={n_query}, "
        f"N_WORKERS={N_WORKERS}, "
        f"target_task_sec={target_task_sec}, "
        f"effective_target_task_sec={target_task_sec * BATCH_SAFETY_FACTOR:.1f}",
        flush=True,
    )

    batch_size_by_reference = {}
    sec_per_query_by_reference = {}

    for ref in reference_list:
        ref_name = ref["reference_name"]
        cache_path = cache_info_by_reference[ref_name]["cache_path"]

        sec_per_query = estimate_sec_per_query_for_reference(
            reference_key=ref_name,
            cache_path=cache_path,
            pilot_queries=pilot_queries,
        )
        sec_per_query_by_reference[ref_name] = sec_per_query

        batch_size_ref, info = choose_batch_size_ref(
            n_query=n_query,
            n_workers=N_WORKERS,
            sec_per_query=sec_per_query,
            target_task_sec=target_task_sec,
        )

        batch_size_by_reference[ref_name] = batch_size_ref

        print(
            f"[{ref_name}] pilot sec/query={sec_per_query:.3f}, "
            f"target_task_sec={info['target_task_sec']}, "
            f"effective_target_sec={info['effective_target_sec']:.1f}, "
            f"batch_size_ref={batch_size_ref}, "
            f"estimated_tasks={info['estimated_tasks']}, "
            f"estimated_waves={info['estimated_waves']:.1f}, "
            f"estimated_task_sec={info['estimated_task_sec']:.1f}s, "
            f"tasks_by_time={info['tasks_by_time']}, "
            f"desired_tasks={info['desired_tasks']}",

            flush=True,
        )

    return batch_size_by_reference, sec_per_query_by_reference
# ---------------------------
# Reference loading helpers
# ---------------------------
def ensure_reference_loaded(reference_key: str, cache_path: str):
    global G_CURRENT_REFERENCE_KEY, G_REF_RECORDS, G_REF_LENGTHS
    if G_CURRENT_REFERENCE_KEY == reference_key:
        return
    G_REF_RECORDS, G_REF_LENGTHS = load_reference_cache(cache_path)
    G_CURRENT_REFERENCE_KEY = reference_key


def get_start_idx_cached(reference_key: str, min_ref_len: int) -> int:
    key = (reference_key, min_ref_len)
    cached = G_START_IDX_CACHE.get(key)
    if cached is not None:
        return cached
    idx = bisect.bisect_left(G_REF_LENGTHS, min_ref_len)
    G_START_IDX_CACHE[key] = idx
    return idx


def close_open_temp_writer():
    global G_OPEN_TEMP_REFERENCE_NAME, G_OPEN_TEMP_PATH, G_OPEN_TEMP_FILE, G_OPEN_TEMP_WRITER
    if G_OPEN_TEMP_FILE is not None:
        try:
            G_OPEN_TEMP_FILE.flush()
        finally:
            G_OPEN_TEMP_FILE.close()
    G_OPEN_TEMP_REFERENCE_NAME = None
    G_OPEN_TEMP_PATH = None
    G_OPEN_TEMP_FILE = None
    G_OPEN_TEMP_WRITER = None


def ensure_temp_writer_for_reference(reference_name: str, temp_dir: str):
    global G_OPEN_TEMP_REFERENCE_NAME, G_OPEN_TEMP_PATH, G_OPEN_TEMP_FILE, G_OPEN_TEMP_WRITER, G_TEMP_PART_INDEX_BY_REFERENCE
    if G_OPEN_TEMP_REFERENCE_NAME == reference_name and G_OPEN_TEMP_FILE is not None:
        return G_OPEN_TEMP_WRITER, G_OPEN_TEMP_PATH, False

    close_open_temp_writer()
    next_part_idx = G_TEMP_PART_INDEX_BY_REFERENCE.get(reference_name, 0) + 1
    G_TEMP_PART_INDEX_BY_REFERENCE[reference_name] = next_part_idx
    worker_pid = os.getpid()
    temp_csv = os.path.join(temp_dir, f"worker_{worker_pid}_part_{next_part_idx:04d}.csv")
    f = open(temp_csv, "w", newline="", encoding="utf-8-sig")
    writer = csv.writer(f)
    writer.writerow(OUTPUT_COLUMNS)
    G_OPEN_TEMP_REFERENCE_NAME = reference_name
    G_OPEN_TEMP_PATH = temp_csv
    G_OPEN_TEMP_FILE = f
    G_OPEN_TEMP_WRITER = writer
    return writer, temp_csv, True


def worker_init():
    global G_CURRENT_REFERENCE_KEY, G_REF_RECORDS, G_REF_LENGTHS, G_START_IDX_CACHE
    global G_OPEN_TEMP_REFERENCE_NAME, G_OPEN_TEMP_PATH, G_OPEN_TEMP_FILE, G_OPEN_TEMP_WRITER, G_TEMP_PART_INDEX_BY_REFERENCE
    global G_QUERY_BUNDLE_CACHE
    ensure_rank_map_initialized()
    G_CURRENT_REFERENCE_KEY = None
    G_REF_RECORDS = None
    G_REF_LENGTHS = None
    G_START_IDX_CACHE = {}
    G_OPEN_TEMP_REFERENCE_NAME = None
    G_OPEN_TEMP_PATH = None
    G_OPEN_TEMP_FILE = None
    G_OPEN_TEMP_WRITER = None
    G_TEMP_PART_INDEX_BY_REFERENCE = {}
    G_QUERY_BUNDLE_CACHE = {}
    atexit.register(close_open_temp_writer)


# ---------------------------
# Query finalization helpers
# ---------------------------
def finalize_best_hits_for_query(query_seq: str, best_hits_compact: List[tuple], bundle: dict, writer) -> int:
    row_count = 0
    if not best_hits_compact:
        writer.writerow(
            build_output_row_list_fast(
                sequence=query_seq,
                reference_entry_name="Unannotated",
                hit=None,
                total_override=UNANNOTATED_TOTAL_SCORE,
            )
        )
        return 1

    for payload_id, rlen, entry_names, fixed_hit in best_hits_compact:
        payload = bundle["payload_meta_by_id"][payload_id]
        for entry_name in entry_names:
            writer.writerow(
                build_output_row_from_fixed_hit(
                    sequence=query_seq,
                    reference_entry_name=entry_name,
                    payload=payload,
                    fixed_hit=fixed_hit,
                    rlen=rlen,
                )
            )
            row_count += 1
    return row_count


# ---------------------------
# Main worker
# ---------------------------
def annotate_reference_superbatch_to_temp_csv(args):
    (reference_name, reference_key, cache_path, batch_ids, query_batches, temp_dir) = args
    global G_OPEN_TEMP_FILE, G_QUERY_BUNDLE_CACHE
    task_start = time.perf_counter()

    ensure_reference_loaded(reference_key=reference_key, cache_path=cache_path)
    writer, temp_csv, is_new_temp_file = ensure_temp_writer_for_reference(reference_name=reference_name, temp_dir=temp_dir)
    row_count = 0

    for _batch_id, query_batch in zip(batch_ids, query_batches):
        for query_seq in query_batch:
            bundle = G_QUERY_BUNDLE_CACHE.get(query_seq)
            if bundle is None:
                bundle = build_query_payload_bundle_splitrank(
                    query_seq=query_seq,
                    rank_map=RANK_MAP,
                    max_ext_total=MAX_EXT_TOTAL,
                    max_mm=MAX_MISMATCH,
                )
                G_QUERY_BUNDLE_CACHE[query_seq] = bundle

            start_idx = get_start_idx_cached(reference_key, bundle["min_ref_len_global"])
            current_best_group_rank = -1
            best_hits_compact: List[tuple] = []

            for _ref_seq, ref_bytes, entry_names, rlen in G_REF_RECORDS[start_idx:]:
                if rlen < bundle["min_ref_len_global"]:
                    continue

                payloads_both_exactlen = bundle["payloads_both_by_reflen"].get(rlen, None)

                best_rank_for_ref, hits = find_best_hits_multi_payloads_cy(
                    ref_bytes=ref_bytes,
                    payloads_both_exactlen=payloads_both_exactlen,
                    payloads_left_by_corelen=bundle["payloads_left_by_corelen"],
                    payloads_right_by_corelen=bundle["payloads_right_by_corelen"],
                    anywhere_by_corelen=bundle["anywhere_by_corelen"],
                    max_mm=bundle["max_mm"],
                    current_best_group_rank=current_best_group_rank,
                    prefix_check_len=PREFIX_CHECK_LEN,
                )

                if best_rank_for_ref < 0 or not hits:
                    continue

                if current_best_group_rank < 0 or best_rank_for_ref < current_best_group_rank:
                    current_best_group_rank = best_rank_for_ref
                    best_hits_compact = []

                if best_rank_for_ref != current_best_group_rank:
                    continue

                for fixed_hit in hits:
                    payload_id = fixed_hit[0]
                    best_hits_compact.append((payload_id, rlen, entry_names, fixed_hit))

            row_count += finalize_best_hits_for_query(query_seq=query_seq, best_hits_compact=best_hits_compact, bundle=bundle, writer=writer)

    if G_OPEN_TEMP_FILE is not None:
        G_OPEN_TEMP_FILE.flush()

    task_elapsed_sec = time.perf_counter() - task_start
    return {
        "reference_name": reference_name,
        "batch_ids": batch_ids,
        "temp_csv": temp_csv,
        "row_count": row_count,
        "task_elapsed_sec": task_elapsed_sec,
        "is_new_temp_file": is_new_temp_file,
    }


# ---------------------------
# Final merge (query-only)
# ---------------------------
def merge_temp_csvs_to_minimum_only(temp_csv_files, minimum_csv_path):
    header = None
    best_rows_by_key = {}
    for temp_csv in temp_csv_files:
        with open(temp_csv, "r", encoding="utf-8-sig") as fin:
            reader = csv.reader(fin)
            temp_header = next(reader, None)
            if temp_header is None:
                continue
            if header is None:
                header = temp_header
                idx_query = header.index("query_sequence")
                idx_total = header.index("total_mismatch_and_extension_count")
            elif temp_header != header:
                raise ValueError(f"Header mismatch detected in temp file: {temp_csv}")

            for row in reader:
                key = row[idx_query]
                try:
                    total = int(row[idx_total])
                except ValueError as exc:
                    raise ValueError(f"Invalid total_mismatch_and_extension_count value: {row[idx_total]!r}") from exc

                if key not in best_rows_by_key:
                    best_rows_by_key[key] = {"min_total": total, "rows": [row]}
                else:
                    current_min = best_rows_by_key[key]["min_total"]
                    if total < current_min:
                        best_rows_by_key[key] = {"min_total": total, "rows": [row]}
                    elif total == current_min:
                        best_rows_by_key[key]["rows"].append(row)

    if header is None:
        raise ValueError("No header found. temp_csv_files may be empty.")

    with open(minimum_csv_path, "w", newline="", encoding="utf-8-sig") as fout:
        writer = csv.writer(fout)
        writer.writerow(header)
        for key in sorted(best_rows_by_key.keys()):
            for row in best_rows_by_key[key]["rows"]:
                writer.writerow(row)


# ---------------------------
# heavy-first scheduler
# ---------------------------
def build_task_specs_by_reference(
    sequence_list,
    temp_dir_by_reference,
    cache_info_by_reference,
    batch_size_by_reference,
    sec_per_query_by_reference,
):
    """
    Create query batches for each reference and retain task specifications by reference.
    Execution order is not determined here.
    """
    total_batches_by_reference = {}
    total_tasks_by_reference = {}
    task_specs_by_reference = {}

    for ref in reference_list:
        ref_name = ref["reference_name"]
        batch_size_ref = batch_size_by_reference[ref_name]
        sec_per_query = sec_per_query_by_reference[ref_name]

        ref_batches = list(chunk_list(sequence_list, batch_size_ref))
        total_batches_by_reference[ref_name] = len(ref_batches)

        grouped_batches = []
        for start in range(0, len(ref_batches), BATCHES_PER_TASK):
            end = min(start + BATCHES_PER_TASK, len(ref_batches))
            batch_ids = list(range(start, end))
            query_batches = ref_batches[start:end]
            grouped_batches.append((batch_ids, query_batches))

        total_tasks_by_reference[ref_name] = len(grouped_batches)

        cache_path = cache_info_by_reference[ref_name]["cache_path"]
        ref_temp_dir = temp_dir_by_reference[ref_name]
        task_specs = []

        for batch_ids, query_batches in grouped_batches:
            n_queries_in_task = sum(len(qb) for qb in query_batches)
            predicted_task_sec = sec_per_query * n_queries_in_task

            task_specs.append({
                "reference_name": ref_name,
                "reference_key": ref_name,
                "cache_path": cache_path,
                "batch_ids": batch_ids,
                "query_batches": query_batches,
                "temp_dir": ref_temp_dir,
                "predicted_task_sec": predicted_task_sec,
                "n_queries_in_task": n_queries_in_task,
            })

        # Within each reference, place heavier tasks first so the final small batch is not scheduled first.
        # When weights are nearly identical, preserve batch ID order.
        task_specs.sort(
            key=lambda x: (
                -x["predicted_task_sec"],
                x["batch_ids"][0],
            )
        )
        task_specs_by_reference[ref_name] = task_specs

    return task_specs_by_reference, total_batches_by_reference, total_tasks_by_reference


def compute_reference_interleave_weights(task_specs_by_reference):
    """
    Create integer weights for weighted interleave from the predicted total workload by reference.

    Example: when hsa has a larger predicted total workload, its weight increases,
    so more hsa tasks are submitted, while other references are also submitted early.
    """
    total_predicted_by_reference = {}
    for ref in reference_list:
        ref_name = ref["reference_name"]
        specs = task_specs_by_reference.get(ref_name, [])
        total_predicted_by_reference[ref_name] = sum(
            max(float(spec.get("predicted_task_sec", 0.0)), 1e-9)
            for spec in specs
        )

    positive_totals = [
        v for v in total_predicted_by_reference.values()
        if v > 0
    ]
    if not positive_totals:
        return {
            ref["reference_name"]: 1
            for ref in reference_list
        }, total_predicted_by_reference

    min_positive = min(positive_totals)
    weights = {}
    for ref in reference_list:
        ref_name = ref["reference_name"]
        total_pred = total_predicted_by_reference.get(ref_name, 0.0)
        if total_pred <= 0:
            weights[ref_name] = 0
            continue

        relative = total_pred / min_positive
        scaled = relative ** WEIGHTED_INTERLEAVE_WEIGHT_EXP
        weight = int(round(scaled))
        weight = clamp(
            weight,
            WEIGHTED_INTERLEAVE_MIN_WEIGHT,
            WEIGHTED_INTERLEAVE_MAX_WEIGHT,
        )
        weights[ref_name] = weight

    return weights, total_predicted_by_reference


def smooth_weighted_interleave_task_specs(
    task_specs_by_reference,
    weights,
    total_predicted_by_reference=None,
):
    """
    Order tasks using dominant-first plus smooth weighted round-robin.

    1) If a dominant reference exists, prioritize it in the first wave.
       This shortens the critical path for the highest-load reference, such as hsa.
    2) Mix the remaining tasks using smooth weighted round-robin.
       This reduces the tail effect in which only mmu/bta/cfa remain at the end.
    """
    reference_order = [ref["reference_name"] for ref in reference_list]
    order_index = {name: i for i, name in enumerate(reference_order)}

    queues = {
        ref_name: deque(specs)
        for ref_name, specs in task_specs_by_reference.items()
    }
    ordered_specs = []

    if total_predicted_by_reference is None:
        total_predicted_by_reference = {
            ref_name: sum(
                max(float(spec.get("predicted_task_sec", 0.0)), 0.0)
                for spec in specs
            )
            for ref_name, specs in task_specs_by_reference.items()
        }

    # First wave: prioritize the dominant reference.
    # Otherwise, leave scheduling to normal weighted interleave.
    if DOMINANT_FIRST_WAVE and total_predicted_by_reference:
        dominant_ref = max(
            reference_order,
            key=lambda ref_name: total_predicted_by_reference.get(ref_name, 0.0),
        )
        total_pred_all = sum(max(v, 0.0) for v in total_predicted_by_reference.values())
        dominant_fraction = (
            total_predicted_by_reference.get(dominant_ref, 0.0) / total_pred_all
            if total_pred_all > 0 else 0.0
        )
        if dominant_fraction >= DOMINANT_FIRST_MIN_FRACTION:
            n_front = min(
                int(DOMINANT_FIRST_MAX_SLOTS),
                len(queues.get(dominant_ref, ())),
            )
            for _ in range(n_front):
                ordered_specs.append(queues[dominant_ref].popleft())

    # Remaining tasks: smooth weighted round-robin
    current = {ref_name: 0 for ref_name in reference_order}

    while True:
        active_refs = [
            ref_name for ref_name in reference_order
            if len(queues.get(ref_name, ())) > 0
        ]
        if not active_refs:
            break

        total_weight = sum(max(int(weights.get(ref_name, 0)), 0) for ref_name in active_refs)
        if total_weight <= 0:
            selected_ref = active_refs[0]
        else:
            for ref_name in active_refs:
                current[ref_name] += max(int(weights.get(ref_name, 0)), 0)

            selected_ref = max(
                active_refs,
                key=lambda ref_name: (
                    current[ref_name],
                    queues[ref_name][0]["predicted_task_sec"],
                    -order_index[ref_name],
                ),
            )
            current[selected_ref] -= total_weight

        ordered_specs.append(queues[selected_ref].popleft())

    return ordered_specs


def heavy_first_order_task_specs(task_specs_by_reference):
    task_specs = []
    for ref in reference_list:
        ref_name = ref["reference_name"]
        task_specs.extend(task_specs_by_reference.get(ref_name, []))

    task_specs.sort(
        key=lambda x: (
            -x["predicted_task_sec"],
            x["reference_name"],
            x["batch_ids"][0],
        )
    )
    return task_specs


def build_ordered_superbatch_tasks(
    sequence_list,
    temp_dir_by_reference,
    cache_info_by_reference,
    batch_size_by_reference,
    sec_per_query_by_reference,
):
    """
    Create task specifications and determine execution order using the specified task-ordering strategy.
    Recommended: TASK_ORDER_STRATEGY = "weighted_interleave"
    """
    (
        task_specs_by_reference,
        total_batches_by_reference,
        total_tasks_by_reference,
    ) = build_task_specs_by_reference(
        sequence_list=sequence_list,
        temp_dir_by_reference=temp_dir_by_reference,
        cache_info_by_reference=cache_info_by_reference,
        batch_size_by_reference=batch_size_by_reference,
        sec_per_query_by_reference=sec_per_query_by_reference,
    )

    order_info_by_reference = {}

    if TASK_ORDER_STRATEGY == "weighted_interleave":
        weights, total_predicted_by_reference = compute_reference_interleave_weights(task_specs_by_reference)
        task_specs = smooth_weighted_interleave_task_specs(
            task_specs_by_reference=task_specs_by_reference,
            weights=weights,
            total_predicted_by_reference=total_predicted_by_reference,
        )
        for ref in reference_list:
            ref_name = ref["reference_name"]
            order_info_by_reference[ref_name] = {
                "weight": weights.get(ref_name, 0),
                "total_predicted_sec": total_predicted_by_reference.get(ref_name, 0.0),
            }
    elif TASK_ORDER_STRATEGY == "heavy_first":
        task_specs = heavy_first_order_task_specs(task_specs_by_reference)
        for ref in reference_list:
            ref_name = ref["reference_name"]
            order_info_by_reference[ref_name] = {
                "weight": None,
                "total_predicted_sec": sum(
                    max(float(spec.get("predicted_task_sec", 0.0)), 0.0)
                    for spec in task_specs_by_reference.get(ref_name, [])
                ),
            }
    else:
        raise ValueError(
            f"Unknown TASK_ORDER_STRATEGY: {TASK_ORDER_STRATEGY!r}. "
            "Use 'weighted_interleave' or 'heavy_first'."
        )

    tasks = []
    for spec in task_specs:
        tasks.append(
            (
                spec["reference_name"],
                spec["reference_key"],
                spec["cache_path"],
                spec["batch_ids"],
                spec["query_batches"],
                spec["temp_dir"],
            )
        )

    return (
        tasks,
        total_batches_by_reference,
        total_tasks_by_reference,
        task_specs,
        order_info_by_reference,
    )


if __name__ == "__main__":
    mp.freeze_support()
    runtime_args = parse_arguments()
    configure_runtime(runtime_args)

    print(f"QUERY_FASTA = {sequence_fasta}")
    print(f"REFERENCES_FILE = {Path(runtime_args.references_file).resolve()}")
    print(f"OUTPUT_DIR = {output_folder}")
    print(f"CORE_NAME = {core_name}")
    print(f"N_WORKERS = {N_WORKERS}")
    if N_WORKERS < 1:
        raise ValueError("N_WORKERS must be >= 1")
    if BATCH_SIZE < 1:
        raise ValueError("BATCH_SIZE must be >= 1")
    if BATCHES_PER_TASK < 1:
        raise ValueError("BATCHES_PER_TASK must be >= 1")
    if MAX_MM > 3:
        raise ValueError("This prototype currently expects MAX_MM <= 3")

    ensure_rank_map_initialized()
    validate_input_paths()

    print(f"BATCH_SIZE_FALLBACK = {BATCH_SIZE}")
    print(f"BATCHES_PER_TASK = {BATCHES_PER_TASK}")
    print(f"UNANNOTATED_TOTAL_SCORE = {UNANNOTATED_TOTAL_SCORE}")
    print(f"PREFIX_CHECK_LEN = {PREFIX_CHECK_LEN}")
    print(f"USE_DYNAMIC_BATCH_SIZE = {USE_DYNAMIC_BATCH_SIZE}")
    print(f"PILOT_QUERY_COUNT = {PILOT_QUERY_COUNT}")
    print(f"MIN_BATCH_SIZE = {MIN_BATCH_SIZE}")
    print(f"MAX_BATCH_SIZE = {MAX_BATCH_SIZE}")
    print(f"BATCH_SAFETY_FACTOR = {BATCH_SAFETY_FACTOR}")
    print(f"CACHE_DIR = {CACHE_DIR}")
    print(f"MIN_WAVES_PER_REFERENCE = {MIN_WAVES_PER_REFERENCE}")
    print(f"MAX_WAVES_PER_REFERENCE = {MAX_WAVES_PER_REFERENCE}")
    print(f"MAX_TASK_SEC_HARD = {MAX_TASK_SEC_HARD}")
    print(f"USE_AUTO_TARGET_TASK_SEC = {USE_AUTO_TARGET_TASK_SEC}")
    print(f"TARGET_TASK_SEC_FALLBACK = {TARGET_TASK_SEC_FALLBACK}")
    print(f"TARGET_TASK_SEC_BASE = {TARGET_TASK_SEC_BASE}")
    print(f"TARGET_TASK_SEC_QUERY_REF = {TARGET_TASK_SEC_QUERY_REF}")
    print(f"TARGET_TASK_SEC_BASE_WORKERS = {TARGET_TASK_SEC_BASE_WORKERS}")
    print(f"TARGET_TASK_SEC_MIN = {TARGET_TASK_SEC_MIN}")
    print(f"TARGET_TASK_SEC_MAX = {TARGET_TASK_SEC_MAX}")
    print(f"TARGET_TASK_SEC_ROUND_TO = {TARGET_TASK_SEC_ROUND_TO}")
    print(f"TARGET_TASK_SEC_QUERY_EXP = {TARGET_TASK_SEC_QUERY_EXP}")
    print(f"TARGET_TASK_SEC_WORKER_EXP = {TARGET_TASK_SEC_WORKER_EXP}")
    print(f"TASK_ORDER_STRATEGY = {TASK_ORDER_STRATEGY}")
    print(f"WEIGHTED_INTERLEAVE_WEIGHT_EXP = {WEIGHTED_INTERLEAVE_WEIGHT_EXP}")
    print(f"WEIGHTED_INTERLEAVE_MIN_WEIGHT = {WEIGHTED_INTERLEAVE_MIN_WEIGHT}")
    print(f"WEIGHTED_INTERLEAVE_MAX_WEIGHT = {WEIGHTED_INTERLEAVE_MAX_WEIGHT}")
    print(f"DOMINANT_FIRST_WAVE = {DOMINANT_FIRST_WAVE}")
    print(f"DOMINANT_FIRST_MIN_FRACTION = {DOMINANT_FIRST_MIN_FRACTION}")
    print(f"DOMINANT_FIRST_MAX_SLOTS = {DOMINANT_FIRST_MAX_SLOTS}")


    sequence_list = read_fasta_sequences(sequence_fasta)

    if len(sequence_list) == 0:
        raise ValueError(f"No query sequences found in FASTA: {sequence_fasta}")

    print(f"Total unique sequences in fasta: {len(sequence_list)}")

    cache_info_by_reference = prepare_all_reference_caches(reference_list)
    
    if USE_DYNAMIC_BATCH_SIZE:
        batch_size_by_reference, sec_per_query_by_reference = build_batch_size_and_cost_by_reference(
            sequence_list=sequence_list,
            cache_info_by_reference=cache_info_by_reference,
        )
    else:
        batch_size_by_reference = {
            ref["reference_name"]: BATCH_SIZE for ref in reference_list
        }
        sec_per_query_by_reference = {
            ref["reference_name"]: 0.0 for ref in reference_list
        }


    run_stamp = time.strftime("%Y%m%d_%H%M%S")
    temp_root = os.path.join(tempfile.gettempdir(), f"{core_name}_tmp_{run_stamp}_{os.getpid()}")
    os.makedirs(temp_root, exist_ok=True)

    temp_dir_by_reference = {}
    for ref in reference_list:
        ref_name = ref["reference_name"]
        ref_temp_dir = os.path.join(temp_root, ref_name)
        os.makedirs(ref_temp_dir, exist_ok=True)
        temp_dir_by_reference[ref_name] = ref_temp_dir

    (
        tasks,
        total_batches_by_reference,
        total_tasks_by_reference,
        task_specs,
        order_info_by_reference,
    ) = build_ordered_superbatch_tasks(
        sequence_list=sequence_list,
        temp_dir_by_reference=temp_dir_by_reference,
        cache_info_by_reference=cache_info_by_reference,
        batch_size_by_reference=batch_size_by_reference,
        sec_per_query_by_reference=sec_per_query_by_reference,
    )

    for ref in reference_list:
        ref_name = ref["reference_name"]
        print(
            f"[{ref_name}] batch_size_ref={batch_size_by_reference[ref_name]} | "
            f"total query batches: {total_batches_by_reference[ref_name]} | "
            f"total tasks: {total_tasks_by_reference[ref_name]}",
            flush=True,
        )
    print("\n[TASK ORDER INFO]", flush=True)
    for ref in reference_list:
        ref_name = ref["reference_name"]
        info = order_info_by_reference.get(ref_name, {})
        weight = info.get("weight", None)
        total_pred = info.get("total_predicted_sec", 0.0)
        if weight is None:
            print(
                f"[{ref_name}] total_predicted_sec={total_pred:.1f}",
                flush=True,
            )
        else:
            print(
                f"[{ref_name}] interleave_weight={weight} | total_predicted_sec={total_pred:.1f}",
                flush=True,
            )

    print(f"\n[{TASK_ORDER_STRATEGY.upper()} TASK ORDER: top {TASK_ORDER_PREVIEW_N}]", flush=True)
    for i, spec in enumerate(task_specs[:TASK_ORDER_PREVIEW_N], start=1):
        print(
            f"{i:02d}. {spec['reference_name']} "
            f"batch_ids={spec['batch_ids'][0]}-{spec['batch_ids'][-1]} "
            f"n_queries={spec['n_queries_in_task']} "
            f"predicted_task_sec={spec['predicted_task_sec']:.1f}",
            flush=True,
        )


    print("\n[GLOBAL POOL START] " + ", ".join([ref["reference_name"] for ref in reference_list]), flush=True)

    global_start = time.perf_counter()
    temp_csv_paths_by_reference = {ref["reference_name"]: set() for ref in reference_list}
    total_rows_by_reference = {ref["reference_name"]: 0 for ref in reference_list}
    sum_task_sec_by_reference = {ref["reference_name"]: 0.0 for ref in reference_list}
    done_tasks_by_reference = {ref["reference_name"]: 0 for ref in reference_list}
    done_batches_by_reference = {ref["reference_name"]: 0 for ref in reference_list}

    with mp.Pool(processes=N_WORKERS, initializer=worker_init) as pool:
        for result in pool.imap_unordered(annotate_reference_superbatch_to_temp_csv, tasks, chunksize=1):
            reference_name = result["reference_name"]
            batch_ids = result["batch_ids"]
            temp_csv = result["temp_csv"]
            row_count = result["row_count"]
            task_sec = result["task_elapsed_sec"]
            is_new_temp_file = result["is_new_temp_file"]

            temp_csv_paths_by_reference[reference_name].add(temp_csv)
            total_rows_by_reference[reference_name] += row_count
            done_tasks_by_reference[reference_name] += 1
            done_batches_by_reference[reference_name] += len(batch_ids)
            sum_task_sec_by_reference[reference_name] += task_sec

            done_task_i = done_tasks_by_reference[reference_name]
            total_task_i = total_tasks_by_reference[reference_name]
            done_batch_i = done_batches_by_reference[reference_name]
            total_batch_i = total_batches_by_reference[reference_name]
            avg_task_sec = sum_task_sec_by_reference[reference_name] / done_task_i
            remaining_tasks = total_task_i - done_task_i
            est_remaining_sec = remaining_tasks * avg_task_sec
            total_elapsed = time.perf_counter() - global_start
            temp_part_note = "new_temp_part" if is_new_temp_file else "append_existing_temp_part"

            print(
                f"[{reference_name}] finished task {done_task_i}/{total_task_i} "
                f"(batches {batch_ids[0]}-{batch_ids[-1]}, done_batches={done_batch_i}/{total_batch_i}, "
                f"task_sec={task_sec:.1f}s, temp_mode={temp_part_note}, global_elapsed={total_elapsed / 60:.1f} min, "
                f"avg_task={avg_task_sec:.1f}s, est_remaining_if_single_worker={est_remaining_sec / 60:.1f} min)",
                flush=True,
            )

    print("\n[GLOBAL POOL DONE]", flush=True)

    close_open_temp_writer()
    for ref in reference_list:
        ref_name = ref["reference_name"]
        temp_csv_files = sorted(temp_csv_paths_by_reference[ref_name])
        total_rows = total_rows_by_reference[ref_name]
        minimum_csv = os.path.join(output_folder, f"{core_name}_{ref_name}_minimumMiss.csv")
        merge_temp_csvs_to_minimum_only(temp_csv_files=temp_csv_files, minimum_csv_path=minimum_csv)
        print(f"[DONE] {ref_name}", flush=True)
        print(f"[{ref_name}] output (minimum only): {minimum_csv}", flush=True)
        print(f"[{ref_name}] temp CSV parts used for merge: {len(temp_csv_files)}", flush=True)
        print(f"[{ref_name}] rows written from temp tasks: {total_rows}", flush=True)

    total_elapsed = time.perf_counter() - global_start
    print(f"\nAll annotation jobs finished. total elapsed: {total_elapsed / 60:.2f} min")

    if not KEEP_TEMP_FILES:
        shutil.rmtree(temp_root, ignore_errors=True)
