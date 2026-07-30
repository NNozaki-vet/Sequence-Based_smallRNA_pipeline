# Annotation Method

## Scope

The workflow compares each unique query sequence independently with each reference FASTA collection. It does not perform genome alignment, probabilistic mapping, read-count normalization, expression analysis, or annotation priority resolution across different reference collections.

RNA sequences are converted to DNA notation before comparison:

```text
U -> T
```

## Query and reference preprocessing

Query FASTA records are:

1. concatenated when line-wrapped
2. converted to uppercase
3. stripped of spaces and line breaks
4. converted from `U` to `T`
5. deduplicated while preserving first-occurrence order

Reference FASTA records are processed similarly. Identical reference sequences are consolidated internally, while all associated FASTA entry names are retained.

## Candidate query splits

For each query, the workflow enumerates terminal query extensions satisfying:

```text
query_5prime_extension_length + query_3prime_extension_length <= 3
```

The remaining query core is compared with a reference sequence. The four placement modes are:

| Query extensions | Placement rule |
|---|---|
| 5′ = 0 and 3′ = 0 | Full query is searched at every valid reference position |
| 5′ > 0 and 3′ = 0 | Query core is fixed at the reference 5′ end |
| 5′ = 0 and 3′ > 0 | Query core is fixed at the reference 3′ end |
| 5′ > 0 and 3′ > 0 | Query core must span the complete reference sequence |

Reference bases outside the aligned query core are reported as reference 5′ or 3′ unmatched lengths.

## Mismatches

The compiled core allows at most three internal mismatches:

```text
MAX_MM = 3
```

Mismatch positions are reported relative to the aligned reference segment and use 1-based coordinates. Base changes are reported in the direction:

```text
reference_base>query_base
```

Insertions and deletions inside the aligned core are not modeled. Terminal differences are represented only by query extension lengths and reference unmatched lengths.

## Ranking

Candidate hits are ranked lexicographically by:

1. lower mismatch count
2. lower total query extension length
3. lower 5′ query extension length

Therefore, mismatch count has priority over total terminal extension. For equal mismatch count and total extension, a candidate with less 5′ extension receives the better rank.

All hits with the best rank for a query within a reference collection are retained.

## Per-reference output

Each reference definition is processed independently and produces a separate CSV file. The workflow does not select a single annotation across different reference files. Any biological annotation hierarchy, such as mature miRNA before precursor or species-specific before cross-species annotation, must be applied in a downstream, explicitly documented step.

For a query with no accepted hit in a reference collection:

```text
reference_entry_name = Unannotated
total_mismatch_and_extension_count = 10
```

The value `10` is a sentinel and is not an observed alignment distance.

## Parallel execution and caching

Parallel scheduling, dynamic batch sizing, and task order are performance mechanisms. They are not intended to change the accepted hit set.

Reference caches contain preprocessed reference records and are reused when the stored absolute path, file size, modification time, and cache format version match. For formal reproduction, use a new output directory or delete `reference_cache` before rerunning.
