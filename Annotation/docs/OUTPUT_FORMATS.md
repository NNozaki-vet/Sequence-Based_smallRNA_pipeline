# Output Formats

## Annotation files

One UTF-8 CSV with a byte-order mark is written for each reference:

```text
<core_name>_<reference_name>_minimumMiss.csv
```

Rows are grouped by query sequence. Multiple rows are retained when multiple reference entries have the best rank.

## Columns

| Column | Meaning |
|---|---|
| `query_sequence` | Normalized query sequence in DNA notation |
| `reference_entry_name` | Complete reference FASTA header without `>`, or `Unannotated` |
| `total_mismatch_and_extension_count` | `mismatch_count + query_5prime_extension_length + query_3prime_extension_length`; `10` is used as the unannotated sentinel |
| `mismatch_count` | Number of internal substitutions in the aligned core |
| `query_5prime_extension_length` | Number of query bases excluded from the 5′ end before core comparison |
| `query_5prime_extension_sequence` | Excluded 5′ query bases |
| `reference_5prime_unmatched_length` | Number of reference bases before the aligned core |
| `reference_mismatch_position1_1based` | First mismatch position within the aligned reference segment |
| `reference_to_query_mismatch1` | First substitution as `reference>query` |
| `reference_mismatch_position2_1based` | Second mismatch position |
| `reference_to_query_mismatch2` | Second substitution |
| `reference_mismatch_position3_1based` | Third mismatch position |
| `reference_to_query_mismatch3` | Third substitution |
| `query_3prime_extension_length` | Number of query bases excluded from the 3′ end before core comparison |
| `query_3prime_extension_sequence` | Excluded 3′ query bases |
| `reference_3prime_unmatched_length` | Number of reference bases after the aligned core |

Empty mismatch fields mean that the corresponding mismatch does not exist.

## Important interpretation rules

`total_mismatch_and_extension_count` is a reporting value, not the complete ranking rule. Best-hit selection first prioritizes mismatch count, then total query extension, then 5′ query extension.

`reference_5prime_unmatched_length` and `reference_3prime_unmatched_length` describe flanking reference sequence. They are not included in `total_mismatch_and_extension_count`.

Files from separate references must not be combined by choosing the smallest reporting value unless a downstream annotation-priority rule has been defined and documented.

## Reference cache

The following directory is created:

```text
<output_dir>\reference_cache\
```

Cache files are implementation artifacts, not primary results. Do not cite or archive them as substitutes for reference FASTA files.

## Temporary files

Worker CSV files are created under the operating-system temporary directory and are deleted when the run completes because:

```text
KEEP_TEMP_FILES = False
```

Retain the final CSV files and run log, not the temporary worker files.
