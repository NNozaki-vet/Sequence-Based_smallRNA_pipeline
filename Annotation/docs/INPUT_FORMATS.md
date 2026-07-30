# Input Formats

## Query FASTA

`--query-fasta` must be an uncompressed, valid FASTA file.

Example:

```fasta
>sequence_1
TGAGGTAGTAGGTTGTATAGTT
>sequence_2
ACCCGTAGATCCGAACTTGTG
```

Requirements:

- UTF-8 text without a byte-order mark is recommended.
- Every sequence must follow a header beginning with `>`.
- Line-wrapped sequences are supported.
- Blank lines are ignored.
- FASTQ is not supported.
- Gzip-compressed FASTA is not supported.
- Sequence headers are not included in the output.
- Duplicate sequences are collapsed.
- Sequence characters are not restricted to `A`, `C`, `G`, `T`, and `N`; therefore, users must validate the alphabet before analysis.
- Counts encoded in FASTA headers are not parsed.

Normalization performed by the script:

```text
uppercase
remove spaces and line breaks
U -> T
```

For reproducibility, provide one FASTA record per unique query sequence and archive the exact input file.

## Reference definition CSV

`--references-file` must contain exactly named columns:

```csv
reference_name,file_path
```

Example:

```csv
reference_name,file_path
cfa_mature_miR,references/cfa-mature_miR.fa
cfa_pre_mir,references/cfa-pre.fa
```

Rules:

- Column names are case-sensitive.
- Additional columns are ignored.
- Empty rows are ignored.
- `reference_name` must be nonempty and unique.
- `file_path` must be nonempty.
- Relative paths are resolved from the CSV directory.
- Reference order in the CSV determines reporting and scheduling order, but does not create cross-reference annotation priority.
- Avoid Windows filename characters in `reference_name`: `< > : " / \ | ? *`.

## Reference FASTA

Each reference file must be an uncompressed FASTA file.

Example:

```fasta
>cfa-miR-example-5p
TGAGGTAGTAGGTTGTATAGTT
```

Requirements:

- UTF-8 text is required; UTF-8 with a byte-order mark is accepted.
- Headers must be nonempty.
- Sequences must be nonempty.
- Line-wrapped sequences are supported.
- Gzip-compressed FASTA is not supported.
- The complete text after `>` is used as `reference_entry_name`.
- Identical reference sequences may have multiple headers; all names are retained.
- Users must validate sequence alphabet and biological orientation before analysis.

## Output directory and core name

`--output-dir` is created when absent.

`--core-name` is optional. When omitted, it is derived from the query FASTA filename. It must not be empty and must not contain:

```text
< > : " / \ | ? *
```

Use a stable, dataset-specific core name and record it with the final command.
