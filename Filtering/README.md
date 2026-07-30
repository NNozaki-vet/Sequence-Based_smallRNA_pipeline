# Exact-Sequence Filtering for Small RNA FASTA Files

This directory contains `filter_sequences.py`, a Windows-oriented, open-source script for filtering small RNA reads at the **exact nucleotide-sequence level**.

The script was generalized from the filtering procedure used for healthy-sample analysis in the associated study. It accepts user-defined sample groups and paths while preserving the original filtering criteria and computational workflow.

No proprietary bioinformatics software and no third-party Python packages are required for this filtering step.

---

## 1. Overview

`filter_sequences.py` performs the following operations:

1. Reads adapter-trimmed FASTA files.
2. Treats each FASTA record as one sequencing read.
3. Counts identical nucleotide sequences within each sample.
4. Calculates counts per million reads (CPM) separately for each sample.
5. Applies the following sample-level criteria:

   - raw count >= 2
   - CPM >= 1

6. Within each biological group, requires the sample-level criteria to be met in at least:

   ```text
   max(2, ceiling(number of samples in the group / 2))
   ```

7. Retains the union of all sequences that pass in at least one group.
8. Writes a raw-count matrix, filtered FASTA file, group-pass summary, sample metadata, and group-specific filtering thresholds.

The script preserves exact sequences rather than collapsing reads into predefined miRNA names. It is therefore suitable for downstream sequence-based annotation of mature miRNAs, isomiRs, miRNA-offset RNAs, precursor-derived fragments, antisense sequences, and other small RNA species.

---

## 2. Important input concept

The script uses the following structure:

```text
One FASTA record = one sequencing read
One FASTA file   = one biological sample
One folder       = one biological group
One CSV row      = one biological group
```

Examples of biological groups include:

- tissue types, such as `Liver`, `Kidney`, or `Brain`
- sample sources, such as `Plasma` or `Urine`
- experimental groups, such as `Control` or `Disease`
- sequencing fractions, such as `Protein_fraction` or `EV_enriched`

The filtering procedure was originally used for healthy tissue and body-fluid groups, but the generalized script can process any user-defined groups.

---

## 3. Requirements

### Operating system

The instructions below are written for:

- Windows 10
- Windows 11
- Anaconda Prompt, Miniforge Prompt, Command Prompt, or PowerShell

The script uses only the Python standard library and may also run on Linux or macOS, although the examples in this README focus on Windows.

### Python

Recommended:

```text
Python 3.9 or later
```

No additional Python packages are required.

Check the installed Python version:

```bat
python --version
```

An optional isolated conda environment can be created as follows:

```bat
conda create -n sequence_filter python=3.11
conda activate sequence_filter
```

---

## 4. Input FASTA requirements

### 4.1 Read-level FASTA is required

Each FASTA record must represent one read.

Correct example:

```fasta
>read_000001
TGAGGTAGTAGGTTGTATAGTT
>read_000002
TGAGGTAGTAGGTTGTATAGTT
>read_000003
TGTAAACATCCCCGACTGGAAG
```

In this example:

- `TGAGGTAGTAGGTTGTATAGTT` is counted twice.
- `TGTAAACATCCCCGACTGGAAG` is counted once.

### 4.2 Collapsed FASTA is not supported

The script does not read abundance information from FASTA headers.

The following collapsed format is not appropriate:

```fasta
>sequence_1 count=250
TGAGGTAGTAGGTTGTATAGTT
```

This record would be counted as one read, not 250 reads.

If the input data are collapsed, they must first be expanded to read-level FASTA or converted using a separate count-aware workflow.

### 4.3 Preprocessing must be completed before filtering

This script does not perform:

- adapter trimming
- quality filtering
- read-length filtering
- FASTQ-to-FASTA conversion
- removal of ambiguous reads
- removal of contaminating RNA classes
- alignment or annotation

The input FASTA files should therefore already represent the reads that are intended to contribute to the sample library size.

### 4.4 Supported file extensions

The following extensions are detected:

```text
.fa
.fasta
.fna
.fa.gz
.fasta.gz
.fna.gz
```

Files are searched only in the top level of each group folder. Subdirectories are not searched recursively.

### 4.5 Sequence handling

During reading, the script:

- converts sequences to uppercase
- joins line-wrapped FASTA sequences
- converts `U` to `T`
- ignores blank lines
- ignores FASTA header contents for counting
- skips empty sequence records

The script does not restrict the nucleotide alphabet. Any nonempty sequence text is counted after uppercasing and `U`-to-`T` conversion. Input quality control should therefore be completed before running this script.

---

## 5. Recommended directory structure

For example, assume that healthy liver, kidney, and brain samples will be analyzed:

```text
D:\smallRNA_example\
│
├─ input_fasta\
│  ├─ Liver\
│  │  ├─ Liver_01.fa
│  │  ├─ Liver_02.fa
│  │  └─ Liver_03.fa
│  │
│  ├─ Kidney\
│  │  ├─ Kidney_01.fa
│  │  ├─ Kidney_02.fa
│  │  ├─ Kidney_03.fa
│  │  └─ Kidney_04.fa
│  │
│  └─ Brain\
│     ├─ Brain_01.fa
│     ├─ Brain_02.fa
│     └─ Brain_03.fa
│
├─ group_directories.csv
└─ filter_sequences.py
```

Each FASTA file is treated as one sample.

The sample name is derived from the FASTA filename. Supported filename suffixes are removed automatically. The following trailing tags are also removed when present:

```text
.sra_trimmed
_sra_trimmed
.trimmed
_trimmed
```

For example:

```text
SRR000001.sra_trimmed.fa
```

becomes:

```text
SRR000001
```

---

## 6. Create the group configuration CSV

Create a CSV file with exactly these required column names:

```csv
group,folder
```

Example using absolute Windows paths:

```csv
group,folder
Liver,D:\smallRNA_example\input_fasta\Liver
Kidney,D:\smallRNA_example\input_fasta\Kidney
Brain,D:\smallRNA_example\input_fasta\Brain
```

The order of rows in the CSV determines:

- the order in which groups are processed
- the order of groups in `pass_groups`
- the group order in `group_filter_thresholds.csv`

### Paths containing commas

Enclose a path in double quotation marks when it contains a comma:

```csv
group,folder
Liver,"D:\smallRNA, example\input_fasta\Liver"
```

### Relative paths

Relative folder paths are interpreted relative to the directory containing the groups CSV.

Example:

```csv
group,folder
Liver,input_fasta\Liver
Kidney,input_fasta\Kidney
Brain,input_fasta\Brain
```

### Group-name requirements

Group names must:

- not be empty
- be unique within the CSV

Simple names containing letters, numbers, hyphens, or underscores are recommended:

```text
Liver
Oral_mucosa
Healthy_plasma
EV-enriched
```

Duplicate group names cause the script to stop with an error.

---

## 7. Run the script

### Basic command

```bat
python filter_sequences.py ^
  --groups-file "D:\smallRNA_example\group_directories.csv" ^
  --output-dir "D:\smallRNA_example\filtering_results"
```

### Specify the number of worker processes

```bat
python filter_sequences.py ^
  --groups-file "D:\smallRNA_example\group_directories.csv" ^
  --output-dir "D:\smallRNA_example\filtering_results" ^
  --threads 8
```

### Specify an output filename prefix

```bat
python filter_sequences.py ^
  --groups-file "D:\smallRNA_example\group_directories.csv" ^
  --output-dir "D:\smallRNA_example\filtering_results" ^
  --output-prefix "healthy_filtered_sequences" ^
  --threads 8
```

### One-line command

```bat
python filter_sequences.py --groups-file "D:\smallRNA_example\group_directories.csv" --output-dir "D:\smallRNA_example\filtering_results" --output-prefix "healthy_filtered_sequences" --threads 8
```

### Run from a conda environment

```bat
conda activate sequence_filter
cd /d "D:\smallRNA_example"
python filter_sequences.py --groups-file "D:\smallRNA_example\group_directories.csv" --output-dir "D:\smallRNA_example\filtering_results" --output-prefix "healthy_filtered_sequences" --threads 8
```

---

## 8. Command-line arguments

| Argument | Required | Default | Description |
|---|---:|---:|---|
| `--groups-file` | Yes | None | CSV containing the columns `group` and `folder` |
| `--output-dir` | Yes | None | Directory in which output files are written |
| `--output-prefix` | No | `filtered_sequences` | Prefix for the main count table, FASTA, and pass-group summary |
| `--threads` | No | `8` | Number of parallel worker processes used for per-sample counting |

Use a positive integer for `--threads`.

To display the built-in help:

```bat
python filter_sequences.py --help
```

---

## 9. Filtering algorithm

### 9.1 Per-sample exact-sequence counting

For every sample, the script counts the number of FASTA records with each exact nucleotide sequence.

For sequence \(i\) in sample \(j\):

```text
raw_count(i,j) = number of FASTA records exactly matching sequence i
```

The total number of reads in sample \(j\) is:

```text
library_size(j) = total number of nonempty FASTA sequence records
```

### 9.2 CPM calculation

CPM is calculated independently for each sample:

```text
CPM(i,j) = raw_count(i,j) / library_size(j) × 1,000,000
```

The denominator is the total number of input FASTA records in that sample, not the number of retained reads after filtering.

### 9.3 Sample-level pass criterion

A sequence passes in a sample only when both conditions are met:

```text
raw count >= 2
CPM >= 1
```

Both comparisons are inclusive.

### 9.4 Group-level required sample number

For a group containing \(n\) samples:

```text
required_samples = max(2, ceiling(n / 2))
```

Examples:

| Samples in group (`n`) | Required passing samples |
|---:|---:|
| 1 | 2 |
| 2 | 2 |
| 3 | 2 |
| 4 | 2 |
| 5 | 3 |
| 6 | 3 |
| 7 | 4 |
| 8 | 4 |
| 9 | 5 |
| 10 | 5 |

A group containing only one sample cannot retain any sequence because the required number is two. At least two biological samples per group are therefore recommended.

### 9.5 Final retained sequence set

A sequence is retained when it passes the group-level criterion in at least one group.

In set notation:

```text
retained sequences = union of sequences passing in each group
```

For example, a sequence that passes only in `Liver` is retained even if it does not pass in `Kidney` or `Brain`.

This design preserves sequences that are reproducibly detected within a specific tissue or biological group rather than requiring expression across all groups.

---

## 10. Worked example

Suppose the `Liver` group contains three samples.

For one sequence:

| Sample | Raw count | Total reads | CPM | Sample pass |
|---|---:|---:|---:|---|
| Liver_01 | 4 | 2,000,000 | 2.0 | Yes |
| Liver_02 | 2 | 1,000,000 | 2.0 | Yes |
| Liver_03 | 1 | 1,000,000 | 1.0 | No |

For three samples:

```text
required_samples = max(2, ceiling(3 / 2)) = 2
```

The sequence passes in two of the three samples and is therefore retained for the `Liver` group.

A sequence with raw count 1 and CPM 5 does not pass because the raw-count criterion is not met.

A sequence with raw count 10 and CPM 0.8 does not pass because the CPM criterion is not met.

---

## 11. Output files

Assume the following option is used:

```text
--output-prefix healthy_filtered_sequences
```

The output directory will contain:

```text
filtering_results\
│
├─ healthy_filtered_sequences_count_table.csv
├─ healthy_filtered_sequences.fa
├─ healthy_filtered_sequences_pass_groups.csv
├─ sample_metadata.csv
├─ group_filter_thresholds.csv
└─ per_sample_sequence_counts\
```

The script recreates the per-sample counts from all input FASTA files each time it is run. Existing output files with the same names may be overwritten.

---

## 12. Main count table

Filename:

```text
<output-prefix>_count_table.csv
```

Example:

```text
healthy_filtered_sequences_count_table.csv
```

Columns:

| Column | Description |
|---|---|
| `sequence` | Exact nucleotide sequence after uppercasing and `U`-to-`T` conversion |
| `length` | Sequence length in nucleotides |
| `total_count` | Sum of raw counts across all included samples |
| `n_pass_groups` | Number of groups in which the sequence passed the group-level criterion |
| `pass_groups` | Semicolon-separated names of passing groups |
| sample columns | Raw count of the sequence in each sample |

Sample columns are ordered by:

1. group order in `group_directories.csv`
2. alphabetical FASTA filename order within each group

Sequences are ordered by:

1. descending `total_count`
2. ascending sequence length
3. alphabetical sequence order

### Duplicate sample names

If two input FASTA files produce the same sample name, the script appends the internal sample index:

```text
sample_name__idx0
sample_name__idx1
```

This prevents duplicate count-table column names.

---

## 13. Filtered FASTA

Filename:

```text
<output-prefix>.fa
```

Example:

```text
healthy_filtered_sequences.fa
```

Each retained exact sequence is written once.

Example record:

```fasta
>seq_000000001 length=22 total_count=1540 n_pass_groups=3 pass_groups=Liver;Kidney;Brain
TGAGGTAGTAGGTTGTATAGTT
```

Header fields:

| Field | Description |
|---|---|
| `seq_#########` | Sequential output identifier |
| `length` | Sequence length |
| `total_count` | Sum of raw counts across all samples |
| `n_pass_groups` | Number of passing groups |
| `pass_groups` | Semicolon-separated passing group names |

This FASTA can be used as input for downstream sequence-based annotation.

---

## 14. Pass-group summary

Filename:

```text
<output-prefix>_pass_groups.csv
```

Columns:

| Column | Description |
|---|---|
| `sequence` | Retained exact sequence |
| `length` | Sequence length |
| `total_count` | Total raw count across all samples |
| `n_pass_groups` | Number of groups in which the sequence passed |
| `pass_groups` | Semicolon-separated passing group names |

This file contains the filtering summary without the sample-level raw-count matrix.

---

## 15. Sample metadata

Filename:

```text
sample_metadata.csv
```

Columns:

| Column | Description |
|---|---|
| `sample_idx` | Internal zero-based sample index |
| `group` | Group assigned in the groups CSV |
| `sample_name` | Sample name derived from the FASTA filename |
| `fasta_path` | Input FASTA path |
| `count_tsv` | Per-sample sequence-count TSV path |
| `total_reads` | Number of nonempty FASTA sequence records |
| `unique_sequences` | Number of distinct exact sequences |
| `group_n_samples` | Number of samples successfully detected in the group |
| `group_required_samples` | Minimum number of passing samples required for that group |
| `low_depth_lt_1M` | `YES` when `total_reads < 1,000,000`, otherwise `NO` |

The `low_depth_lt_1M` column is informational only. Low-depth samples are not automatically removed.

Always inspect this file before downstream interpretation.

---

## 16. Group threshold summary

Filename:

```text
group_filter_thresholds.csv
```

Columns:

| Column | Description |
|---|---|
| `group` | Group name |
| `n_samples` | Number of FASTA samples found in that group |
| `required_samples` | Number of samples required to pass the sample-level criteria |

This file documents the exact group-specific thresholds used in the run.

---

## 17. Per-sample count files

Directory:

```text
per_sample_sequence_counts\
```

Each sample receives one tab-delimited count file:

```text
0000__Liver__Liver_01.seq_counts.tsv
```

Columns:

```text
sequence
count
```

These files contain all exact sequences detected in the corresponding sample before group-level filtering.

---

## 18. Healthy-tissue analysis example

To reproduce a healthy-tissue-style analysis, place each tissue in a separate group folder.

Example:

```csv
group,folder
Thymus,D:\data\Trimmed_fa\Thymus
Bone_Marrow,D:\data\Trimmed_fa\Bone_Marrow
Brain,D:\data\Trimmed_fa\Brain
Colon,D:\data\Trimmed_fa\Colon
Heart,D:\data\Trimmed_fa\Heart
Kidney,D:\data\Trimmed_fa\Kidney
Liver,D:\data\Trimmed_fa\Liver
Lung,D:\data\Trimmed_fa\Lung
Plasma,D:\data\Trimmed_fa\Plasma
```

Run:

```bat
python filter_sequences.py ^
  --groups-file "D:\data\healthy_groups.csv" ^
  --output-dir "D:\data\healthy_filtered_sequences" ^
  --output-prefix "healthy_filtered_sequences" ^
  --threads 8
```

A sequence passing in any individual healthy tissue or body-fluid group will be included in the final retained sequence set.

---

## 19. Using all healthy samples as one group

To require detection across the combined healthy cohort rather than within individual tissues, place all FASTA files in one folder or assign all files to one group folder.

Example CSV:

```csv
group,folder
Healthy,D:\data\all_healthy_samples
```

In that design, the required sample number is calculated from the total number of samples in the combined `Healthy` group.

This is methodologically different from tissue-wise filtering. Choose the grouping strategy before analysis and report it clearly.

---

## 20. Warnings and input validation

### Missing group folder

If a folder listed in the CSV does not exist, the script prints:

```text
[WARNING] Directory not found: ...
```

and continues with the remaining groups.

### No FASTA files in a group folder

If no supported FASTA files are detected, the script prints:

```text
[WARNING] No FASTA files in: ...
```

and continues.

Because skipped groups change the effective sample design, review all warnings and inspect:

```text
sample_metadata.csv
group_filter_thresholds.csv
```

before using the results.

### No samples found

The script stops when no FASTA samples are found in any valid group.

### Duplicate groups

Duplicate group names in the groups CSV cause the script to stop.

### Empty group or folder field

A partially empty CSV row causes the script to stop and reports the row number.

---

## 21. Troubleshooting

### `Groups file not found`

Confirm the full path passed to `--groups-file`.

```bat
dir "D:\smallRNA_example\group_directories.csv"
```

### `The groups file is missing required columns`

The header must contain exactly the required names:

```csv
group,folder
```

Capitalization matters.

### `No FASTA samples were found`

Check:

- whether the folders exist
- whether files are directly inside the folders
- whether extensions are supported
- whether the CSV paths are correct
- whether the files are not hidden inside subdirectories

### Unexpectedly low read counts

Confirm that the FASTA is read-level rather than collapsed.

Also confirm that adapter trimming and FASTQ-to-FASTA conversion did not merge identical reads.

### A one-sample group retains no sequences

This is expected.

For `n = 1`:

```text
max(2, ceiling(1 / 2)) = 2
```

A minimum of two samples is therefore required.

### Output contains `T` instead of `U`

This is expected. The script converts RNA-style `U` to DNA-style `T` for sequence consistency.

### Windows multiprocessing error

Run the script as a file from a terminal:

```bat
python filter_sequences.py ...
```

Do not paste the entire script interactively into the Python interpreter.

### High memory use

Memory use increases with:

- the number of unique sequences
- the number of samples
- the number of retained sequences

Reduce `--threads` if simultaneous per-sample counting uses too much memory:

```bat
--threads 4
```

---

## 22. Reproducibility recommendations

For every analysis, retain:

- the exact `filter_sequences.py` version
- the groups CSV
- the complete command used
- `sample_metadata.csv`
- `group_filter_thresholds.csv`
- the filtered FASTA
- the raw-count matrix
- preprocessing details for the input FASTA files

Record the following in the Methods section:

- read preprocessing criteria
- grouping strategy
- raw-count threshold
- CPM threshold
- group-level required sample rule
- whether sequences were represented using `T` or `U`

For this script, the fixed filtering criteria are:

```text
raw count >= 2
CPM >= 1
required samples per group = max(2, ceiling(n / 2))
```

---

## 23. Scope and limitations

This script is designed for transparent exact-sequence filtering. It is not a complete small RNA-seq analysis pipeline by itself.

It does not:

- identify miRNAs
- distinguish miRNAs from other small RNAs
- assign genomic coordinates
- perform differential expression analysis
- correct batch effects
- normalize count matrices for statistical modeling
- parse abundance values from collapsed FASTA headers
- recursively search nested directories
- automatically exclude low-depth samples
- automatically stop when one listed group is missing

The output count matrix contains raw counts and is suitable as an input for downstream statistical methods that require raw integer counts, provided that the biological design and library-size handling are appropriate.

---

## 24. Method summary for manuscripts

A concise description suitable for adaptation in a Methods section is:

> Adapter-trimmed reads were analyzed at the exact-sequence level. Identical sequences were counted independently in each sample, and counts per million reads (CPM) were calculated using the total number of FASTA records in each sample as the library size. A sequence was considered detected in a sample when its raw count was at least 2 and its CPM was at least 1. Within each biological group, sequences meeting both criteria in at least the larger of two samples or one-half of the group size rounded upward were retained. The final sequence set was defined as the union of sequences passing this criterion in at least one group.

The wording should be adjusted to match the preprocessing and sample grouping used in each study.

---

## 25. Citation

When using this script, cite:

1. the associated research article describing the sequence-based analysis
2. the archived software release DOI
3. any reference databases or downstream tools used after filtering

Release-specific citation information should be provided in the repository-level `CITATION.cff` file.

---

## 26. License

This script is distributed under the license provided in the root directory of the repository.

---

## 27. File location in this repository

Recommended repository structure:

```text
v2_pipeline\
└─ filtering\
   ├─ filter_sequences.py
   └─ README.md
```

The groups CSV and input FASTA data are user-supplied and are not hard-coded into the script.
