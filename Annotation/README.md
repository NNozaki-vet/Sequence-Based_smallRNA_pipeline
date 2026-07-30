# Sequence-Based Annotation

A Cython-accelerated workflow for annotating exact small RNA sequences against user-defined reference FASTA files.

This directory contains the version 2 annotation implementation of the [Sequence-Based Small RNA Pipeline](../README.md). It is the recommended annotation workflow for new analyses. The historical version 1 workflow, which aggregates annotation tables exported from CLC Genomics Workbench, is retained separately in [`../Original_code/`](../Original_code/) for manuscript reproduction.

## Overview

`annotate_sequences.py` compares each distinct query sequence with one or more reference FASTA files while retaining sequence-level information.

The workflow:

* accepts an exact-sequence query FASTA file
* accepts any number of user-defined reference FASTA files
* uses FASTA headers as reference entry names
* supports multiple names associated with an identical reference sequence
* collapses duplicate query sequences while preserving first-occurrence order
* evaluates mismatches and terminal query extensions
* retains tied best hits under the implemented ranking rules
* writes one minimum-score CSV file per reference
* caches preprocessed reference data
* uses multiprocessing
* uses the compiled `annotation_align_core` extension for the computational search

Reference databases are not included. Users must obtain, document, and cite the reference sequences used in their analyses.

## Recommended use

Use this version 2 module when:

* starting a new sequence-based small RNA analysis
* annotating against custom or updated reference FASTA files
* comparing the same query sequences against multiple reference databases
* avoiding dependence on CLC Genomics Workbench
* requiring documented command-line inputs and outputs
* requiring a fixed end-to-end regression test

Use the version 1 annotation utilities in [`../Original_code/`](../Original_code/) only when reproducing the original manuscript workflow based on CLC Genomics Workbench exports.

## Directory structure

```text
Annotation/
├─ README.md
├─ annotate_sequences.py
├─ annotation_align_core.pyx
├─ annotation_align_core.c
├─ setup.py
├─ pyproject.toml
├─ requirements-build.txt
├─ requirements-build-c.txt
├─ Cpython_3.13_Windows_64-bit/
│  ├─ README.md
│  └─ annotation_align_core.cp313-win_amd64.pyd
├─ docs/
│  ├─ ANNOTATION_METHOD.md
│  ├─ CONFIGURATION.md
│  ├─ INPUT_FORMATS.md
│  ├─ OUTPUT_FORMATS.md
│  ├─ REPRODUCIBILITY.md
│  └─ TROUBLESHOOTING.md
└─ tests/
   ├─ README.md
   ├─ requirements-test.txt
   ├─ test_annotation.py
   ├─ data/
   └─ expected/
```

## Requirements

### Prebuilt Windows extension

The included prebuilt extension is intended only for:

* 64-bit Microsoft Windows
* 64-bit CPython 3.13
* the matching file:
  `annotation_align_core.cp313-win_amd64.pyd`

No compiler is required when the compatible prebuilt extension is used.

### Source build

Building from source requires:

* a compatible Python implementation
* a C compiler compatible with that Python installation
* `setuptools`
* `wheel`
* Cython when building from `annotation_align_core.pyx`

The generated `annotation_align_core.c` file can be used to build without Cython, but a compatible C compiler is still required.

The released Windows binary was tested with:

```text
64-bit Microsoft Windows
64-bit CPython 3.13
Cython 3.2.4
setuptools 82.0.1
wheel 0.46.3
```

## Quick start with the prebuilt Windows extension

The compiled extension must be importable as:

```python
from annotation_align_core import find_best_hits_multi_payloads_cy
```

The simplest method is to copy the compatible `.pyd` file beside `annotate_sequences.py`.

From the repository root:

```bat
copy Annotation\Cpython_3.13_Windows_64-bit\annotation_align_core.cp313-win_amd64.pyd Annotation\
```

Confirm the imported file:

```bat
python -c "import annotation_align_core; print(annotation_align_core.__file__)"
```

The displayed path should point to the intended released `.pyd`, not an older development build.

No hard-coded `sys.path` modification is required.

## Input files

The workflow requires:

1. one query FASTA file
2. one reference-definition CSV file
3. one or more reference FASTA files
4. one output directory

### Query FASTA

The query FASTA contains the sequences to annotate.

Example:

```fasta
>query_1
TGAGGTAGTAGGTTGTATAGTT
>query_2
TAGCTTATCAGACTGATGTTGA
```

Important behavior:

* FASTA headers are not used as query identifiers
* query sequences are converted to uppercase
* `U` is converted to `T`
* spaces and line breaks within sequences are removed
* duplicate query sequences are collapsed
* first-occurrence order is retained
* line-wrapped sequences are supported
* plain-text FASTA is required
* compressed `.gz` input is not supported

The query FASTA may be produced by the version 2 [`Filtering`](../Filtering/) module or by another validated preprocessing workflow.

### Reference-definition CSV

The CSV must contain these exact column names:

```csv
reference_name,file_path
```

Example:

```csv
reference_name,file_path
mirna,references\mirna_reference.fa
other_ncrna,references\other_ncrna_reference.fa
```

Rules:

* `reference_name` must be unique
* `file_path` may be absolute or relative
* relative paths are resolved from the directory containing the CSV file
* blank rows are ignored
* UTF-8 and UTF-8 with BOM are supported

Each `reference_name` becomes part of the corresponding output filename.

### Reference FASTA

Example:

```fasta
>cfa-miR-example-1
TGAGGTAGTAGGTTGTATAGTT
>cfa-miR-example-2
TAGCTTATCAGACTGATGTTGA
```

Important behavior:

* each header becomes a `reference_entry_name`
* empty headers are rejected
* empty sequences are rejected
* sequence data before the first header is rejected
* line-wrapped sequences are supported
* identical reference sequences may have multiple entry names
* sequences are converted to uppercase
* `U` is converted to `T`
* plain-text FASTA is required
* compressed `.gz` input is not supported

For complete specifications, see [`docs/INPUT_FORMATS.md`](docs/INPUT_FORMATS.md).

## Command-line usage

```text
python annotate_sequences.py
    --query-fasta QUERY_FASTA
    --references-file REFERENCES_CSV
    --output-dir OUTPUT_DIRECTORY
    [--core-name OUTPUT_PREFIX]
    [--workers NUMBER_OF_WORKERS]
```

### Required arguments

| Argument            | Description                                           |
| ------------------- | ----------------------------------------------------- |
| `--query-fasta`     | Query FASTA containing exact sequences to annotate    |
| `--references-file` | CSV defining reference names and FASTA paths          |
| `--output-dir`      | Directory for annotation outputs and reference caches |

### Optional arguments

| Argument      |                                       Default | Description                       |
| ------------- | --------------------------------------------: | --------------------------------- |
| `--core-name` | Query FASTA filename without its FASTA suffix | Prefix used for output filenames  |
| `--workers`   |                                          `16` | Number of multiprocessing workers |

`--workers` must be at least 1.

## Minimal example

From the repository root:

```bat
python Annotation\annotate_sequences.py ^
  --query-fasta "filtering_results\project_filtered_sequences.fa" ^
  --references-file "references.csv" ^
  --output-dir "annotation_results" ^
  --core-name "project_sequences" ^
  --workers 16
```

PowerShell example:

```powershell
python .\Annotation\annotate_sequences.py `
  --query-fasta ".\filtering_results\project_filtered_sequences.fa" `
  --references-file ".\references.csv" `
  --output-dir ".\annotation_results" `
  --core-name "project_sequences" `
  --workers 16
```

## Output files

One CSV file is written for each reference definition:

```text
<core_name>_<reference_name>_minimumMiss.csv
```

Example:

```text
project_sequences_mirna_minimumMiss.csv
project_sequences_other_ncrna_minimumMiss.csv
```

The output records:

* query sequence
* reference entry name
* total mismatch and extension count
* mismatch count
* 5-prime query extension length and sequence
* unmatched 5-prime reference length
* up to three mismatch positions and substitutions
* 3-prime query extension length and sequence
* unmatched 3-prime reference length

Queries without a retained hit are written as:

```text
reference_entry_name = Unannotated
total_mismatch_and_extension_count = 10
```

For exact column definitions, see [`docs/OUTPUT_FORMATS.md`](docs/OUTPUT_FORMATS.md).

## Matching and ranking

The default search limits are:

```text
maximum mismatch count: 3
maximum total query extension: 3
```

Candidate relationships are ranked by:

1. mismatch count
2. total extension length (`5-prime extension + 3-prime extension`)
3. 5-prime extension length

Lower values are preferred.

All tied best hits retained under these rules may be reported. An annotation result therefore does not necessarily represent a unique reference entry.

The score reflects sequence similarity under the implemented rules. It is not, by itself, evidence of RNA biogenesis, biological identity, or function.

For the complete method, see [`docs/ANNOTATION_METHOD.md`](docs/ANNOTATION_METHOD.md).

## Reference cache

Preprocessed references are stored under:

```text
<output-dir>\reference_cache\
```

Cache validity is checked using:

* cache format version
* absolute source path
* source modification time
* source file size

A cache is rebuilt automatically when these values no longer match.

Delete the corresponding cache file when a forced rebuild is required.

## Building the extension

### Build from the Cython source

From the `Annotation` directory:

```bat
python -m pip install -r requirements-build.txt
python setup.py build_ext --inplace
```

This route compiles:

```text
annotation_align_core.pyx
```

### Build from the generated C source

From the `Annotation` directory in Windows Command Prompt:

```bat
python -m pip install -r requirements-build-c.txt
set ANNOTATION_ALIGN_BUILD_FROM_C=1
python setup.py build_ext --inplace
set ANNOTATION_ALIGN_BUILD_FROM_C=
```

PowerShell:

```powershell
python -m pip install -r requirements-build-c.txt
$env:ANNOTATION_ALIGN_BUILD_FROM_C = "1"
python .\setup.py build_ext --inplace
Remove-Item Env:ANNOTATION_ALIGN_BUILD_FROM_C
```

This route compiles:

```text
annotation_align_core.c
```

Cython is not required for the generated-C route, but a compatible C compiler is required.

After either build, confirm the imported extension:

```bat
python -c "import annotation_align_core; print(annotation_align_core.__file__)"
```

See also [`Cpython_3.13_Windows_64-bit/README.md`](Cpython_3.13_Windows_64-bit/README.md).

## End-to-end regression test

The included regression test runs the annotation workflow with fixed query and reference FASTA files and compares the generated CSV files with manually reviewed expected outputs.

From the `Annotation` directory, make the included prebuilt extension importable:

### Windows Command Prompt

```bat
set PYTHONPATH=%CD%\Cpython_3.13_Windows_64-bit
python tests\test_annotation.py
set PYTHONPATH=
```

### PowerShell

```powershell
$env:PYTHONPATH = "$PWD\Cpython_3.13_Windows_64-bit"
python .\tests\test_annotation.py
Remove-Item Env:PYTHONPATH
```

A successful run ends with:

```text
Annotation integration test: PASSED
```

The test verifies:

* import and use of the compiled extension
* execution of `annotate_sequences.py`
* loading of the included input files
* creation of one output CSV per reference
* exact output-header agreement
* output-row agreement while preserving duplicate counts
* immunity to output-row order alone

The test does not independently establish biological validity, performance at production scale, or compatibility with untested platforms.

See [`tests/README.md`](tests/README.md).

## Documentation

| Document                                                                         | Purpose                                                      |
| -------------------------------------------------------------------------------- | ------------------------------------------------------------ |
| [`docs/ANNOTATION_METHOD.md`](docs/ANNOTATION_METHOD.md)                         | Matching logic, ranking, tied hits, and method details       |
| [`docs/INPUT_FORMATS.md`](docs/INPUT_FORMATS.md)                                 | Query FASTA, reference CSV, and reference FASTA requirements |
| [`docs/OUTPUT_FORMATS.md`](docs/OUTPUT_FORMATS.md)                               | Output filenames, columns, and interpretation                |
| [`docs/CONFIGURATION.md`](docs/CONFIGURATION.md)                                 | Runtime settings and configuration                           |
| [`docs/TROUBLESHOOTING.md`](docs/TROUBLESHOOTING.md)                             | Common import, build, path, and runtime problems             |
| [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md)                             | Information to record for reproducible analyses              |
| [`tests/README.md`](tests/README.md)                                             | Regression-test design and execution                         |
| [`Cpython_3.13_Windows_64-bit/README.md`](Cpython_3.13_Windows_64-bit/README.md) | Prebuilt Windows binary compatibility and use                |

## Scope and limitations

This module does not perform:

* adapter trimming
* FASTQ quality filtering
* FASTQ-to-FASTA conversion
* abundance filtering
* reference database download
* automatic biological classification
* RNA secondary-structure prediction
* differential-expression testing
* multiple-testing correction
* pathway analysis
* target prediction
* automatic integration of outputs from multiple references

Annotation depends on the supplied reference FASTA files. A sequence absent from a reference cannot receive an annotation from that reference.

Runtime and memory use depend on:

* number of query sequences
* query lengths
* reference size
* reference sequence lengths
* number of references
* number of workers
* frequency of near matches and tied hits

The included prebuilt `.pyd` is not portable across unsupported Python versions, operating systems, or architectures.

## Reproducibility

For each analysis, record:

* repository release, branch, and commit hash
* operating system
* Python implementation and version
* compiled extension filename and build route
* worker count
* query FASTA-generation procedure
* query FASTA checksum
* reference database name and version
* reference download date
* reference FASTA checksum
* reference-definition CSV
* command-line arguments
* output directory and core name
* downstream integration and statistical criteria

For publication, use a tagged release rather than an unpinned development branch whenever a version 2 release is available.

## License

This software is distributed under the repository-level [Apache License 2.0](../LICENSE).

Reference databases and third-party software remain subject to their own licenses and terms of use.

## Citation

A release-specific citation and archived DOI will be provided at the repository level.

Until then, report the exact repository branch and commit hash used. Also cite:

1. the associated research article
2. the software release or exact commit
3. every reference database used
4. relevant preprocessing and downstream statistical software

See the repository-level [`README.md`](../README.md) for current citation information.
