# Sequence-Based Small RNA Pipeline

A sequence-centered workflow for filtering and annotating small RNA-seq reads while preserving exact nucleotide sequences.

> **Recommended implementation: version 2**
>
> New analyses should use the generalized workflows in [`Filtering/`](Filtering/) and [`Annotation/`](Annotation/) from the [`v2-development`](https://github.com/NNozaki-vet/Sequence-Based_smallRNA_pipeline/tree/v2-development) branch.
>
> The [`Original_code/`](Original_code/) directory is retained as the version 1 research implementation for manuscript reproduction and historical reference. It is not the recommended starting point for new projects.

## Overview

Conventional small RNA-seq workflows often summarize reads under predefined RNA or miRNA names early in the analysis. This repository instead preserves each distinct read sequence as an analytical unit.

The approach is designed to:

* retain exact nucleotide sequences
* distinguish sequence variants that may otherwise be collapsed under the same annotation
* reduce dependence on a single predefined annotation set
* support user-selected reference FASTA files
* retain canonical miRNAs and noncanonical or precursor-derived sequences
* capture fragments derived from other small or noncoding RNAs, including tRNA, rRNA, Y RNA, snoRNA, snRNA, and lncRNA
* generate sequence-level outputs suitable for downstream statistical analysis

Potential applications include the analysis of:

* mature miRNAs
* isomiRs
* seed-shifted miRNA sequences
* miRNA-offset RNAs (moRNAs)
* pre-miRNA-derived fragments
* antisense or mirror-derived sequences
* loop-derived fragments
* tRNA-, rRNA-, Y RNA-, snoRNA-, snRNA-, and lncRNA-derived small RNAs
* unannotated small RNA sequences

The repository is associated with the manuscript:

> “Exact-sequence profiling reveals conserved small RNA fragment remodeling in hepatocellular carcinoma across dogs, humans, and mice”
> Manuscript link forthcoming.

## Recommended workflow

```text
Adapter-trimmed, quality-controlled read-level FASTA files
                         |
                         v
                 Filtering/filter_sequences.py
                         |
                         +--> exact-sequence raw-count matrix
                         +--> filtered exact-sequence FASTA
                         +--> group-pass and sample metadata files
                         |
                         v
                 Annotation/annotate_sequences.py
                         |
                         +--> one minimum-score annotation CSV per reference
                         |
                         v
        Downstream annotation integration and statistical analysis
```

Adapter trimming, read-quality control, FASTQ-to-FASTA conversion, differential-expression testing, and biological interpretation are outside the scope of the version 2 scripts and must be performed separately.

## Repository structure

```text
Sequence-Based_smallRNA_pipeline/
├─ README.md
├─ LICENSE
├─ Filtering/
│  ├─ README.md
│  └─ filter_sequences.py
├─ Annotation/
│  ├─ annotate_sequences.py
│  ├─ annotation_align_core.c
│  ├─ setup.py
│  ├─ pyproject.toml
│  ├─ requirements-build.txt
│  ├─ requirements-build-c.txt
│  ├─ Cpython_3.13_Windows_64-bit/
│  │  ├─ README.md
│  │  └─ annotation_align_core.cp313-win_amd64.pyd
│  ├─ docs/
│  └─ tests/
└─ Original_code/
```

## Version 1 and version 2

| Feature              | Version 1: `Original_code/`                                              | Version 2: `Filtering/` and `Annotation/`                                                       |
| -------------------- | ------------------------------------------------------------------------ | ----------------------------------------------------------------------------------------------- |
| Primary purpose      | Reproduce the original manuscript analysis                               | Provide a generalized workflow for new datasets                                                 |
| Status               | Archived research implementation                                         | Recommended implementation                                                                      |
| Input paths          | Manuscript- and environment-specific assumptions may remain              | Supplied through command-line arguments and CSV configuration files                             |
| Filtering            | Original study-specific scripts and utilities                            | Stand-alone exact-sequence filtering with user-defined biological groups                        |
| Annotation           | Aggregation of CLC Genomics Workbench exports                            | Direct annotation against user-defined reference FASTA files                                    |
| Proprietary software | CLC Genomics Workbench v22 is required for the original annotation route | CLC Genomics Workbench is not required                                                          |
| Reference selection  | Tied to the original analytical workflow                                 | User-defined reference names, FASTA files, and processing order                                 |
| Computational core   | Original implementation                                                  | Cython-accelerated sequence comparison                                                          |
| Documentation        | Primarily manuscript-reproduction oriented                               | Input, output, method, configuration, build, troubleshooting, and reproducibility documentation |
| Tests                | No fixed end-to-end regression suite                                     | Fixed, manually reviewed end-to-end regression test                                             |
| Portability          | Limited                                                                  | Generalized command-line interface; platform limitations are documented                         |

Version 1 remains available because it records the implementation used during the original research. Version 2 was reorganized and generalized to improve reuse, transparency, portability, and reproducibility.

## Version 2 module 1: exact-sequence filtering

The [`Filtering/`](Filtering/) module reads adapter-trimmed FASTA files in which:

```text
one FASTA record = one sequencing read
one FASTA file   = one biological sample
one folder       = one biological group
one CSV row      = one biological group
```

For each sample, `filter_sequences.py`:

1. counts identical nucleotide sequences
2. calculates counts per million reads (CPM)
3. applies the sample-level criteria:

   * raw count >= 2
   * CPM >= 1
4. requires these criteria to be met in at least:

```text
max(2, ceiling(number of samples in the group / 2))
```

5. retains the union of sequences that pass in at least one biological group
6. writes a filtered FASTA file, a raw-count matrix, and supporting metadata

The group-level criterion is calculated separately for each group. Groups may contain different numbers of samples.

### Filtering requirements

* Python 3.9 or later is recommended
* only the Python standard library is required
* Windows instructions are provided
* the script may also run on Linux or macOS
* input must be read-level FASTA, not abundance-collapsed FASTA

### Minimal filtering example

Create a CSV file such as:

```csv
group,folder
Control,input_fasta\Control
Disease,input_fasta\Disease
```

Then run:

```bat
python Filtering\filter_sequences.py ^
  --groups-file "group_directories.csv" ^
  --output-dir "filtering_results" ^
  --output-prefix "project_filtered_sequences" ^
  --threads 8
```

For complete input requirements, command-line options, filtering equations, output definitions, and troubleshooting, see:

* [`Filtering/README.md`](Filtering/README.md)

## Version 2 module 2: sequence-based annotation

The [`Annotation/`](Annotation/) module compares each filtered query sequence with one or more user-defined reference FASTA files.

The annotation workflow:

* preserves query sequences as distinct analytical units
* accepts any number of named reference FASTA files
* resolves relative reference paths from a configuration CSV
* supports line-wrapped reference FASTA records
* caches preprocessed reference data
* uses multiprocessing
* uses a compiled Cython extension for the computational search
* reports all tied best hits retained under the ranking rules
* writes one minimum-score output CSV for each reference definition

The matching procedure evaluates sequence relationships using:

* mismatch count
* total 5-prime and 3-prime query extension
* 5-prime extension as the final ranking component

The default search limits are:

```text
maximum mismatch count: 3
maximum total terminal extension: 3
```

Detailed definitions are provided in:

* [`Annotation/docs/ANNOTATION_METHOD.md`](Annotation/docs/ANNOTATION_METHOD.md)
* [`Annotation/docs/INPUT_FORMATS.md`](Annotation/docs/INPUT_FORMATS.md)
* [`Annotation/docs/OUTPUT_FORMATS.md`](Annotation/docs/OUTPUT_FORMATS.md)
* [`Annotation/docs/CONFIGURATION.md`](Annotation/docs/CONFIGURATION.md)
* [`Annotation/docs/BUILD.md`](Annotation/docs/BUILD.md)
* [`Annotation/docs/TROUBLESHOOTING.md`](Annotation/docs/TROUBLESHOOTING.md)
* [`Annotation/docs/REPRODUCIBILITY.md`](Annotation/docs/REPRODUCIBILITY.md)

### Reference configuration

Create a CSV file with the required columns:

```csv
reference_name,file_path
mirna,references\mirna_reference.fa
other_ncrna,references\other_ncrna_reference.fa
```

Reference FASTA headers are used as annotation entry names.

Reference databases are not bundled with this repository. Users are responsible for:

* obtaining reference sequences from appropriate sources
* complying with database licenses and terms of use
* recording database names, versions, release dates, and download dates
* citing all reference databases used in the analysis

### Prebuilt Windows extension

A prebuilt extension is supplied at:

```text
Annotation/Cpython_3.13_Windows_64-bit/
annotation_align_core.cp313-win_amd64.pyd
```

It is intended only for:

* Microsoft Windows
* 64-bit CPython 3.13
* 64-bit operating systems

For the simplest setup, copy the `.pyd` file into the directory containing `annotate_sequences.py`.

```bat
copy Annotation\Cpython_3.13_Windows_64-bit\annotation_align_core.cp313-win_amd64.pyd Annotation\
```

Confirm the import:

```bat
python -c "import annotation_align_core; print(annotation_align_core.__file__)"
```

Users of other Python versions or operating systems must build the extension for their own environment. See:

* [`Annotation/docs/BUILD.md`](Annotation/docs/BUILD.md)
* [`Annotation/Cpython_3.13_Windows_64-bit/README.md`](Annotation/Cpython_3.13_Windows_64-bit/README.md)

### Minimal annotation example

```bat
python Annotation\annotate_sequences.py ^
  --query-fasta "filtering_results\project_filtered_sequences.fa" ^
  --references-file "references.csv" ^
  --output-dir "annotation_results" ^
  --core-name "project_sequences" ^
  --workers 16
```

Output files follow this pattern:

```text
<core_name>_<reference_name>_minimumMiss.csv
```

## Regression test

Version 2 includes a fixed end-to-end regression test for the annotation workflow.

The test:

* runs `annotate_sequences.py` using included query and reference FASTA files
* checks creation of all expected output files
* compares output headers exactly
* compares output rows as multisets
* preserves duplicate-row counts
* ignores row-order differences alone
* uses fixed, manually reviewed expected CSV files

On 64-bit Windows with CPython 3.13, make the prebuilt extension importable and run:

```bat
set PYTHONPATH=%CD%\Annotation\Cpython_3.13_Windows_64-bit
python Annotation\tests\test_annotation.py
```

A successful run ends with:

```text
Annotation integration test: PASSED
```

The current `v2-development` branch has also been validated from a fresh Git clone using the included prebuilt extension and test data.

See:

* [`Annotation/tests/README.md`](Annotation/tests/README.md)

## Scope and limitations

Version 2 does not perform:

* sequencing-adapter trimming
* FASTQ quality filtering
* FASTQ-to-FASTA conversion
* read-length filtering unless performed before input
* removal of contaminating RNA classes
* de novo RNA discovery
* RNA secondary-structure prediction
* differential-expression testing
* multiple-testing correction
* pathway analysis
* target prediction
* automatic integration of outputs from multiple references
* biological classification of every annotated sequence

Annotation results depend on the reference FASTA files and database versions supplied by the user. A sequence absent from a reference cannot receive that reference's annotation.

The annotation score describes sequence similarity under the implemented mismatch and terminal-extension rules. It is not, by itself, evidence of RNA biogenesis, function, or biological identity.

## Reproducibility recommendations

For each analysis, record:

* repository version, release tag, branch, and commit hash
* operating system
* Python implementation and version
* number of worker processes
* filtering configuration
* input sample groups
* input FASTA-generation procedure
* reference database names and versions
* exact reference FASTA files
* reference-definition CSV
* output filename prefix or core name
* any downstream filtering or statistical criteria

For formal publications, use a tagged release rather than an unpinned development branch whenever a version 2 release is available.

## Which version should I use?

Use version 2 when:

* starting a new analysis
* working with new sample groups or paths
* avoiding dependence on CLC Genomics Workbench
* annotating against custom or updated reference FASTA files
* requiring documented command-line inputs and outputs
* requiring the included regression test

Use version 1 only when:

* reproducing the original manuscript workflow
* inspecting the historical implementation
* reproducing CLC Genomics Workbench export aggregation used in that workflow

In brief:

```text
New analysis                 -> use Filtering/ and Annotation/ from version 2
Original manuscript replay   -> use Original_code/ from version 1
```

## Development status

Version 2 is currently available on the `v2-development` branch while final release preparation is completed:

* https://github.com/NNozaki-vet/Sequence-Based_smallRNA_pipeline/tree/v2-development

After the version 2 release is published, users should prefer the corresponding tagged release over the development branch.

## License

This repository is distributed under the Apache License 2.0.

See:

* [`LICENSE`](LICENSE)

Reference databases and third-party tools remain subject to their own licenses and terms of use.

## Citation

A release-specific citation and archived DOI will be provided for version 2.

Until then, analyses using the development branch should record and report the exact commit hash. The associated manuscript link will be added when available.

When using this workflow, cite:

1. the associated research article
2. the archived software release or exact repository commit
3. every reference database used for annotation
4. relevant third-party preprocessing or statistical software

## Issues and contributions

Questions, bug reports, and reproducible feature requests may be submitted through the repository's GitHub Issues page.

When reporting a problem, include:

* operating system
* Python version
* exact command
* complete error message
* repository branch, tag, or commit hash
* minimal example input when it can be shared
