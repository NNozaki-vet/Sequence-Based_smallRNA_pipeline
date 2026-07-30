# Annotation End-to-End Regression Test

This directory contains a minimal end-to-end regression test for `annotate_sequences.py`. The `test_annotation.py` script runs the annotation workflow using the included query and reference FASTA files, then compares the generated CSV files with fixed, manually reviewed expected outputs.

## Directory structure

```text
tests/
├─ README.md
├─ requirements-test.txt
├─ test_annotation.py
├─ data/
│  ├─ queries.fa
│  ├─ references.csv
│  └─ references/
│     ├─ reference_a.fa
│     ├─ reference_b.fa
│     └─ reference_c.fa
└─ expected/
   ├─ test_sequences_reference_a_minimumMiss.csv
   ├─ test_sequences_reference_b_minimumMiss.csv
   └─ test_sequences_reference_c_minimumMiss.csv
```

## What the test verifies

The test verifies that:

- the compiled `annotation_align_core` extension can be imported and used
- `annotate_sequences.py` can be executed through its command-line interface
- the query FASTA, reference-definition CSV, and reference FASTA files are loaded correctly
- one output CSV file is created for each reference definition
- output column names match the expected files exactly
- all output rows match the fixed expected results
- duplicate-row counts are preserved
- differences in output-row order alone do not cause failure

CSV rows are compared as multisets. Missing rows, unexpected rows, or incorrect duplicate counts cause the test to fail.

## Scope and limitations

This test checks software behavior and output reproducibility for the included test data.

It does not independently verify:

- biological validity of the sequences or annotations
- performance or memory usage
- compatibility with untested operating systems or Python versions
- equivalence of multiprocessing behavior across different worker counts
- successful compilation from `annotation_align_core.pyx` or `annotation_align_core.c`

Build procedures and platform compatibility must be validated separately.

## Requirements

The test uses only the Python standard library. No third-party testing framework is required.

The active Python interpreter must be able to import:

```python
from annotation_align_core import find_best_hits_multi_payloads_cy
```

For the provided prebuilt binary, the supported environment is:

```text
64-bit Microsoft Windows
64-bit CPython 3.13
annotation_align_core.cp313-win_amd64.pyd
```

Place the compatible `.pyd` file beside `annotate_sequences.py`, or install the extension into the active Python environment.

Confirm the extension path before running the test:

```bat
python -c "import annotation_align_core; print(annotation_align_core.__file__)"
```

The displayed path should identify the intended released or locally built extension, not an older development build.

## Test data

The test uses:

- `data/queries.fa`
- `data/references.csv`
- `data/references/reference_a.fa`
- `data/references/reference_b.fa`
- `data/references/reference_c.fa`

The reference-definition CSV must contain:

```csv
reference_name,file_path
reference_a,references/reference_a.fa
reference_b,references/reference_b.fa
reference_c,references/reference_c.fa
```

Relative reference paths are resolved from the directory containing `references.csv`.

The included sequences are test fixtures for validating the annotation workflow. They are not a biological reference database.

## Expected outputs

The `expected` directory contains the fixed reference outputs:

```text
test_sequences_reference_a_minimumMiss.csv
test_sequences_reference_b_minimumMiss.csv
test_sequences_reference_c_minimumMiss.csv
```

These filenames correspond to:

```text
core_name = test_sequences
reference_name = reference_a, reference_b, or reference_c
workers = 1
```

The expected files must be manually reviewed before inclusion. They must not be regenerated automatically during routine testing, because doing so could incorrectly redefine a software defect as expected behavior.

## Run

From the `Annotation` directory:

```bat
python tests\test_annotation.py
```

The script can also be executed by absolute path from another working directory.

During execution, `test_annotation.py`:

1. validates the required test files
2. creates a temporary output directory
3. runs `annotate_sequences.py` as a subprocess
4. confirms that all expected output files are created
5. compares generated CSV headers and rows with the fixed expected files
6. removes the temporary output directory automatically

A successful run ends with:

```text
Annotation integration test: PASSED
```

Any missing file, import failure, command failure, header difference, missing row, unexpected row, or duplicate-count difference causes a nonzero exit status.

## Updating expected outputs

Update the expected CSV files only after an intentional change to the annotation method or output format.

Before replacing them:

1. inspect every changed row
2. manually verify mismatch, extension, unmatched-reference, tied-hit, and unannotated cases
3. document the reason for the change
4. update the relevant method and output documentation
5. rerun the test with a clean temporary output directory
6. commit the source and expected-output changes together

A change caused only by a different row order does not require replacing the expected files.

## Release check

Before publishing a release:

1. clear any temporary `PYTHONPATH` setting that points to a development build
2. confirm the imported `annotation_align_core` path
3. run this test with the released prebuilt binary
4. run it again with a clean source-built extension when source-build support is claimed
5. record the Python version, platform, extension path, release tag, and commit hash

This test should pass before the release is archived or cited.
