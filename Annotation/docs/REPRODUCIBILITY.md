# Reproducibility

## Tested scope

The current release was tested only with:

```text
64-bit Windows
CPython 3.13
Cython 3.2.4
setuptools 82.0.1
wheel 0.46.3
```

Do not describe Python 3.11 or 3.12 as tested until the complete build, import test, toy test, and output comparison have been performed.

## Files that must be archived

Archive the exact versions of:

```text
annotate_sequences.py
annotation_align_core.pyx
annotation_align_core.c
setup.py
pyproject.toml
requirements-build.txt
requirements-build-c.txt
references.csv
all query FASTA files
all reference FASTA files
final annotation CSV files
run log
```

For the binary route, also archive:

```text
annotation_align_core.cp313-win_amd64.pyd
```

Do not treat `reference_cache` or temporary worker CSV files as primary reproducibility materials.

## Required run record

Record:

- software release tag and commit hash
- operating system edition and architecture
- Python implementation and full version
- whether the extension came from `.pyx`, `.c`, or prebuilt `.pyd`
- compiler and build tools when compiled
- Cython, setuptools, and wheel versions
- full command line
- worker count
- query and reference database provenance
- reference database release/version and retrieval date
- any transformation used to create sense, antisense, precursor, loop, or other reference FASTA files
- start and completion dates
- all warnings and errors

Capture the run log:

```bat
python annotate_sequences.py ^
  --query-fasta "D:\project\queries.fa" ^
  --references-file "D:\project\references.csv" ^
  --output-dir "D:\project\annotation_output" ^
  --core-name "project_sequences" ^
  --workers 16 > "D:\project\annotation_output\run.log" 2>&1
```

## Environment record

```bat
python --version
python -c "import sys, platform; print(sys.executable); print(platform.platform()); print(platform.architecture())"
python -c "import annotation_align_core; print(annotation_align_core.__file__)"
python -m pip freeze > environment_pip_freeze.txt
```

For a conda environment:

```bat
conda env export --from-history > environment_from_history.yml
conda list --explicit > environment_explicit.txt
```

## SHA-256 checksums

Create checksums before analysis and after collecting outputs:

```bat
certutil -hashfile annotate_sequences.py SHA256
certutil -hashfile annotation_align_core.pyx SHA256
certutil -hashfile annotation_align_core.c SHA256
certutil -hashfile annotation_align_core.cp313-win_amd64.pyd SHA256
certutil -hashfile queries.fa SHA256
certutil -hashfile reference.fa SHA256
certutil -hashfile output.csv SHA256
```

Store the values in a plain-text manifest. A filename alone is not evidence of identical content.

## Formal reproduction procedure

1. Create a clean environment.
2. Check out the exact release or commit.
3. verify all input and source checksums.
4. Build or select the documented compiled core.
5. confirm the extension import path.
6. use a new output directory.
7. execute the archived command.
8. preserve the complete log.
9. compare output row content and SHA-256 checksums.

If byte-identical output is not obtained, compare parsed CSV rows after preserving all tied rows. Investigate environment, input, binary, and cache differences before concluding that the algorithm changed.

## Validation before release

A release should not claim broader compatibility until it passes:

- clean installation
- extension import test
- toy query/reference test
- expected-output comparison
- full analysis smoke test
- source-build and prebuilt-binary comparison
- repeated-run comparison with a fresh cache

## Citation and versioning

Cite the archived release that contains the exact code used. Do not cite an older release DOI as though it archives a later development version. Record the release tag, commit hash, and archive DOI together.
