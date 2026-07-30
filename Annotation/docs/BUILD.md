# Building `annotation_align_core`

## Tested environment

The research implementation and the supplied prebuilt binary were tested with:

- Windows 64-bit
- CPython 3.13
- Cython 3.2.4
- setuptools 82.0.1
- wheel 0.46.3

The only prebuilt binary supplied by this project is:

```text
CPython_3.13_Windows_64-bit/annotation_align_core.cp313-win_amd64.pyd
```

This binary is intended only for 64-bit CPython 3.13 on Windows. It must not
be presented as compatible with CPython 3.11, CPython 3.12, 32-bit Python,
Linux, or macOS.

Source builds may work with other CPython versions that satisfy the package
metadata, but those versions have not been validated unless they are
explicitly listed as tested in the release documentation.

## Required source filenames

Use the following public filenames:

```text
annotation_align_core.pyx
annotation_align_core.c
setup.py
pyproject.toml
```

Both build routes create an importable extension named:

```python
annotation_align_core
```

The annotation script therefore imports:

```python
from annotation_align_core import find_best_hits_multi_payloads_cy
```

## Important step when renaming the original sources

Rename:

```text
annotation_align_core_core_precomputed.pyx
```

to:

```text
annotation_align_core.pyx
```

Then regenerate `annotation_align_core.c` with Cython 3.2.4. Regeneration is
preferred to simply renaming the old generated C file because it ensures that
the generated module initialization symbol and source metadata correspond to
the public module name `annotation_align_core`.

## Route A: standard build from `.pyx`

Install or activate CPython and Microsoft C++ Build Tools, then run:

```bat
python -m pip install --upgrade pip
python -m pip install -r requirements-build.txt
python -m pip install .
```

The `[build-system]` table in `pyproject.toml` also declares the isolated
build requirements used by `pip`.

Import test:

```bat
python -c "from annotation_align_core import find_best_hits_multi_payloads_cy; print('annotation_align_core: OK')"
```

## Regenerate the distributed C source

With Cython 3.2.4 installed:

```bat
python -m cython annotation_align_core.pyx -3 -o annotation_align_core.c
```

Confirm the version used:

```bat
python -c "import Cython; print(Cython.__version__)"
```

The expected result for this release is:

```text
3.2.4
```

Commit both `annotation_align_core.pyx` and the regenerated
`annotation_align_core.c`.

## Route B: build from generated C without Cython

This route still requires:

- a compatible CPython installation
- Microsoft C++ Build Tools on Windows
- setuptools
- wheel

It does not require Cython.

### Windows Command Prompt

```bat
python -m pip install -r requirements-build-c.txt
set ANNOTATION_ALIGN_BUILD_FROM_C=1
python -m pip install --no-build-isolation .
```

### PowerShell

```powershell
python -m pip install -r requirements-build-c.txt
$env:ANNOTATION_ALIGN_BUILD_FROM_C = "1"
python -m pip install --no-build-isolation .
```

Import test:

```bat
python -c "from annotation_align_core import find_best_hits_multi_payloads_cy; print('annotation_align_core: OK')"
```

`--no-build-isolation` is required for this route because the normal
`pyproject.toml` build environment intentionally includes Cython 3.2.4 for
the standard `.pyx` build.

## Route C: use the prebuilt CPython 3.13 Windows binary

For 64-bit CPython 3.13 on Windows only, copy:

```text
annotation_align_core.cp313-win_amd64.pyd
```

into the same directory as `annotate_sequences.py`, or install/copy it into
a directory on the active Python environment's import path.

Test:

```bat
python -c "import sys; print(sys.version); print(sys.maxsize > 2**32)"
python -c "from annotation_align_core import find_best_hits_multi_payloads_cy; print('annotation_align_core: OK')"
```

The first command should report Python 3.13 and `True` for 64-bit Python.

## Run the public annotation script

Prepare a reference-definition CSV with:

```csv
reference_name,file_path
hsa_mature_miR,references/hsa-mature_miR.fa
cfa_mature_miR,references/cfa-mature_miR.fa
```

Then run:

```bat
python annotate_sequences.py ^
  --query-fasta "D:\project\filtered_sequences.fa" ^
  --references-file "D:\project\references.csv" ^
  --output-dir "D:\project\annotation_results" ^
  --core-name "project_sequences" ^
  --workers 16
```

Relative reference paths are resolved relative to the directory containing
the reference-definition CSV.
