# Prebuilt Cython Extension for CPython 3.13 on Windows 64-bit

This directory provides a prebuilt binary extension for the sequence-based small RNA annotation workflow.

Included file:

```text
annotation_align_core.cp313-win_amd64.pyd
```

The extension provides the compiled function:

```python
from annotation_align_core import find_best_hits_multi_payloads_cy
```

It is used by `annotate_sequences.py` to accelerate exact-sequence comparison between query small RNA sequences and reference sequences.

---

## Supported environment

This prebuilt binary is intended only for:

- Microsoft Windows
- 64-bit operating systems
- CPython 3.13
- 64-bit CPython installation

The filename tag indicates the supported environment:

```text
cp313       = CPython 3.13
win_amd64   = 64-bit Windows
```

This binary is not expected to work with:

- Python 3.11
- Python 3.12
- Python 3.14 or later
- 32-bit Python
- Linux
- macOS
- PyPy
- other Python implementations

The current release was tested with CPython 3.13 on 64-bit Windows.

---

## No compiler or Cython installation is required

When this prebuilt `.pyd` file is used, users do not need to install:

- Cython
- Microsoft C++ Build Tools
- setuptools for compilation
- wheel for compilation

The `.pyd` file is already compiled.

Python 3.13 for 64-bit Windows is still required.

---

## Recommended directory structure

Place the `.pyd` file in the same directory as `annotate_sequences.py`.

```text
Annotation/
├─ annotate_sequences.py
├─ annotation_align_core.cp313-win_amd64.pyd
├─ references.csv
└─ CPython_3.13_Windows_64-bit/
   └─ README.md
```

Alternatively, the distributed directory may initially contain:

```text
CPython_3.13_Windows_64-bit/
├─ annotation_align_core.cp313-win_amd64.pyd
└─ README.md
```

In that case, copy:

```text
annotation_align_core.cp313-win_amd64.pyd
```

into the directory containing:

```text
annotate_sequences.py
```

The final working directory should contain both files:

```text
annotate_sequences.py
annotation_align_core.cp313-win_amd64.pyd
```

---

## Confirm the Python version

Open Command Prompt, Anaconda Prompt, or Miniforge Prompt and run:

```bat
python --version
```

The output must indicate Python 3.13, for example:

```text
Python 3.13.x
```

Confirm that the Python installation is 64-bit:

```bat
python -c "import sys; print(sys.maxsize > 2**32)"
```

The expected result is:

```text
True
```

You can also display the full platform information:

```bat
python -c "import platform; print(platform.python_implementation()); print(platform.python_version()); print(platform.architecture())"
```

The expected environment is:

```text
CPython
3.13.x
('64bit', 'WindowsPE')
```

---

## Confirm that the extension can be imported

Change to the directory containing both `annotate_sequences.py` and the `.pyd` file.

Example:

```bat
cd /d "D:\smallRNA\cfa_miRNA_moRNA_mirror\v2_pipeline_For_Github\Annotation"
```

Then run:

```bat
python -c "from annotation_align_core import find_best_hits_multi_payloads_cy; print('annotation_align_core import: OK')"
```

Expected result:

```text
annotation_align_core import: OK
```

A more detailed check is:

```bat
python -c "import annotation_align_core; print(annotation_align_core.__file__); print(hasattr(annotation_align_core, 'find_best_hits_multi_payloads_cy'))"
```

Expected output includes the path to:

```text
annotation_align_core.cp313-win_amd64.pyd
```

and:

```text
True
```

---

## Run the annotation workflow

After the import test succeeds, run the annotation script.

Example:

```bat
python annotate_sequences.py ^
  --query-fasta "D:\project\filtered_sequences.fa" ^
  --references-file "D:\project\references.csv" ^
  --output-dir "D:\project\annotation_results" ^
  --core-name "project_sequences" ^
  --workers 16
```

Use:

```bat
python annotate_sequences.py --help
```

to display all available command-line options.

---

## Import behavior

The annotation script imports the compiled extension directly:

```python
from annotation_align_core import find_best_hits_multi_payloads_cy
```

No hard-coded path such as the following is required:

```python
sys.path.append(r"D:\cython_align")
```

Python searches the current script directory and installed module locations automatically.

Keeping the `.pyd` file in the same directory as `annotate_sequences.py` is the simplest and recommended approach.

---

## Troubleshooting

### `ModuleNotFoundError: No module named 'annotation_align_core'`

Possible causes:

- the `.pyd` file is not in the same directory as `annotate_sequences.py`
- the command is being run from an unexpected directory
- the `.pyd` filename was changed
- the binary was not copied from this directory

Confirm that the following file exists:

```text
annotation_align_core.cp313-win_amd64.pyd
```

and place it next to:

```text
annotate_sequences.py
```

Then rerun the import test.

---

### `ImportError: DLL load failed while importing annotation_align_core`

Possible causes include:

- Python is not version 3.13
- Python is 32-bit
- the operating system is not 64-bit Windows
- the binary is damaged or incomplete
- a required Microsoft runtime component is unavailable
- another incompatible file named `annotation_align_core` is being imported first

Check the environment:

```bat
python --version
python -c "import sys; print(sys.maxsize > 2**32)"
python -c "import platform; print(platform.platform()); print(platform.architecture())"
```

Also check which files are present in the working directory:

```bat
dir annotation_align_core*
```

There should normally be only one active prebuilt extension with this module name.

---

### `ImportError: dynamic module does not define module export function`

This usually indicates that:

- the `.pyd` file was renamed incorrectly
- the binary module name and filename no longer match
- an unrelated `.pyd` file is being used

Do not rename:

```text
annotation_align_core.cp313-win_amd64.pyd
```

to an unrelated module name.

The Python import name must remain:

```python
import annotation_align_core
```

---

### The user has Python 3.11 or Python 3.12

This prebuilt binary cannot be used directly.

Use the source distribution instead:

```text
annotation_align_core.pyx
annotation_align_core.c
setup.py
pyproject.toml
```

The extension must be compiled for the user's own Python version.

Compatibility with Python versions other than 3.13 is not claimed unless separately tested.

---

### The user has no compiler and a different Python version

The supported prebuilt binary requires CPython 3.13 on 64-bit Windows.

Without a compiler, the practical option is to install a compatible 64-bit CPython 3.13 environment and use this `.pyd` file.

---

## Binary and source relationship

This file:

```text
annotation_align_core.cp313-win_amd64.pyd
```

is the compiled form of the Cython source:

```text
annotation_align_core.pyx
```

The generated C source is:

```text
annotation_align_core.c
```

The three files have different roles:

| File | Role |
|---|---|
| `annotation_align_core.pyx` | Human-readable Cython source |
| `annotation_align_core.c` | C source generated by Cython |
| `annotation_align_core.cp313-win_amd64.pyd` | Compiled Python extension for CPython 3.13 on 64-bit Windows |

The `.pyd` file is provided for convenience. The `.pyx` and `.c` files remain the primary source materials for reproducibility and rebuilding.

---

## Reproducibility information

The corresponding Cython source was processed using:

```text
Cython 3.2.4
```

The recorded build environment also included:

```text
setuptools 82.0.1
wheel 0.46.3
```

These build tools are not required merely to use this prebuilt `.pyd` file.

The prebuilt binary should be treated as environment-specific. For long-term reproducibility, retain:

- `annotation_align_core.pyx`
- `annotation_align_core.c`
- `setup.py`
- `pyproject.toml`
- `requirements-build.txt`
- this prebuilt `.pyd` file
- the exact Python and operating-system information used for testing

---

## Security and provenance

Only use the `.pyd` file obtained from the official project repository or archived project release.

A `.pyd` file is executable native code. Do not replace it with a binary obtained from an unverified source.

For formal releases, it is recommended to publish a checksum for the binary, for example:

```bat
certutil -hashfile annotation_align_core.cp313-win_amd64.pyd SHA256
```

The resulting SHA-256 value can be recorded in the GitHub Release notes or a checksum file.

---

## License

This binary is distributed under the same license as the source code in the main project repository.

Refer to the repository-level `LICENSE` file for the applicable terms.

---

## Citation

When this binary is used as part of the annotation workflow, cite:

1. the associated research article
2. the archived software release
3. the reference databases used for annotation

Release-specific citation information should be provided in the repository-level `CITATION.cff` file.
