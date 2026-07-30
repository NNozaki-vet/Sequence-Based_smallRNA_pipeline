# Troubleshooting

## `ModuleNotFoundError: No module named 'annotation_align_core'`

The compiled extension is unavailable to the active Python interpreter.

Check:

```bat
python -c "import sys; print(sys.executable)"
dir annotation_align_core*
```

For the prebuilt binary, place:

```text
annotation_align_core.cp313-win_amd64.pyd
```

beside `annotate_sequences.py` and use 64-bit CPython 3.13 on Windows.

## `ImportError: DLL load failed`

The `.pyd` is incompatible or a required runtime component is unavailable.

Verify:

```bat
python --version
python -c "import platform; print(platform.python_implementation()); print(platform.architecture())"
```

The prebuilt binary requires CPython 3.13 and `64bit`. Rebuild from source for other environments.

## Reference CSV errors

Required headers are exactly:

```csv
reference_name,file_path
```

Common causes:

- capitalization differs
- duplicate `reference_name`
- empty name or path
- relative path interpreted from the wrong directory

Use absolute paths temporarily to diagnose path resolution.

## Query or reference FASTA not found

Confirm the exact path:

```bat
dir "D:\path\file.fa"
```

The workflow does not read gzip-compressed FASTA files.

## Empty or malformed FASTA

Reference FASTA files raise errors for an empty header, empty sequence, or sequence before the first header. Query files should also be valid FASTA even though query parsing is less strict.

Validate all files before a production run.

## Unexpected `Unannotated` results

Check:

- query and reference orientation
- DNA/RNA notation; the script converts `U` to `T`
- sequence alphabet
- reference version and species
- maximum of three internal mismatches
- maximum total terminal query extension of three
- absence of internal indel support

## Results differ between reruns

Use:

1. identical source files
2. identical query and reference files
3. identical command
4. a fresh output directory
5. the same tested Python environment
6. the same compiled core

Compare SHA-256 checksums. Do not rely only on filenames or modification dates.

## Stale reference cache suspected

Delete:

```text
<output_dir>\reference_cache\
```

or use a new output directory, then rerun.

## Run stops or memory is exhausted

Reduce workers:

```bat
--workers 4
```

The reference cache and per-worker query bundles can require substantial memory. Record the changed worker count in the run log.

## Output CSV opens incorrectly in spreadsheet software

The CSV is UTF-8 with a byte-order mark. Import it explicitly as UTF-8 and comma-delimited. Do not resave the primary output before checksum calculation.
