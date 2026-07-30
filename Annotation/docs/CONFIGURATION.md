# Configuration

## Command-line options

| Option | Required | Default | Description |
|---|---:|---:|---|
| `--query-fasta` | yes | — | Uncompressed query FASTA |
| `--references-file` | yes | — | CSV containing `reference_name,file_path` |
| `--output-dir` | yes | — | Output and cache directory |
| `--core-name` | no | query FASTA basename | Output filename prefix |
| `--workers` | no | `16` | Multiprocessing worker count; must be at least 1 |

Show the installed interface:

```bat
python annotate_sequences.py --help
```

## Analysis-defining constants

The following values affect accepted hits or reported results and must not be changed without creating a separately versioned analysis:

```text
MAX_EXT_TOTAL = 3
MAX_MM = 3
UNANNOTATED_TOTAL_SCORE = 10
```

The split-rank ordering is implemented in `build_split_rank_map()` and must be treated as part of the method.

## Performance constants

The following settings control batching, scheduling, caching, or temporary output:

```text
PREFIX_CHECK_LEN = 6
USE_DYNAMIC_BATCH_SIZE = True
PILOT_QUERY_COUNT = 16
MIN_BATCH_SIZE = 40
MAX_BATCH_SIZE = 100000
TASK_ORDER_STRATEGY = weighted_interleave
DOMINANT_FIRST_WAVE = True
KEEP_TEMP_FILES = False
REFERENCE_CACHE_FORMAT_VERSION = 2
```

These settings are intended to preserve results while changing runtime behavior. Nevertheless, formal reproducibility requires retaining the exact script rather than reconstructing settings from documentation.

## Worker count

`--workers` changes parallelism and memory use. It should not change the accepted hit set. Record the worker count because it affects runtime, task allocation, and log output.

Use a lower value when memory is limited. Do not exceed the number of logical processors without benchmarking.

## Cache policy

A cache is reused when the following match:

- cache format version
- absolute reference path
- reference file modification time
- reference file size

A checksum is not used for cache validation. Consequently, use a fresh output directory or delete the cache when:

- reproducing a published analysis
- replacing a reference file
- moving between releases
- uncertain whether file metadata accurately reflects content

## Source and binary configuration

The extension import name is:

```python
annotation_align_core
```

The source files are:

```text
annotation_align_core.pyx
annotation_align_core.c
```

The prebuilt binary for the tested environment is:

```text
annotation_align_core.cp313-win_amd64.pyd
```

Do not mix extension binaries from different Python versions or architectures in the same working directory.
