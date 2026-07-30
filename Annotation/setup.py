#!/usr/bin/env python3
"""
Build configuration for the annotation_align_core extension.

Default:
    Build from annotation_align_core.pyx using Cython.

Build from the generated C source without Cython:
    Windows Command Prompt:
        set ANNOTATION_ALIGN_BUILD_FROM_C=1
        python -m pip install --no-build-isolation .
"""

import os
from pathlib import Path

from setuptools import Extension, setup


PROJECT_ROOT = Path(__file__).resolve().parent

BUILD_FROM_C = (
    os.environ.get("ANNOTATION_ALIGN_BUILD_FROM_C", "")
    .strip()
    .lower()
    in {"1", "true", "yes", "on"}
)

EXTENSION_NAME = "annotation_align_core"

# Preserve the optimization setting used in the original analysis.
if os.name == "nt":
    extra_compile_args = ["/O2"]
else:
    extra_compile_args = ["-O3"]


if BUILD_FROM_C:
    source_file = PROJECT_ROOT / "annotation_align_core.c"

    if not source_file.is_file():
        raise FileNotFoundError(
            "C-source build was requested, but the generated C file "
            f"was not found: {source_file}"
        )

    ext_modules = [
        Extension(
            name=EXTENSION_NAME,
            sources=[str(source_file)],
            extra_compile_args=extra_compile_args,
        )
    ]

else:
    source_file = PROJECT_ROOT / "annotation_align_core.pyx"

    if not source_file.is_file():
        raise FileNotFoundError(
            "Cython build was requested, but the Cython source file "
            f"was not found: {source_file}"
        )

    try:
        from Cython.Build import cythonize
    except ImportError as exc:
        raise RuntimeError(
            "Cython is required for the default build. Install the "
            "build requirements, or set "
            "ANNOTATION_ALIGN_BUILD_FROM_C=1 to build from C."
        ) from exc

    extensions = [
        Extension(
            name=EXTENSION_NAME,
            sources=[str(source_file)],
            extra_compile_args=extra_compile_args,
        )
    ]

    ext_modules = cythonize(
        extensions,
        compiler_directives={
            "language_level": 3,
            "boundscheck": False,
            "wraparound": False,
            "initializedcheck": False,
            "cdivision": True,
        },
        annotate=False,
    )


setup(
    name="sequence-based-smallrna-annotation-core",
    version="2.0.0",
    description=(
        "Compiled core for exact-sequence small RNA reference annotation"
    ),
    python_requires=">=3.11",
    py_modules=["annotate_sequences"],
    ext_modules=ext_modules,
)