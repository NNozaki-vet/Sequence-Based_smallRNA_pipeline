#!/usr/bin/env python3
"""
Run the annotation integration test without third-party test packages.

Execute from any working directory:

    python tests/test_annotation.py
"""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path


TESTS_DIR = Path(__file__).resolve().parent
ANNOTATION_DIR = TESTS_DIR.parent

ANNOTATION_SCRIPT = ANNOTATION_DIR / "annotate_sequences.py"
QUERY_FASTA = TESTS_DIR / "data" / "queries.fa"
REFERENCES_CSV = TESTS_DIR / "data" / "references.csv"
EXPECTED_DIR = TESTS_DIR / "expected"

CORE_NAME = "test_sequences"
WORKERS = "1"


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def validate_required_files() -> list[str]:
    required = [
        ANNOTATION_SCRIPT,
        QUERY_FASTA,
        REFERENCES_CSV,
        TESTS_DIR / "data" / "references" / "Reference_A.fa",
        TESTS_DIR / "data" / "references" / "Reference_B.fa",
        TESTS_DIR / "data" / "references" / "Reference_C.fa",
    ]

    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        fail("Required test files are missing:\n  " + "\n  ".join(missing))

    reference_names: list[str] = []

    with REFERENCES_CSV.open(
        "r",
        newline="",
        encoding="utf-8-sig",
    ) as handle:
        reader = csv.DictReader(handle)

        if reader.fieldnames is None:
            fail(f"Reference CSV has no header: {REFERENCES_CSV}")

        required_columns = {"reference_name", "file_path"}
        if not required_columns.issubset(reader.fieldnames):
            fail(
                "Reference CSV must contain the columns "
                "'reference_name' and 'file_path'."
            )

        for row_number, row in enumerate(reader, start=2):
            reference_name = str(row["reference_name"]).strip()
            file_path = str(row["file_path"]).strip()

            if not reference_name and not file_path:
                continue

            if not reference_name or not file_path:
                fail(f"Incomplete reference definition at row {row_number}.")

            reference_names.append(reference_name)

    if not reference_names:
        fail("No reference definitions were found.")

    if len(reference_names) != len(set(reference_names)):
        fail("Duplicate reference_name values were found.")

    expected_files = [
        EXPECTED_DIR / f"{CORE_NAME}_{name}_minimumMiss.csv"
        for name in reference_names
    ]

    missing_expected = [
        str(path)
        for path in expected_files
        if not path.is_file()
    ]

    if missing_expected:
        fail(
            "Expected CSV files are missing. Generate and manually review "
            "them before running the test:\n  "
            + "\n  ".join(missing_expected)
        )

    return reference_names


def read_csv(path: Path) -> tuple[tuple[str, ...], Counter[tuple[str, ...]]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        header = next(reader, None)

        if header is None:
            fail(f"CSV file is empty: {path}")

        rows = Counter(tuple(row) for row in reader)

    return tuple(header), rows


def compare_csv(expected_path: Path, actual_path: Path) -> None:
    expected_header, expected_rows = read_csv(expected_path)
    actual_header, actual_rows = read_csv(actual_path)

    if actual_header != expected_header:
        fail(
            f"Header mismatch: {actual_path.name}\n"
            f"Expected: {expected_header}\n"
            f"Actual:   {actual_header}"
        )

    if actual_rows != expected_rows:
        missing_rows = expected_rows - actual_rows
        unexpected_rows = actual_rows - expected_rows

        details = [f"CSV content mismatch: {actual_path.name}"]

        if missing_rows:
            details.append("Missing rows:")
            for row, count in missing_rows.items():
                details.append(f"  {count} x {row}")

        if unexpected_rows:
            details.append("Unexpected rows:")
            for row, count in unexpected_rows.items():
                details.append(f"  {count} x {row}")

        fail("\n".join(details))


def run_test() -> None:
    reference_names = validate_required_files()

    with tempfile.TemporaryDirectory(
        prefix="annotation_integration_test_"
    ) as temporary_directory:
        output_dir = Path(temporary_directory) / "output"

        command = [
            sys.executable,
            str(ANNOTATION_SCRIPT),
            "--query-fasta",
            str(QUERY_FASTA),
            "--references-file",
            str(REFERENCES_CSV),
            "--output-dir",
            str(output_dir),
            "--core-name",
            CORE_NAME,
            "--workers",
            WORKERS,
        ]

        completed = subprocess.run(
            command,
            cwd=ANNOTATION_DIR,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )

        if completed.returncode != 0:
            fail(
                "Annotation command failed.\n"
                f"Command: {' '.join(command)}\n"
                "Output:\n"
                f"{completed.stdout}"
            )

        for reference_name in reference_names:
            filename = (
                f"{CORE_NAME}_{reference_name}_minimumMiss.csv"
            )
            expected_path = EXPECTED_DIR / filename
            actual_path = output_dir / filename

            if not actual_path.is_file():
                fail(
                    f"Expected output was not created: {actual_path}\n"
                    f"Annotation output:\n{completed.stdout}"
                )

            compare_csv(expected_path, actual_path)

    print("Annotation integration test: PASSED")


if __name__ == "__main__":
    run_test()
