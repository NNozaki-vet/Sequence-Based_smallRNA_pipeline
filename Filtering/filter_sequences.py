#!/usr/bin/env python3
r"""
Exact-sequence filtering for small RNA FASTA files.

This script preserves the filtering procedure used in the original analysis:

1. Count every exact sequence in each sample FASTA file.
2. Calculate CPM using the total number of reads in each sample.
3. Within each group, retain a sequence when:
       raw count >= 2
       CPM >= 1
   in at least max(2, ceiling(n / 2)) samples, where n is the number
   of samples in that group.
4. Retain the union of sequences passing the criterion in at least one group.
5. Export the raw-count matrix, filtered FASTA, group-pass summary,
   sample metadata, and group-specific thresholds.

Input directories and the output directory are supplied at runtime. The
filtering criteria and computational procedure are unchanged.

Example:
    python filter_sequences.py ^
        --groups-file "D:\project\group_directories.csv" ^
        --output-dir "D:\project\filtered_sequences" ^
        --threads 8

The groups CSV must contain the following header:

    group,folder

Example rows:

    Thymus,D:\project\Trimmed_fa\Thymus
    Brain,D:\project\Trimmed_fa\Brain
"""

import os
import re
import csv
import gzip
import math
import argparse
from pathlib import Path
from collections import Counter, defaultdict
from multiprocessing import Pool, cpu_count

# ============================================================
# Fixed analysis settings
# ============================================================

THREADS = 8

RAW_COUNT_MIN = 2
CPM_MIN = 1.0

# If FASTA contains U, convert to T for DNA-style sequence consistency.
CONVERT_U_TO_T = True

FASTA_EXTENSIONS = (
    ".fa", ".fasta", ".fna",
    ".fa.gz", ".fasta.gz", ".fna.gz"
)


# ============================================================
# Input configuration
# ============================================================

def load_group_dirs(groups_file):
    """
    Read group names and FASTA directories from a CSV file.

    Required columns:
        group
        folder

    Relative folder paths are interpreted relative to the directory
    containing the groups CSV file. Row order is preserved and is used
    as the group order in output files.
    """
    groups_file = Path(groups_file)

    if not groups_file.exists():
        raise FileNotFoundError(
            f"Groups file not found: {groups_file}"
        )

    group_dirs = []

    with open(groups_file, "r", newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)

        if reader.fieldnames is None:
            raise ValueError(
                f"The groups file has no header: {groups_file}"
            )

        cleaned_fieldnames = [
            str(name).replace("\ufeff", "").strip()
            for name in reader.fieldnames
        ]
        reader.fieldnames = cleaned_fieldnames

        required_columns = {"group", "folder"}
        missing_columns = required_columns.difference(reader.fieldnames)

        if missing_columns:
            raise ValueError(
                "The groups file is missing required columns: "
                + ", ".join(sorted(missing_columns))
            )

        for row_number, row in enumerate(reader, start=2):
            group = str(row["group"]).strip()
            folder_text = str(row["folder"]).strip()

            if not group and not folder_text:
                continue

            if not group:
                raise ValueError(
                    f"Empty group name at row {row_number}: {groups_file}"
                )

            if not folder_text:
                raise ValueError(
                    f"Empty folder path at row {row_number}: {groups_file}"
                )

            folder = Path(folder_text).expanduser()

            if not folder.is_absolute():
                folder = groups_file.parent / folder

            group_dirs.append((group, folder))

    if not group_dirs:
        raise ValueError(
            f"No group directories were found in: {groups_file}"
        )

    group_names = [group for group, _folder in group_dirs]

    if len(group_names) != len(set(group_names)):
        duplicated = sorted({
            group
            for group in group_names
            if group_names.count(group) > 1
        })
        raise ValueError(
            "Duplicate group names were found in the groups file: "
            + ", ".join(duplicated)
        )

    return group_dirs


# ============================================================
# Functions
# ============================================================

def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def sample_name_from_file(path):
    name = Path(path).name

    for suf in [".fasta.gz", ".fna.gz", ".fa.gz", ".fasta", ".fna", ".fa"]:
        if name.endswith(suf):
            name = name[:-len(suf)]
            break

    for tag in [".sra_trimmed", "_sra_trimmed", ".trimmed", "_trimmed"]:
        if name.endswith(tag):
            name = name[:-len(tag)]

    return name


def sanitize_filename(s):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", s)


def iter_fasta_sequences(path):
    """
    Read FASTA records safely even when sequences are line-wrapped.
    """
    seq_parts = []

    with open_text(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if seq_parts:
                    seq = "".join(seq_parts).upper()
                    if CONVERT_U_TO_T:
                        seq = seq.replace("U", "T")
                    yield seq
                    seq_parts = []
            else:
                seq_parts.append(line)

        if seq_parts:
            seq = "".join(seq_parts).upper()
            if CONVERT_U_TO_T:
                seq = seq.replace("U", "T")
            yield seq


def count_one_sample(task):
    sample_idx = task["sample_idx"]
    group = task["group"]
    sample_name = task["sample_name"]
    fasta_path = Path(task["fasta_path"])
    per_sample_dir = Path(task["per_sample_dir"])

    counter = Counter()
    total_reads = 0

    for seq in iter_fasta_sequences(fasta_path):
        if not seq:
            continue
        counter[seq] += 1
        total_reads += 1

    out_tsv = (
        per_sample_dir
        / f"{sample_idx:04d}__{sanitize_filename(group)}__"
          f"{sanitize_filename(sample_name)}.seq_counts.tsv"
    )

    with open(out_tsv, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["sequence", "count"])
        for seq, count in counter.items():
            w.writerow([seq, count])

    return {
        "sample_idx": sample_idx,
        "group": group,
        "sample_name": sample_name,
        "fasta_path": str(fasta_path),
        "count_tsv": str(out_tsv),
        "total_reads": total_reads,
        "unique_sequences": len(counter),
    }


def collect_samples(group_dirs, per_sample_dir):
    samples = []
    sample_idx = 0

    for group, folder in group_dirs:
        if not folder.exists():
            print(f"[WARNING] Directory not found: {folder}")
            continue

        files = []
        for ext in FASTA_EXTENSIONS:
            files.extend(folder.glob(f"*{ext}"))

        files = sorted(set(files), key=lambda p: p.name)

        if not files:
            print(f"[WARNING] No FASTA files in: {folder}")
            continue

        for fp in files:
            sample_name = sample_name_from_file(fp)

            samples.append({
                "sample_idx": sample_idx,
                "group": group,
                "sample_name": sample_name,
                "fasta_path": str(fp),
                "per_sample_dir": str(per_sample_dir),
            })
            sample_idx += 1

    return samples


def min_required_samples(n):
    return max(2, math.ceil(n / 2))


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Filter exact small RNA sequences by raw count and CPM "
            "within user-defined sample groups."
        )
    )
    parser.add_argument(
        "--groups-file",
        required=True,
        help=(
            "CSV file containing the columns 'group' and 'folder'. "
            "Each folder must contain sample FASTA files for one group."
        ),
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory in which all output files will be written.",
    )
    parser.add_argument(
        "--output-prefix",
        default="filtered_sequences",
        help=(
            "Prefix used for the count table, FASTA, and pass-group "
            "summary filenames. Default: filtered_sequences"
        ),
    )
    parser.add_argument("--threads", type=int, default=THREADS)
    parser.add_argument(
        "--force-recount",
        action="store_true",
        help=(
            "Recount FASTA files even if per-sample count files "
            "already exist."
        ),
    )
    args = parser.parse_args()

    group_dirs = load_group_dirs(args.groups_file)

    outdir = Path(args.output_dir)
    per_sample_dir = outdir / "per_sample_sequence_counts"

    outdir.mkdir(parents=True, exist_ok=True)
    per_sample_dir.mkdir(parents=True, exist_ok=True)

    samples = collect_samples(
        group_dirs=group_dirs,
        per_sample_dir=per_sample_dir,
    )

    if not samples:
        raise RuntimeError("No FASTA samples were found.")

    print("============================================================")
    print("Exact-sequence filtering")
    print("Condition:")
    print(f"  raw count >= {RAW_COUNT_MIN}")
    print(f"  CPM >= {CPM_MIN}")
    print("  in max(2, ceiling(n/2)) samples within each group")
    print("Output:")
    print(f"  {outdir}")
    print("============================================================")
    print(f"Number of samples: {len(samples)}")
    print(f"Threads: {args.threads}")
    print("")

    # --------------------------------------------------------
    # 1. Count sequences per sample
    # --------------------------------------------------------

    results_by_idx = {}

    tasks_to_count = []
    for s in samples:
        expected_tsv = (
            per_sample_dir
            / f"{s['sample_idx']:04d}__{sanitize_filename(s['group'])}__"
              f"{sanitize_filename(s['sample_name'])}.seq_counts.tsv"
        )
        if expected_tsv.exists() and not args.force_recount:
            # Existing count files are retained in the same manner as
            # in the original analysis script.
            pass
        tasks_to_count.append(s)

    print("[Step 1] Counting sequences in each FASTA file...")
    with Pool(processes=args.threads) as pool:
        for res in pool.imap_unordered(count_one_sample, tasks_to_count):
            results_by_idx[res["sample_idx"]] = res
            print(
                f"  counted: {res['sample_name']} "
                f"({res['group']}) "
                f"reads={res['total_reads']} unique={res['unique_sequences']}"
            )

    sample_results = [results_by_idx[i] for i in sorted(results_by_idx)]

    # --------------------------------------------------------
    # 2. Write sample metadata and group thresholds
    # --------------------------------------------------------

    sample_metadata_csv = outdir / "sample_metadata.csv"

    group_to_sample_indices = defaultdict(list)
    for r in sample_results:
        group_to_sample_indices[r["group"]].append(r["sample_idx"])

    group_thresholds = {}
    for group, idxs in group_to_sample_indices.items():
        n = len(idxs)
        group_thresholds[group] = min_required_samples(n)

    with open(sample_metadata_csv, "w", newline="") as out:
        w = csv.writer(out)
        w.writerow([
            "sample_idx",
            "group",
            "sample_name",
            "fasta_path",
            "count_tsv",
            "total_reads",
            "unique_sequences",
            "group_n_samples",
            "group_required_samples",
            "low_depth_lt_1M",
        ])

        for r in sample_results:
            group = r["group"]
            w.writerow([
                r["sample_idx"],
                group,
                r["sample_name"],
                r["fasta_path"],
                r["count_tsv"],
                r["total_reads"],
                r["unique_sequences"],
                len(group_to_sample_indices[group]),
                group_thresholds[group],
                "YES" if r["total_reads"] < 1_000_000 else "NO",
            ])

    group_summary_csv = outdir / "group_filter_thresholds.csv"
    with open(group_summary_csv, "w", newline="") as out:
        w = csv.writer(out)
        w.writerow(["group", "n_samples", "required_samples"])
        for group in [
            group
            for group, _folder in group_dirs
            if group in group_to_sample_indices
        ]:
            w.writerow([
                group,
                len(group_to_sample_indices[group]),
                group_thresholds[group],
            ])

    print("")
    print("[Step 2] Sample metadata written:")
    print(f"  {sample_metadata_csv}")
    print(f"  {group_summary_csv}")

    # --------------------------------------------------------
    # 3. Identify sequences passing the group-level filter
    # --------------------------------------------------------

    print("")
    print("[Step 3] Finding sequences passing group-level filter...")

    pass_count_by_seq_group = defaultdict(lambda: defaultdict(int))

    total_reads_by_idx = {
        r["sample_idx"]: r["total_reads"] for r in sample_results
    }
    group_by_idx = {
        r["sample_idx"]: r["group"] for r in sample_results
    }

    for r in sample_results:
        sample_idx = r["sample_idx"]
        group = group_by_idx[sample_idx]
        libsize = total_reads_by_idx[sample_idx]

        if libsize <= 0:
            continue

        with open(r["count_tsv"], "r", newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                seq = row["sequence"]
                count = int(row["count"])
                cpm = count / libsize * 1_000_000

                if count >= RAW_COUNT_MIN and cpm >= CPM_MIN:
                    pass_count_by_seq_group[seq][group] += 1

    retained = set()
    pass_groups_by_seq = {}

    group_order = {
        group: index
        for index, (group, _folder) in enumerate(group_dirs)
    }

    for seq, gdict in pass_count_by_seq_group.items():
        passed_groups = []
        for group, pass_n in gdict.items():
            if pass_n >= group_thresholds[group]:
                passed_groups.append(group)

        if passed_groups:
            retained.add(seq)
            pass_groups_by_seq[seq] = sorted(
                passed_groups,
                key=lambda group: group_order[group],
            )

    print(f"  retained sequences: {len(retained)}")

    # --------------------------------------------------------
    # 4. Build the count matrix for retained sequences
    # --------------------------------------------------------

    print("")
    print("[Step 4] Building count table for retained sequences...")

    n_samples = len(sample_results)
    counts_by_seq = {
        seq: [0] * n_samples for seq in retained
    }
    total_count_by_seq = defaultdict(int)

    for r in sample_results:
        sample_idx = r["sample_idx"]

        with open(r["count_tsv"], "r", newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                seq = row["sequence"]
                if seq not in retained:
                    continue

                count = int(row["count"])
                counts_by_seq[seq][sample_idx] = count
                total_count_by_seq[seq] += count

    ordered_sequences = sorted(
        retained,
        key=lambda s: (-total_count_by_seq[s], len(s), s)
    )

    # --------------------------------------------------------
    # 5. Write the output count table
    # --------------------------------------------------------

    count_table_csv = (
        outdir
        / f"{args.output_prefix}_count_table.csv"
    )

    sample_names = [r["sample_name"] for r in sample_results]

    # Add sample indices when duplicate sample names are detected.
    if len(sample_names) != len(set(sample_names)):
        print(
            "[WARNING] Duplicate sample names detected. "
            "Adding sample_idx to column names."
        )
        sample_names = [
            f"{r['sample_name']}__idx{r['sample_idx']}"
            for r in sample_results
        ]

    with open(count_table_csv, "w", newline="") as out:
        w = csv.writer(out)
        w.writerow([
            "sequence",
            "length",
            "total_count",
            "n_pass_groups",
            "pass_groups",
            *sample_names
        ])

        for seq in ordered_sequences:
            w.writerow([
                seq,
                len(seq),
                total_count_by_seq[seq],
                len(pass_groups_by_seq[seq]),
                ";".join(pass_groups_by_seq[seq]),
                *counts_by_seq[seq]
            ])

    # --------------------------------------------------------
    # 6. Write filtered sequences in FASTA format
    # --------------------------------------------------------

    fasta_out = outdir / f"{args.output_prefix}.fa"

    with open(fasta_out, "w") as out:
        for i, seq in enumerate(ordered_sequences, start=1):
            out.write(
                f">seq_{i:09d} "
                f"length={len(seq)} "
                f"total_count={total_count_by_seq[seq]} "
                f"n_pass_groups={len(pass_groups_by_seq[seq])} "
                f"pass_groups={';'.join(pass_groups_by_seq[seq])}\n"
            )
            out.write(seq + "\n")

    # --------------------------------------------------------
    # 7. Write the pass-group summary
    # --------------------------------------------------------

    pass_summary_csv = (
        outdir
        / f"{args.output_prefix}_pass_groups.csv"
    )

    with open(pass_summary_csv, "w", newline="") as out:
        w = csv.writer(out)
        w.writerow([
            "sequence",
            "length",
            "total_count",
            "n_pass_groups",
            "pass_groups",
        ])

        for seq in ordered_sequences:
            w.writerow([
                seq,
                len(seq),
                total_count_by_seq[seq],
                len(pass_groups_by_seq[seq]),
                ";".join(pass_groups_by_seq[seq]),
            ])

    print("")
    print("============================================================")
    print("Done.")
    print("Output files:")
    print(f"  Count table : {count_table_csv}")
    print(f"  FASTA       : {fasta_out}")
    print(f"  Pass groups : {pass_summary_csv}")
    print(f"  Metadata    : {sample_metadata_csv}")
    print(f"  Thresholds  : {group_summary_csv}")
    print("============================================================")


if __name__ == "__main__":
    main()
