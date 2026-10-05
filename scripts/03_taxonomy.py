#!/usr/bin/env python3
"""
03_taxonomy.py
Windows-safe ITS taxonomy step for the metabarcoding workbench.

Supports selecting reference dataset regions: ITS1, ITS2, or full ITS.
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO


PROJECT_ROOT = Path(__file__).resolve().parents[1]
ASV_FASTA = PROJECT_ROOT / "dnabarcoder_input" / "asvs.fasta"
REFERENCE_DIR = PROJECT_ROOT / "references"
RESULTS_DIR = PROJECT_ROOT / "results" / "dnabarcoder"

# Reference dataset mapping for UNITE+INSD 2025
REFERENCE_CONFIG = {
    "its1": {
        "label": "unite2025ITS1",
        "fasta": "unite2025ITS1.fasta",
        "classification": "unite2025ITS1.classification",
        "cutoffs": "unite2025ITS1.unique.cutoffs.best.json",
    },
    "its2": {
        "label": "unite2025ITS2",
        "fasta": "unite2025ITS2.fasta",
        "classification": "unite2025ITS2.classification",
        "cutoffs": "unite2025ITS2.unique.cutoffs.best.json",
    },
    "its": {
        "label": "unite2025ITS",
        "fasta": "unite2025ITS.fasta",
        "classification": "unite2025ITS.classification",
        "cutoffs": "unite2025ITS.unique.cutoffs.best.json",
    },
}

DNABARCODER_CLASSIFY_CANDIDATES = [
    PROJECT_ROOT / "dnabarcoder" / "classification" / "classify.py",
    PROJECT_ROOT / "yeastBarcoder" / "dnabarcoder" / "classification" / "classify.py",
]


def find_classify_py() -> Path | None:
    for path in DNABARCODER_CLASSIFY_CANDIDATES:
        if path.is_file():
            return path
    return None


def check_required_files(region: str) -> dict[str, Path]:
    if not ASV_FASTA.is_file():
        raise FileNotFoundError(
            f"ASV FASTA not found:\n{ASV_FASTA}\n"
            "Run the DADA2 step first."
        )

    region_config = REFERENCE_CONFIG[region]
    refs: dict[str, Path] = {}
    missing: list[str] = []

    for key in ["fasta", "classification", "cutoffs"]:
        filename = region_config[key]
        path = REFERENCE_DIR / filename
        if path.is_file():
            refs[key] = path
        else:
            missing.append(filename)

    if missing:
        raise FileNotFoundError(
            f"Required reference files for region '{region}' are missing from "
            f"{REFERENCE_DIR}:\n  - " + "\n  - ".join(missing)
        )

    return refs


def run_command(args: list[str], cwd: Path = PROJECT_ROOT) -> None:
    print("\n$ " + " ".join(f'"{x}"' if " " in str(x) else str(x) for x in args))
    subprocess.run(args, cwd=cwd, check=True)


def fasta_records(path: Path):
    records = list(SeqIO.parse(path, "fasta"))
    if not records:
        raise RuntimeError(f"No FASTA records found in {path}")
    return records


def make_indexed_query(query: Path, destination: Path) -> tuple[dict[int, str], list]:
    records = fasta_records(query)
    id_map: dict[int, str] = {}

    with destination.open("w", encoding="utf-8", newline="\n") as out:
        for i, record in enumerate(records):
            id_map[i] = record.id
            out.write(f">{i}|{record.id}\n")
            sequence = str(record.seq).strip()
            for start in range(0, len(sequence), 80):
                out.write(sequence[start:start + 80] + "\n")

    return id_map, records


def get_base(path: Path) -> Path:
    return path.with_suffix("")


def bestmatch_search(
    query: Path,
    reference: Path,
    output: Path,
    min_coverage: int = 50,
    ncpus: int = 4,
) -> None:
    indexed_query = output.with_suffix(".indexed.fasta")
    blast_output = output.with_suffix(".blastoutput")

    id_map, _ = make_indexed_query(query, indexed_query)

    db = get_base(reference).with_name(get_base(reference).name + ".blastdb")

    db_nsq = Path(str(db) + ".nsq")
    if not db_nsq.exists():
        run_command([
            "makeblastdb",
            "-in", str(reference),
            "-dbtype", "nucl",
            "-out", str(db),
        ])
    else:
        print(f"Using existing BLAST database: {db}")

    blast_args = [
        "blastn",
        "-query", str(indexed_query),
        "-db", str(db),
        "-task", "blastn-short",
        "-outfmt", "6",
        "-out", str(blast_output),
        "-num_threads", str(ncpus),
    ]
    run_command(blast_args)

    best_score = {i: 0.0 for i in id_map}
    best_sim = {i: 0.0 for i in id_map}
    best_coverage = {i: 0 for i in id_map}
    best_ref = {i: "" for i in id_map}

    with blast_output.open("r", encoding="utf-8") as handle:
        for line in handle:
            words = line.rstrip("\n").split("\t")
            if len(words) < 8:
                continue

            try:
                query_id = int(words[0].split("|", 1)[0])
                qstart = int(words[6])
                qend = int(words[7])
                identity = float(words[2])
            except (ValueError, IndexError):
                continue

            if query_id not in id_map:
                continue

            similarity = identity / 100.0
            coverage = abs(qend - qstart)
            score = similarity

            if coverage < min_coverage:
                score = (score * coverage) / min_coverage

            if (
                score > best_score[query_id]
                or (
                    score == best_score[query_id]
                    and coverage > best_coverage[query_id]
                )
            ):
                best_score[query_id] = score
                best_sim[query_id] = similarity
                best_coverage[query_id] = coverage
                best_ref[query_id] = words[1]

    output.parent.mkdir(parents=True, exist_ok=True)

    with output.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "ID",
            "ReferenceID",
            "BLAST score",
            "BLAST sim",
            "BLAST coverage",
        ])
        for i in range(len(id_map)):
            writer.writerow([
                id_map[i],
                best_ref[i],
                best_score[i],
                best_sim[i],
                best_coverage[i],
            ])

    if indexed_query.exists():
        indexed_query.unlink()

    if blast_output.exists():
        blast_output.unlink()

    matched = sum(1 for value in best_ref.values() if value)
    print(f"Best matches written: {matched}/{len(id_map)}")


def run_classification(
    classify_py: Path,
    bestmatch: Path,
    classification: Path,
    cutoffs: Path,
    output: Path,
) -> Path:
    output_dir = output.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    run_command(
        [
            sys.executable,
            str(classify_py),
            "-i", bestmatch.name,
            "-c", str(classification),
            "-cutoffs", str(cutoffs),
            "-o", ".",
        ],
        cwd=output_dir,
    )

    candidates = sorted(
        output_dir.glob(f"{bestmatch.stem}*.classified"),
        key=lambda p: p.stat().st_mtime,
        reverse=True,
    )

    if not candidates:
        candidates = sorted(
            output_dir.glob("*.classified"),
            key=lambda p: p.stat().st_mtime,
            reverse=True,
        )

    if not candidates:
        raise RuntimeError(
            "dnabarcoder classification completed but no .classified "
            f"file was found in:\n{output_dir}"
        )

    generated = candidates[0]

    if generated.resolve() != output.resolve():
        if output.exists():
            output.unlink()
        generated.replace(output)

    return output


def process_reference(
    label: str,
    fasta: Path,
    classification: Path,
    cutoffs: Path,
    classify_py: Path,
    ncpus: int,
) -> Path:
    print("\n" + "=" * 70)
    print(f"{label}: BLAST SEARCH")
    print("=" * 70)

    bestmatch = RESULTS_DIR / f"{ASV_FASTA.stem}.{fasta.stem}_BLAST.bestmatch"

    bestmatch_search(
        query=ASV_FASTA,
        reference=fasta,
        output=bestmatch,
        min_coverage=50,
        ncpus=ncpus,
    )

    print("\n" + "=" * 70)
    print(f"{label}: DNABARCODER CLASSIFICATION")
    print("=" * 70)

    classified = RESULTS_DIR / f"{ASV_FASTA.stem}.{fasta.stem}_BLAST.classified"

    run_classification(
        classify_py=classify_py,
        bestmatch=bestmatch,
        classification=classification,
        cutoffs=cutoffs,
        output=classified,
    )

    return classified


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Windows-safe ITS dnabarcoder taxonomy workflow."
    )
    parser.add_argument(
        "--region", "-r",
        type=str.lower,
        choices=["its1", "its2", "its"],
        default="its2",
        help="Select reference marker region: 'its1', 'its2', or 'its' (default: its2).",
    )
    parser.add_argument(
        "--ncpus",
        type=int,
        default=4,
        help="BLAST CPU threads (default: 4).",
    )
    args = parser.parse_args()

    if args.ncpus < 1:
        parser.error("--ncpus must be at least 1")

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    ref_info = REFERENCE_CONFIG[args.region]
    ref_label = ref_info["label"]

    print("=" * 70)
    print(f"03 — DNA BARCODER TAXONOMY ({args.region.upper()})")
    print("=" * 70)
    print(f"Project:   {PROJECT_ROOT}")
    print(f"ASVs:      {ASV_FASTA}")
    print(f"Reference: {ref_label}")

    refs = check_required_files(args.region)

    classify_py = find_classify_py()
    if classify_py is None:
        raise FileNotFoundError(
            "dnabarcoder classification.py was not found. Expected one of:\n"
            + "\n".join(f"  {p}" for p in DNABARCODER_CLASSIFY_CANDIDATES)
        )

    if shutil.which("blastn") is None:
        raise FileNotFoundError(
            "blastn was not found on PATH. In PowerShell, confirm with:\n"
            "  blastn -version"
        )

    if shutil.which("makeblastdb") is None:
        raise FileNotFoundError(
            "makeblastdb was not found on PATH. In PowerShell, confirm with:\n"
            "  makeblastdb -help"
        )

    print(f"\nClassification code:\n  {classify_py}")

    classified = process_reference(
        ref_label,
        refs["fasta"],
        refs["classification"],
        refs["cutoffs"],
        classify_py,
        args.ncpus,
    )

    print("\n" + "=" * 70)
    print("TAXONOMY COMPLETE")
    print("=" * 70)
    print(f"Classification: {classified}")
    print("\nFiles in the dnabarcoder results folder:")
    for path in sorted(RESULTS_DIR.iterdir()):
        if path.is_file():
            print(f"  {path.name}")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except subprocess.CalledProcessError as exc:
        print(f"\nERROR: external command failed with exit code {exc.returncode}")
        raise SystemExit(exc.returncode or 1)
    except KeyboardInterrupt:
        print("\nInterrupted.")
        raise SystemExit(130)
    except Exception as exc:
        print(f"\nERROR: {exc}")
        raise SystemExit(1)