#!/usr/bin/env python3
"""Download paired-end FASTQ files for any BioProject/study accession from ENA.

Usage examples:
    python download_fastq.py PRJNA934949
    python download_fastq.py PRJNA356769 --library-name ITS2
    python download_fastq.py PRJNA934949 --fastq-dir D:/fastq --dry-run

Only Python 3 (standard library) is required. Paired files are saved as
<run>_1.fastq.gz and <run>_2.fastq.gz.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ENA_API = "https://www.ebi.ac.uk/ena/portal/api/search"
ACCESSION_RE = re.compile(r"^(PRJ(NA|EB|DB)\d+|[SED]RP\d+)$", re.IGNORECASE)
FIELDS = [
    "run_accession", "study_accession", "experiment_accession",
    "sample_accession", "library_name", "library_strategy",
    "library_layout", "fastq_ftp",
]
CHUNK = 1024 * 1024
TIMEOUT = 60
RETRIES = 3
USER_AGENT = "metabarcoding-workbench/1.0"


def parse_args():
    p = argparse.ArgumentParser(
        description="Download paired-end FASTQ files for a BioProject from ENA."
    )
    p.add_argument("bioproject",
                   help="BioProject/study accession, e.g. PRJNA934949 (also PRJEB..., PRJDB..., SRP...)")
    p.add_argument("--fastq-dir", type=Path, default=PROJECT_ROOT / "data" / "fastq",
                   help="Where to save FASTQ files (default: <project>/data/fastq)")
    p.add_argument("--results-dir", type=Path, default=PROJECT_ROOT / "results",
                   help="Where to save the manifest (default: <project>/results)")
    p.add_argument("--library-name", default=None,
                   help="Optional: only keep runs whose library_name contains this text (case-insensitive), e.g. ITS2")
    p.add_argument("--library-strategy", default=None,
                   help="Optional: only keep runs with this library_strategy, e.g. AMPLICON")
    p.add_argument("--dry-run", action="store_true",
                   help="Write the manifest and list selected runs, but do not download")
    p.add_argument("--limit", type=int, default=None,
                   help="Optional: maximum number of samples/runs to download")
    
    return p.parse_args()


def get_bytes(url):
    last = None
    for attempt in range(1, RETRIES + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
                return r.read()
        except Exception as e:
            last = e
            print(f"Request failed ({attempt}/{RETRIES}): {e}")
            if attempt < RETRIES:
                time.sleep(2 * attempt)
    raise RuntimeError(last)


def pair_urls(row):
    """Return (R1_url, R2_url) chosen by filename, or None if either is missing.

    ENA may list an extra combined/unpaired file (RUN.fastq.gz) before
    RUN_1.fastq.gz and RUN_2.fastq.gz, so we can't just take the first two.
    """
    urls = [u.strip() for u in (row.get("fastq_ftp") or "").split(";") if u.strip()]
    r1 = next((u for u in urls if u.endswith("_1.fastq.gz")), None)
    r2 = next((u for u in urls if u.endswith("_2.fastq.gz")), None)
    return (to_url(r1), to_url(r2)) if r1 and r2 else None


def metadata(accession, library_name=None, library_strategy=None):
    params = {
        "result": "read_run",
        "query": f'study_accession="{accession}"',
        "fields": ",".join(FIELDS),
        "format": "tsv",
    }
    url = ENA_API + "?" + urllib.parse.urlencode(params)
    rows = list(csv.DictReader(io.StringIO(get_bytes(url).decode()), delimiter="\t"))
    if not rows:
        raise RuntimeError(f"ENA returned no runs for {accession}.")

    selected = []
    for row in rows:
        layout = (row.get("library_layout") or "").strip().upper()
        name = (row.get("library_name") or "").strip().lower()
        strategy = (row.get("library_strategy") or "").strip().upper()

        if layout != "PAIRED" or pair_urls(row) is None:
            continue
        if library_name and library_name.lower() not in name:
            continue
        if library_strategy and library_strategy.upper() != strategy:
            continue
        selected.append(row)
    return rows, selected


def to_url(path):
    path = path.strip()
    if path.startswith("ftp://"):
        path = path[6:]
    if path.startswith("http://") or path.startswith("https://"):
        return path
    return "https://" + path.lstrip("/")


def valid_fastq(p):
    if not p.exists() or p.stat().st_size == 0:
        return False
    try:
        with gzip.open(p, "rt", encoding="utf-8", errors="replace") as f:
            rec, n = [], 0
            for line in f:
                rec.append(line.rstrip("\r\n"))
                if len(rec) == 4:
                    if (not rec[0].startswith("@") or not rec[2].startswith("+")
                            or len(rec[1]) != len(rec[3])):
                        return False
                    n += 1
                    rec = []
            return n > 0 and not rec
    except (OSError, EOFError, gzip.BadGzipFile):
        return False


def download(u, dest):
    part = dest.with_name(dest.name + ".part")
    if part.exists():
        part.unlink()
    last = None
    for attempt in range(1, RETRIES + 1):
        try:
            req = urllib.request.Request(u, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=TIMEOUT) as r, open(part, "wb") as out:
                total = int(r.headers.get("Content-Length") or 0)
                done = 0
                while True:
                    chunk = r.read(CHUNK)
                    if not chunk:
                        break
                    out.write(chunk)
                    done += len(chunk)
                    if total:
                        print(f"\r    {dest.name}: {done / total * 100:6.2f}%", end="")
            print()
            part.replace(dest)
            return
        except Exception as e:
            last = e
            print(f"\n    Download failed ({attempt}/{RETRIES}): {e}")
            if part.exists():
                part.unlink()
            if attempt < RETRIES:
                time.sleep(2 * attempt)
    raise RuntimeError(last)


def main():
    args = parse_args()
    accession = args.bioproject.strip().upper()
    if not ACCESSION_RE.match(accession):
        print(f"ERROR: '{args.bioproject}' doesn't look like a valid accession "
              "(expected e.g. PRJNA934949, PRJEB12345, PRJDB1234, SRP123456).")
        return 2

    fastq_dir, results_dir = args.fastq_dir, args.results_dir
    fastq_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

    filters = []
    if args.library_name:
        filters.append(f"library_name~{args.library_name}")
    if args.library_strategy:
        filters.append(f"strategy={args.library_strategy}")
    print(f"Querying ENA: {accession} | PAIRED" + "".join(f" | {f}" for f in filters))

    rows, selected = metadata(accession, args.library_name, args.library_strategy)
    print(f"Runs returned: {len(rows)}")
    print(f"Paired-end runs selected: {len(selected)}")
    if not selected:
        raise RuntimeError("No paired-end runs matched. Try removing --library-name / --library-strategy.")

    if args.limit and args.limit > 0:
        selected = selected[:args.limit]
        print(f"Limited run selection to first {len(selected)} samples.")

    manifest = results_dir / f"{accession}_download_manifest.tsv"
    with open(manifest, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(selected)
    print(f"Metadata: {manifest}")

    if args.dry_run:
        for row in selected:
            print(f"  {row['run_accession']}  {row.get('library_name', '')}")
        print("Dry run: nothing downloaded.")
        return 0

    complete = 0
    for i, row in enumerate(selected, 1):
        run = row["run_accession"]
        urls = pair_urls(row)
        dests = [fastq_dir / f"{run}_1.fastq.gz", fastq_dir / f"{run}_2.fastq.gz"]
        print(f"\n[{i}/{len(selected)}] {run}")
        ok = True
        for u, d in zip(urls, dests):
            if valid_fastq(d):
                print(f"  OK: {d.name}")
                continue
            if d.exists():
                d.unlink()
            try:
                download(u, d)
            except Exception as e:
                print(f"  ERROR: {e}")
                ok = False
                break
            if not valid_fastq(d):
                print(f"  ERROR: integrity check failed: {d.name}")
                d.unlink(missing_ok=True)
                ok = False
                break
        if ok and all(valid_fastq(d) for d in dests):
            complete += 1
            print("  Pair ready.")
        else:
            for d in dests:
                if d.exists() and not valid_fastq(d):
                    d.unlink()
            print("  Pair skipped/incomplete.")

    print(f"\nComplete paired samples: {complete}/{len(selected)}")
    print(f"FASTQ directory: {fastq_dir}")
    return 0 if complete else 1


if __name__ == "__main__":
    sys.exit(main())