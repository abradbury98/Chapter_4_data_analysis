#!/usr/bin/env python3
"""
hamronize_and_dedup.py

Scans the Desktop for abricate_reads_* folders, processes each as one
enrichment type. Within each folder, barcodes are ranked by number so
the lowest barcode = CH01 and the highest = CH06. Produces one
deduplicated TSV report per enrichment folder.

Output files land in: ~/Desktop/hamronized/<ENRICHMENT>/
"""

import os
import re
import sys
import subprocess
import pandas as pd
from glob import glob
from pathlib import Path

# --- Configuration ---
DESKTOP      = Path("/Volumes/Alice/2_POULTRY_PROCESSING/Desktop")
OUTPUT_DIR   = DESKTOP / "hamronized"
ABRICATE_VERSION = "1.0.0"
DB_VERSION   = "1"
OVERLAP_THRESHOLD = 0.5
DB_PRIORITY  = ["ncbi", "resfinder", "card", "megares",
                "vfdb", "ecoli_vf", "plasmidfinder"]

DB_COL     = "reference_database_name"
BARCODE_RE = re.compile(r"all_barcode(\d+)_porechop_nanofilt$")
TSV_DB_RE  = re.compile(r"_reads_(.+)\.tsv$")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def reciprocal_overlap(s1, e1, s2, e2):
    overlap = max(0, min(e1, e2) - max(s1, s2))
    shorter = min(e1 - s1, e2 - s2)
    return overlap / shorter if shorter > 0 else 0


def db_rank(db):
    try:
        return DB_PRIORITY.index(db)
    except ValueError:
        return len(DB_PRIORITY)


def cluster_read_hits(rows):
    """Greedy single-linkage clustering by coordinate overlap."""
    rows = sorted(rows, key=lambda r: r["_s"])
    clusters = []
    for row in rows:
        placed = False
        for cluster in clusters:
            for member in cluster:
                if reciprocal_overlap(row["_s"], row["_e"],
                                      member["_s"], member["_e"]) >= OVERLAP_THRESHOLD:
                    cluster.append(row)
                    placed = True
                    break
            if placed:
                break
        if not placed:
            clusters.append([row])
    return clusters


def deduplicate(df, sample_name):
    df = df.copy()
    df["input_gene_start"] = pd.to_numeric(df["input_gene_start"], errors="coerce")
    df["input_gene_stop"]  = pd.to_numeric(df["input_gene_stop"],  errors="coerce")
    df = df.dropna(subset=["input_gene_start", "input_gene_stop"])
    df["input_gene_start"] = df["input_gene_start"].astype(int)
    df["input_gene_stop"]  = df["input_gene_stop"].astype(int)
    df["_s"] = df[["input_gene_start", "input_gene_stop"]].min(axis=1)
    df["_e"] = df[["input_gene_start", "input_gene_stop"]].max(axis=1)

    by_read = {}
    for r in df.to_dict("records"):
        by_read.setdefault(r["input_sequence_id"], []).append(r)

    dedup_rows = []
    cluster_id = 0
    for hits in by_read.values():
        for cluster in cluster_read_hits(hits):
            cluster_id += 1
            dbs = sorted(set(r[DB_COL] for r in cluster), key=db_rank)
            gene_names = sorted(set(
                str(r["gene_symbol"]) for r in cluster
                if pd.notna(r["gene_symbol"])
            ))
            rep = sorted(cluster, key=lambda r: db_rank(r[DB_COL]))[0].copy()
            rep["cluster_id"]         = cluster_id
            rep["databases_detected"] = ";".join(dbs)
            rep["n_databases"]        = len(dbs)
            rep["all_gene_names"]     = ";".join(gene_names)
            rep["sample"]             = sample_name
            dedup_rows.append(rep)

    out = pd.DataFrame(dedup_rows).drop(columns=["_s", "_e"])
    front = ["sample", "cluster_id", "n_databases", "databases_detected",
             "all_gene_names", "gene_symbol", "input_sequence_id",
             "input_gene_start", "input_gene_stop", DB_COL,
             "antimicrobial_agent", "coverage_percentage", "sequence_identity"]
    rest = [c for c in out.columns if c not in front]
    return out[front + rest]


# ---------------------------------------------------------------------------
# Discover enrichment folders  (abricate_reads_*)
# ---------------------------------------------------------------------------
OUTPUT_DIR.mkdir(exist_ok=True)

enrich_dirs = sorted(DESKTOP.glob("abricate_reads_*"))
if not enrich_dirs:
    sys.exit(f"No abricate_reads_* folders found on Desktop: {DESKTOP}")

print(f"Found {len(enrich_dirs)} enrichment folder(s):")
for d in enrich_dirs:
    print(f"  {d.name}")
print()

# ---------------------------------------------------------------------------
# Process each enrichment folder
# ---------------------------------------------------------------------------
for enrich_dir in enrich_dirs:
    # e.g. "abricate_reads_EVI_BFBB" -> "EVI_BFBB"
    enrichment = enrich_dir.name.replace("abricate_reads_", "")
    enrich_out = OUTPUT_DIR / enrichment
    enrich_out.mkdir(exist_ok=True)

    print(f"{'='*60}")
    print(f"Enrichment: {enrichment}")

    # Find barcode subdirectories and sort by barcode number
    barcode_dirs = []
    for d in enrich_dir.iterdir():
        if d.is_dir():
            m = BARCODE_RE.match(d.name)
            if m:
                barcode_dirs.append((int(m.group(1)), d))

    if not barcode_dirs:
        print(f"  No barcode directories found, skipping.\n")
        continue

    barcode_dirs.sort(key=lambda x: x[0])
    n = len(barcode_dirs)
    print(f"  {n} barcode(s): "
          + ", ".join(f"barcode{b:02d}->CH{r:02d}" for r, (b, _) in enumerate(barcode_dirs, 1)))
    print()

    all_sample_frames = []

    for rank, (barcode_num, barcode_dir) in enumerate(barcode_dirs, start=1):
        sample_name = f"CH{rank:02d}_{enrichment}"
        sample_out  = enrich_out / sample_name
        sample_out.mkdir(exist_ok=True)

        print(f"  --- {sample_name} (barcode{barcode_num:02d}) ---")

        tsv_files = sorted(barcode_dir.glob("*.tsv"))
        if not tsv_files:
            print(f"    No TSV files found, skipping.")
            continue

        # Step 1: Hamronize each database TSV
        json_files = []
        for tsv in tsv_files:
            m = TSV_DB_RE.search(tsv.name)
            db_name = m.group(1) if m else tsv.stem

            with open(tsv) as f:
                data_lines = [l for l in f if not l.startswith("#") and l.strip()]
            if not data_lines:
                print(f"    SKIP (empty): {db_name}")
                continue

            json_out = sample_out / f"{db_name}.json"
            cmd = [
                "hamronize", "abricate", str(tsv),
                "--analysis_software_version", ABRICATE_VERSION,
                "--reference_database_version", DB_VERSION,
                "--output", str(json_out),
            ]
            print(f"    Hamronizing: {db_name}")
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                print(f"      ERROR: {result.stderr.strip()}")
                continue
            json_files.append(str(json_out))

        if not json_files:
            print(f"    No files hamronized successfully, skipping.")
            continue

        # Step 2: Summarize databases -> combined TSV
        combined_tsv = sample_out / f"{sample_name}_combined.tsv"
        cmd = ["hamronize", "summarize", "--summary_type", "tsv",
               "--output", str(combined_tsv)] + json_files
        print(f"    Summarizing {len(json_files)} database(s)")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"    ERROR summarizing: {result.stderr.strip()}")
            continue

        # Step 3: Deduplicate
        df = pd.read_csv(combined_tsv, sep="\t", low_memory=False)
        req = [DB_COL, "input_sequence_id", "input_gene_start", "input_gene_stop"]
        missing = [c for c in req if c not in df.columns]
        if missing:
            print(f"    ERROR: missing columns {missing}")
            continue

        print(f"    Deduplicating {len(df)} hits")
        dedup_df = deduplicate(df, sample_name)

        dedup_out = sample_out / f"{sample_name}_deduplicated.tsv"
        dedup_df.to_csv(dedup_out, sep="\t", index=False)

        multi = (dedup_df["n_databases"] > 1).sum()
        print(f"    Total: {len(df)}  |  Clusters: {len(dedup_df)}  |  Multi-DB: {multi}")
        print(f"    -> {dedup_out.name}\n")

        all_sample_frames.append(dedup_df)

    # Step 4: Merge all samples for this enrichment into one report
    if all_sample_frames:
        merged = pd.concat(all_sample_frames, ignore_index=True)
        merged_out = enrich_out / f"{enrichment}_all_samples_deduplicated.tsv"
        merged.to_csv(merged_out, sep="\t", index=False)
        print(f"  Enrichment report: {merged_out.name}")
        print(f"  Total clusters across all samples: {len(merged)}\n")

print("Done.")
