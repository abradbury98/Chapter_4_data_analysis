#!/usr/bin/env python3
"""
normalize_abricate.py  (updated for harmonized outputs)

Reads hAMRonization outputs from /Volumes/Alice/2_POULTRY_PROCESSING/Desktop/hamronized/
and computes BPKG (Bases Per Kilobase per Gigabase) per gene per sample.

Per-database BPKG : reads *_combined.tsv  (all hits, all databases)
Unified BPKG      : reads *_deduplicated.tsv (cross-db deduplicated hits)

NOTE: all_samples_total_bases.tsv must use the new sample naming
      (e.g. CH01_EVI_BFBB, not CH01_EVI_CAT). Update that file first.

Outputs to /Volumes/Alice/2_POULTRY_PROCESSING/Desktop/Reads/:
  all_samples_BPKG_long.tsv        – per-database BPKG, long format
  all_samples_BPKG_dedup_long.tsv  – unified (deduplicated) BPKG, long format
  BPKG_wide_<database>.tsv         – wide format per database
"""

import sys
import pandas as pd
from pathlib import Path

# --- Configuration ---
BASE             = Path("/Volumes/Alice/2_POULTRY_PROCESSING/Desktop")
HAMRONIZED_DIR   = BASE / "hamronized"
OUTPUT_DIR       = BASE / "Reads"
BASES_FILE       = BASE / "all_samples_total_bases.tsv"
MIN_COVERAGE_PCT = 80.0
MIN_IDENTITY_PCT = 90.0


def load_bases(path):
    df = pd.read_csv(path, sep="\t")
    # bpm_sample uses the new BFBB naming that matches hamronized output; fall back to abricate_sample or first col
    if "bpm_sample" in df.columns:
        col = "bpm_sample"
    elif "abricate_sample" in df.columns:
        col = "abricate_sample"
    else:
        col = df.columns[0]
    return dict(zip(df[col], df["total_bases"]))


def first_valid(s):
    vals = s.dropna()
    return vals.iloc[0] if len(vals) > 0 else pd.NA


def compute_bpkg(df, sample, total_gbp, group_by_db):
    """Filter reads, compute bases_on_gene, aggregate to gene level, compute BPKG."""
    df = df.copy()
    for col in ["coverage_percentage", "sequence_identity",
                "reference_gene_length", "input_gene_start", "input_gene_stop"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df[
        (df["coverage_percentage"] >= MIN_COVERAGE_PCT) &
        (df["sequence_identity"]   >= MIN_IDENTITY_PCT)
    ].copy()
    if df.empty:
        return pd.DataFrame()

    # Gene length: use reference_gene_length; fall back to read-span if missing
    read_span          = (df["input_gene_stop"] - df["input_gene_start"]).abs() + 1
    df["gene_len"]     = df["reference_gene_length"].fillna(read_span)
    df["bases_on_gene"] = df["gene_len"] * df["coverage_percentage"] / 100

    # Resistance: antimicrobial_agent preferred; fall back to drug_class
    if "antimicrobial_agent" in df.columns and "drug_class" in df.columns:
        df["_resistance"] = df["antimicrobial_agent"].fillna(df["drug_class"])
    elif "drug_class" in df.columns:
        df["_resistance"] = df["drug_class"]
    else:
        df["_resistance"] = df.get("antimicrobial_agent", pd.NA)

    group_cols = (["gene_symbol", "reference_database_name"]
                  if group_by_db else ["gene_symbol"])

    stats = df.groupby(group_cols, dropna=False).agg(
        n_reads      = ("input_sequence_id", "nunique"),
        n_hits       = ("input_sequence_id", "count"),
        bases_mapped = ("bases_on_gene",     "sum"),
        gene_length  = ("gene_len",          "first"),
        resistance   = ("_resistance",        first_valid),
        product      = ("gene_name",          first_valid),
    ).reset_index()

    stats = stats[stats["gene_length"] > 0].copy()
    stats["BPKG"]   = (stats["bases_mapped"] /
                       (stats["gene_length"] / 1000) /
                       total_gbp)
    stats["sample"] = sample
    return stats


# ---------------------------------------------------------------------------
OUTPUT_DIR.mkdir(exist_ok=True)
bases_lookup = load_bases(BASES_FILE)

combined_tsvs = sorted(p for p in HAMRONIZED_DIR.glob("*/*/*_combined.tsv")
                       if not p.name.startswith("._"))
dedup_tsvs    = sorted(p for p in HAMRONIZED_DIR.glob("*/*/*_deduplicated.tsv")
                       if not p.name.startswith("._"))

if not combined_tsvs:
    sys.exit(f"No *_combined.tsv files found under {HAMRONIZED_DIR}")

per_db_frames = []
dedup_frames  = []

for tsv in combined_tsvs:
    sample = tsv.stem.replace("_combined", "")
    total_bases = bases_lookup.get(sample)
    if total_bases is None:
        print(f"SKIP (no total_bases for '{sample}'): check {BASES_FILE.name}")
        continue
    total_gbp = total_bases / 1e9
    print(f"Per-DB  : {sample}")
    df = pd.read_csv(tsv, sep="\t", low_memory=False)
    result = compute_bpkg(df, sample, total_gbp, group_by_db=True)
    if not result.empty:
        per_db_frames.append(result)
    else:
        print(f"  EMPTY after filters")

for tsv in dedup_tsvs:
    sample = tsv.stem.replace("_deduplicated", "")
    total_bases = bases_lookup.get(sample)
    if total_bases is None:
        continue
    total_gbp = total_bases / 1e9
    print(f"Unified : {sample}")
    df = pd.read_csv(tsv, sep="\t", low_memory=False)
    if "all_gene_names" in df.columns:
        df["gene_symbol"] = df["gene_symbol"].fillna(df["all_gene_names"])
    result = compute_bpkg(df, sample, total_gbp, group_by_db=False)
    if not result.empty:
        result["reference_database_name"] = "deduplicated"
        dedup_frames.append(result)

# Rename to Rmd-compatible column names
RENAME = {"gene_symbol": "GENE", "reference_database_name": "database"}

if per_db_frames:
    per_db = pd.concat(per_db_frames, ignore_index=True).rename(columns=RENAME)
    out = OUTPUT_DIR / "all_samples_BPKG_long.tsv"
    per_db.to_csv(out, sep="\t", index=False)
    print(f"\nPer-database : {len(per_db)} rows  ->  {out.name}")
    for db, grp in per_db.groupby("database"):
        wide = grp.pivot_table(
            index=["GENE", "gene_length", "resistance", "product"],
            columns="sample", values="BPKG",
            fill_value=0, aggfunc="sum",
            dropna=False
        ).reset_index()
        wide.columns.name = None
        wide.to_csv(OUTPUT_DIR / f"BPKG_wide_{db}.tsv", sep="\t", index=False)
        print(f"  {db}: {len(wide)} genes")

if dedup_frames:
    dedup = pd.concat(dedup_frames, ignore_index=True).rename(columns=RENAME)
    out = OUTPUT_DIR / "all_samples_BPKG_dedup_long.tsv"
    dedup.to_csv(out, sep="\t", index=False)
    print(f"Unified      : {len(dedup)} rows  ->  {out.name}")

print("\nDone.")
