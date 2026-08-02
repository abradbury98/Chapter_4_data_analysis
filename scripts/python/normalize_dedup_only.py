#!/usr/bin/env python3
"""
normalize_dedup_only.py
Regenerates ONLY the deduplicated BPKG output (all_samples_BPKG_dedup_long.tsv).
Use when the combined file is already up to date but the dedup file needs updating.
"""

import sys
import pandas as pd
from pathlib import Path

BASE             = Path("/Volumes/Alice/2_POULTRY_PROCESSING/Desktop")
HAMRONIZED_DIR   = BASE / "hamronized"
OUTPUT_DIR       = BASE / "Reads"
BASES_FILE       = BASE / "all_samples_total_bases.tsv"
MIN_COVERAGE_PCT = 80.0
MIN_IDENTITY_PCT = 90.0


def load_bases(path):
    df = pd.read_csv(path, sep="\t")
    col = "bpm_sample" if "bpm_sample" in df.columns else df.columns[0]
    return dict(zip(df[col], df["total_bases"]))


def first_valid(s):
    vals = s.dropna()
    return vals.iloc[0] if len(vals) > 0 else pd.NA


def compute_bpkg(df, sample, total_gbp):
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

    read_span       = (df["input_gene_stop"] - df["input_gene_start"]).abs() + 1
    df["gene_len"]  = df["reference_gene_length"].fillna(read_span)
    df["bases_on_gene"] = df["gene_len"] * df["coverage_percentage"] / 100

    # Resistance: antimicrobial_agent preferred; fall back to drug_class
    if "antimicrobial_agent" in df.columns and "drug_class" in df.columns:
        df["_resistance"] = df["antimicrobial_agent"].fillna(df["drug_class"])
    elif "drug_class" in df.columns:
        df["_resistance"] = df["drug_class"]
    else:
        df["_resistance"] = df.get("antimicrobial_agent", pd.NA)

    # Fill gene_symbol from all_gene_names if missing (dedup files)
    if "all_gene_names" in df.columns:
        df["gene_symbol"] = df["gene_symbol"].fillna(df["all_gene_names"])

    stats = df.groupby(["gene_symbol"], dropna=False).agg(
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
    stats["database"] = "deduplicated"
    return stats


OUTPUT_DIR.mkdir(exist_ok=True)
bases_lookup = load_bases(BASES_FILE)

dedup_tsvs = sorted(p for p in HAMRONIZED_DIR.glob("*/*/*_deduplicated.tsv")
                    if not p.name.startswith("._"))

print(f"Found {len(dedup_tsvs)} deduplicated TSVs")

frames = []
for tsv in dedup_tsvs:
    sample = tsv.stem.replace("_deduplicated", "")
    total_bases = bases_lookup.get(sample)
    if total_bases is None:
        print(f"  SKIP (no total_bases for '{sample}')")
        continue
    total_gbp = total_bases / 1e9
    print(f"  {sample}")
    df = pd.read_csv(tsv, sep="\t", low_memory=False)
    result = compute_bpkg(df, sample, total_gbp)
    if not result.empty:
        frames.append(result)

if frames:
    dedup = pd.concat(frames, ignore_index=True).rename(
        columns={"gene_symbol": "GENE"}
    )
    out = OUTPUT_DIR / "all_samples_BPKG_dedup_long.tsv"
    dedup.to_csv(out, sep="\t", index=False)
    print(f"\nDone: {len(dedup)} rows -> {out.name}")
    res_count = (dedup["resistance"].notna() & (dedup["resistance"] != "")).sum()
    print(f"Rows with resistance class: {res_count}")
else:
    print("No data written.")
