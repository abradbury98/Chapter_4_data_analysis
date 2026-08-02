#!/bin/bash
#SBATCH --job-name=abricate_kraken2
#SBATCH --account=massey04083
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=200GB
#SBATCH --output=logs/abricate_kraken2_%j.log
#SBATCH --error=logs/abricate_kraken2_%j.err

set -euo pipefail

# ── User-configurable paths ───────────────────────────────────────────────────
READS_DIR=/home/abradbur/00_nesi_projects/00_nesi_projects/massey04083_nobackup/method_development_chapter1/Nanopore_spiked_run_for_thesis/host_removed_reads_2
OUT_DIR=/home/abradbur/00_nesi_projects/00_nesi_projects/massey04083_nobackup/method_development_chapter1/Nanopore_spiked_run_for_thesis/abricate_kraken2_results

# Kraken2 NT database — update this path to the actual NeSI location
DB=/opt/nesi/db/Kraken2/nt

THREADS=16

# Add or remove barcodes as needed
BARCODES=(
    barcode01
    barcode02
    barcode03
    barcode04
    barcode05
    barcode06
)

# ABRicate databases to run
ABRICATE_DBS=(card ncbi resfinder)

# ── Modules ───────────────────────────────────────────────────────────────────
module load ABRicate     # update version if needed, e.g. ABRicate/1.0.1-gimkl-2020a
module load Kraken2      # update version if needed, e.g. Kraken2/2.1.2-GCC-11.3.0

# ── Directory setup ───────────────────────────────────────────────────────────
mkdir -p logs
mkdir -p "${OUT_DIR}/abricate"
mkdir -p "${OUT_DIR}/kraken2"

# ── Main loop ─────────────────────────────────────────────────────────────────
for BARCODE in "${BARCODES[@]}"; do
    READS="${READS_DIR}/all_${BARCODE}_porechop_host_removed.fastq.gz"

    if [ ! -f "${READS}" ]; then
        echo "WARNING: reads file not found for ${BARCODE}: ${READS} — skipping"
        continue
    fi

    echo ""
    echo "====== ${BARCODE} ======"

    # ── ABRicate ─────────────────────────────────────────────────────────────
    for DB in "${ABRICATE_DBS[@]}"; do
        ABRICATE_OUT="${OUT_DIR}/abricate/${BARCODE}_abricate_${DB}.tsv"

        if [ -f "${ABRICATE_OUT}" ]; then
            echo "  ABRicate (${DB}): output exists, skipping"
            continue
        fi

        echo "  ABRicate (${DB}): running..."
        abricate \
            --db "${DB}" \
            --threads "${THREADS}" \
            "${READS}" \
            > "${ABRICATE_OUT}" \
            2> "${OUT_DIR}/abricate/${BARCODE}_abricate_${DB}.log"
        echo "  ABRicate (${DB}): done"
    done

    # ── Kraken2 ──────────────────────────────────────────────────────────────
    KRAKEN2_OUT="${OUT_DIR}/kraken2/${BARCODE}_kraken2"
    KRAKEN2_REPORT="${KRAKEN2_OUT}.report.txt"
    KRAKEN2_OUTPUT="${KRAKEN2_OUT}.output.txt"

    if [ -f "${KRAKEN2_REPORT}" ]; then
        echo "  Kraken2: output exists, skipping"
    else
        echo "  Kraken2: running..."
        kraken2 \
            --db "${DB}" \
            --threads "${THREADS}" \
            --report "${KRAKEN2_REPORT}" \
            --output "${KRAKEN2_OUTPUT}" \
            --gzip-compressed \
            "${READS}" \
            2> "${KRAKEN2_OUT}.log"
        echo "  Kraken2: done"
    fi

done

# ── Combine ABRicate results into per-database summary files ──────────────────
echo ""
echo "====== Combining ABRicate results ======"
for DB in "${ABRICATE_DBS[@]}"; do
    SUMMARY="${OUT_DIR}/abricate/summary_${DB}.tsv"
    abricate --summary "${OUT_DIR}/abricate/"*"_abricate_${DB}.tsv" > "${SUMMARY}"
    echo "  ${DB} summary: ${SUMMARY}"
done

echo ""
echo "====== All done ======"
echo "ABRicate results: ${OUT_DIR}/abricate/"
echo "Kraken2 results:  ${OUT_DIR}/kraken2/"
