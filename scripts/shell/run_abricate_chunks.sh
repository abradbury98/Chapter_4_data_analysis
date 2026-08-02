#!/bin/bash
# Splits a fastq into chunks and runs abricate sequentially
# Usage: bash run_abricate_chunks.sh

FASTQ="/Users/alicebradbury/Downloads/all_barcode01_porechop_nanofilt_filtered.fastq"
OUTDIR="/Users/alicebradbury/Downloads/abricate_chunks"
DB="ncbi"
THREADS=2
CHUNKS=4

mkdir -p "${OUTDIR}/chunks"

echo "Counting reads..."
TOTAL_LINES=$(wc -l < "${FASTQ}")
TOTAL_READS=$(( TOTAL_LINES / 4 ))
echo "Total reads: ${TOTAL_READS}"

# Lines per chunk — must be divisible by 4 (each fastq read = 4 lines)
READS_PER_CHUNK=$(( (TOTAL_READS + CHUNKS - 1) / CHUNKS ))
LINES_PER_CHUNK=$(( READS_PER_CHUNK * 4 ))

echo "Splitting into ${CHUNKS} chunks (~${READS_PER_CHUNK} reads each)..."
split -l "${LINES_PER_CHUNK}" "${FASTQ}" "${OUTDIR}/chunks/chunk_"
echo "Done splitting."
echo ""

# ── Activate conda env ────────────────────────────────────────────────────────
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate abricate_x86

# ── Run abricate on each chunk sequentially ───────────────────────────────────
for CHUNK in "${OUTDIR}"/chunks/chunk_*; do
    CHUNK_NAME=$(basename "${CHUNK}")
    echo "=============================="
    echo "Running: ${CHUNK_NAME}"
    echo "Started: $(date)"
    echo "=============================="

    abricate --db "${DB}" --threads "${THREADS}" --minid 80 --mincov 80 --noheader "${CHUNK}" \
        > "${OUTDIR}/chunks/${CHUNK_NAME}_results.tab"

    echo "Finished: $(date)"
    echo "Cooling down 30 seconds..."
    sleep 30
    echo ""
done

# ── Combine results ───────────────────────────────────────────────────────────
echo "Combining results..."
echo -e "#FILE\tSEQUENCE\tSTART\tEND\tSTRAND\tGENE\tCOVERAGE\tCOVERAGE_MAP\tGAPS\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT\tRESISTANCE" \
    > "${OUTDIR}/barcode01_${DB}_combined.tab"
cat "${OUTDIR}"/chunks/chunk_*_results.tab >> "${OUTDIR}/barcode01_${DB}_combined.tab"

echo ""
echo "Done! Combined results: ${OUTDIR}/barcode01_${DB}_combined.tab"
echo "Total hits: $(grep -c '' ${OUTDIR}/barcode01_${DB}_combined.tab)"

# ── Clean up chunks ───────────────────────────────────────────────────────────
read -p "Delete chunk files? (y/n): " CLEAN
if [[ "$CLEAN" == "y" ]]; then
    rm -rf "${OUTDIR}/chunks"
    echo "Chunks deleted."
fi
