#!/bin/bash
# Run abricate (ecoh + ecoli_vf) on all MAGs in a directory
# Usage: bash run_abricate_MAGs.sh

INDIR="/Users/alicebradbury/Desktop/E.coli_NZ_CAT"
OUTDIR="/Users/alicebradbury/Desktop/E.coli_NZ_CAT/abricate_results"
DATABASES=("ecoh" "ecoli_vf" "megares")
THREADS=4

mkdir -p "${OUTDIR}"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate abricate_x86

for DB in "${DATABASES[@]}"; do
    echo "=============================="
    echo "Database: ${DB}"
    echo "=============================="

    DB_OUTDIR="${OUTDIR}/${DB}"
    mkdir -p "${DB_OUTDIR}"

    for FILE in "${INDIR}"/*.fa "${INDIR}"/*.fasta; do
        [[ -f "$FILE" ]] || continue
        SAMPLE=$(basename "${FILE}" | sed 's/\.\(fa\|fasta\)$//')
        echo "  Processing: ${SAMPLE}"
        abricate --db "${DB}" --threads "${THREADS}" --minid 80 --mincov 80 \
            "${FILE}" > "${DB_OUTDIR}/${SAMPLE}_${DB}.tab"
    done

    echo "  Generating summary..."
    abricate --summary "${DB_OUTDIR}"/*_${DB}.tab > "${OUTDIR}/summary_${DB}.tab"
    echo "  Summary: ${OUTDIR}/summary_${DB}.tab"
    echo ""
done

echo "All done! Results in: ${OUTDIR}/"
