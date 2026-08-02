#!/bin/bash
#SBATCH --job-name=checkm2_gtdbtk_bins
#SBATCH --account=massey04083
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=100G
#SBATCH --output=logs/checkm2_gtdbtk_%j.log
#SBATCH --error=logs/checkm2_gtdbtk_%j.err

set -euo pipefail

# ── User-configurable paths and databases ─────────────────────────────────────
BASE=/home/abradbur/00_nesi_projects/00_nesi_projects/massey04083_nobackup/poultry_processing_chapter2/ASC_BFBB
OUT_DIR="${BASE}/all_MAGs_checkm2_gtdbtk"
THREADS=32

GTDBTK_DATA_PATH=/opt/nesi/db/GTDBTk/R232

BARCODES=(barcode01 barcode02 barcode03 barcode04 barcode05 barcode06)
TOOLS=(concoct maxbin2 metabat2)

# Maps input folder barcode names to output sample names
declare -A SAMPLE_NAMES=(
    [barcode01]="CH01_ASC_BFBB"
    [barcode02]="CH02_ASC_BFBB"
    [barcode03]="CH03_ASC_BFBB"
    [barcode04]="CH04_ASC_BFBB"
    [barcode05]="CH05_ASC_BFBB"
    [barcode06]="CH06_ASC_BFBB"
)

# ── Modules ───────────────────────────────────────────────────────────────────
module load CheckM2/1.0.1-Miniconda3
module load GTDB-Tk/2.7.1-foss-2023a-Python-3.11.6

export GTDBTK_DATA_PATH

# ── Directory setup ───────────────────────────────────────────────────────────
mkdir -p "${OUT_DIR}/checkm2"
mkdir -p "${OUT_DIR}/gtdbtk"
mkdir -p "${OUT_DIR}/tmp_bins"

# ── Helper: resolve bin directory and extension per tool ──────────────────────
get_bin_dir_and_ext () {
    local barcode=$1 tool=$2
    local base_dir="${BASE}/bins_${tool}_all_${barcode}_porechop_nanofilt"
    case ${tool} in
        concoct)
            echo "${base_dir}/fasta_bins" "fa"
            ;;
        maxbin2)
            echo "${base_dir}" "fasta"
            ;;
        metabat2)
            echo "${base_dir}" "fa"
            ;;
    esac
}

# ── Main loop: CheckM2 then GTDB-Tk per sample × tool ─────────────────────────
for BARCODE in "${BARCODES[@]}"; do
    for TOOL in "${TOOLS[@]}"; do
        SAMPLE="${SAMPLE_NAMES[$BARCODE]}"
        LABEL="${SAMPLE}_${TOOL}"
        echo ""
        echo "====== ${LABEL} ======"

        read -r BIN_DIR EXT < <(get_bin_dir_and_ext "${BARCODE}" "${TOOL}")

        if [ ! -d "${BIN_DIR}" ]; then
            echo "  WARNING: bin directory not found: ${BIN_DIR} — skipping"
            continue
        fi

        N_BINS=$(find "${BIN_DIR}" -maxdepth 1 -name "*.${EXT}" | wc -l)
        if [ "${N_BINS}" -eq 0 ]; then
            echo "  WARNING: no *.${EXT} files in ${BIN_DIR} — skipping"
            continue
        fi
        echo "  Found ${N_BINS} bins (*.${EXT})"

        # ── CheckM2 ──────────────────────────────────────────────────────────
        CHECKM2_OUT="${OUT_DIR}/checkm2/${LABEL}"
        mkdir -p "${CHECKM2_OUT}"

        if [ -f "${CHECKM2_OUT}/quality_report.tsv" ]; then
            echo "  CheckM2: output exists, skipping"
        else
            echo "  CheckM2: running..."
            checkm2 predict \
                --input "${BIN_DIR}" \
                --output-directory "${CHECKM2_OUT}" \
                --extension "${EXT}" \
                --threads "${THREADS}" \
                --force \
                > "${CHECKM2_OUT}/checkm2.log" 2>&1 || true
            echo "  CheckM2: done"
        fi

        # ── GTDB-Tk ──────────────────────────────────────────────────────────
        GTDBTK_OUT="${OUT_DIR}/gtdbtk/${LABEL}"
        mkdir -p "${GTDBTK_OUT}"

        # GTDB-Tk requires a single extension per run; symlink bins into a
        # temp dir so we can normalise without copying large files.
        TMP_BINS="${OUT_DIR}/tmp_bins/${LABEL}"
        mkdir -p "${TMP_BINS}"
        find "${BIN_DIR}" -maxdepth 1 -name "*.${EXT}" | while read -r BIN; do
            ln -sf "${BIN}" "${TMP_BINS}/$(basename "${BIN}")"
        done

        # In GTDB-Tk 2.x, summary files are written to classify/ subdirectory
        GTDBTK_CLASSIFY_DIR="${GTDBTK_OUT}/classify"
        GTDBTK_DONE=0
        [ -f "${GTDBTK_CLASSIFY_DIR}/${LABEL}.bac120.summary.tsv" ] && GTDBTK_DONE=1
        [ -f "${GTDBTK_CLASSIFY_DIR}/${LABEL}.ar53.summary.tsv"   ] && GTDBTK_DONE=1

        if [ "${GTDBTK_DONE}" -eq 1 ]; then
            echo "  GTDB-Tk: output exists, skipping"
        else
            echo "  GTDB-Tk: running classify_wf..."
            gtdbtk classify_wf \
                --genome_dir "${TMP_BINS}" \
                --out_dir "${GTDBTK_OUT}" \
                --cpus "${THREADS}" \
                --pplacer_cpus 1 \
                --extension "${EXT}" \
                --prefix "${LABEL}" \
                > "${GTDBTK_OUT}/gtdbtk.log" 2>&1 || true
            echo "  GTDB-Tk: done"
        fi

        # Clean up temp symlinks for this sample (keep disk tidy)
        rm -rf "${TMP_BINS}"
    done
done

# ── Summary TSV ───────────────────────────────────────────────────────────────
SUMMARY="${OUT_DIR}/summary_checkm2_gtdbtk.tsv"
echo "======"
echo "Building summary: ${SUMMARY}"

printf 'sample\ttool\tbin\tcompleteness\tcontamination\tgtdbtk_classification\tgtdbtk_source\n' > "${SUMMARY}"

for BARCODE in "${BARCODES[@]}"; do
    for TOOL in "${TOOLS[@]}"; do
        SAMPLE="${SAMPLE_NAMES[$BARCODE]}"
        LABEL="${SAMPLE}_${TOOL}"
        CHECKM2_REPORT="${OUT_DIR}/checkm2/${LABEL}/quality_report.tsv"

        [ -f "${CHECKM2_REPORT}" ] || continue

        # Parse CheckM2 report (cols: Name, Completeness, Contamination, ...)
        tail -n +2 "${CHECKM2_REPORT}" | while IFS=$'\t' read -r BIN_NAME COMPLETENESS CONTAMINATION _REST; do

            GTDBTK_CLASS="N/A"
            GTDBTK_SOURCE="N/A"

            # Search bac120 then ar53 summary for this bin (GTDB-Tk 2.x writes to classify/)
            for DOMAIN in bac120 ar53; do
                GTDBTK_TSV="${OUT_DIR}/gtdbtk/${LABEL}/classify/${LABEL}.${DOMAIN}.summary.tsv"
                [ -f "${GTDBTK_TSV}" ] || continue

                MATCH=$(awk -F'\t' -v b="${BIN_NAME}" 'NR>1 && $1==b {print $2; exit}' "${GTDBTK_TSV}")
                if [ -n "${MATCH}" ]; then
                    GTDBTK_CLASS="${MATCH}"
                    GTDBTK_SOURCE="${DOMAIN}"
                    break
                fi
            done

            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "${SAMPLE}" "${TOOL}" "${BIN_NAME}" \
                "${COMPLETENESS}" "${CONTAMINATION}" \
                "${GTDBTK_CLASS}" "${GTDBTK_SOURCE}"
        done
    done
done >> "${SUMMARY}"

echo "======"
echo "All done. Summary: ${SUMMARY}"
echo "Rows: $(wc -l < "${SUMMARY}")"
