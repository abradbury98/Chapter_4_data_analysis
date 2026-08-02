#!/bin/bash -e

#SBATCH --job-name=filter_contaminants
#SBATCH --account=massey04083   
#SBATCH --time=24:00:00                # Increased for 6 samples
#SBATCH --mem=32G                    
#SBATCH --cpus-per-task=16            
#SBATCH --output=logs/%j_filter.out
#SBATCH --error=logs/%j_filter.err

# ============================================================================
# Filter Contaminating Reads from Nanopore Sequencing Data
# ============================================================================
# This script removes:
#   1. DCS (DNA Control Sequence) reads from Nanopore sequencing
#   2. Human reads (host contamination)
#   3. Chicken broiler reads (host DNA)
#
# Processes all .fastq.gz files in the input directory
# ============================================================================

# Load required modules on NeSI Mahuika
module purge
module load minimap2/2.28-GCC-12.3.0
module load SAMtools/1.22-GCC-12.3.0
module load seqtk/1.5-GCC-12.3.0

# ============================================================================
# USER CONFIGURATION
# ============================================================================

# Input directory containing FASTQ files (gzipped)
INPUT_DIR="/home/abradbur/00_nesi_projects/00_nesi_projects/massey04083_nobackup/poultry_processing_chapter2/EVI_RVS/nanofilt_2"

WD=$(pwd)
# Output directory
OUTDIR="$WD/host_removed_reads"

# Reference genomes
DCS_REF="/home/abradbur/00_nesi_projects/00_nesi_projects/massey04083_nobackup/poultry_processing_chapter2/ASC_BFBB/DCS.fasta"
HUMAN_REF="/home/abradbur/00_nesi_projects/00_nesi_projects/massey04083_nobackup/poultry_processing_chapter2/ASC_BFBB/GCF_000001405.40_GRCh38.p14_genomic.fna"
CHICKEN_REF="/home/abradbur/00_nesi_projects/00_nesi_projects/massey04083_nobackup/poultry_processing_chapter2/ASC_BFBB/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna"

# ============================================================================
# SCRIPT SETUP
# ============================================================================

THREADS=${SLURM_CPUS_PER_TASK:-16}

mkdir -p "${OUTDIR}" logs

# Temporary directory for intermediate files
TMPDIR="${OUTDIR}/tmp_${SLURM_JOB_ID}"
mkdir -p "${TMPDIR}"

# ============================================================================
# FUNCTIONS
# ============================================================================

log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

check_file() {
    if [[ ! -f "$1" ]]; then
        echo "ERROR: File not found: $1"
        exit 1
    fi
}

check_dir() {
    if [[ ! -d "$1" ]]; then
        echo "ERROR: Directory not found: $1"
        exit 1
    fi
}

# ============================================================================
# VALIDATION
# ============================================================================

log_message "Starting contamination filtering pipeline"
log_message "Input directory: ${INPUT_DIR}"
log_message "Output directory: ${OUTDIR}"
log_message "Using ${THREADS} threads"

check_dir "${INPUT_DIR}"
check_file "${DCS_REF}"
check_file "${HUMAN_REF}"
check_file "${CHICKEN_REF}"

# Count input files
FASTQ_FILES=( $(ls ${INPUT_DIR}/*.fastq.gz ${INPUT_DIR}/*.fq.gz 2>/dev/null || true) )
NUM_FILES=${#FASTQ_FILES[@]}

if [[ ${NUM_FILES} -eq 0 ]]; then
    echo "ERROR: No .fastq.gz or .fq.gz files found in ${INPUT_DIR}"
    exit 1
fi

log_message "Found ${NUM_FILES} FASTQ files to process"

# ============================================================================
# STEP 1: Index each reference separately (once)
# ============================================================================

log_message "STEP 1: Indexing reference genomes with minimap2"

DCS_INDEX="${TMPDIR}/dcs.mmi"
HUMAN_INDEX="${TMPDIR}/human.mmi"
CHICKEN_INDEX="${TMPDIR}/chicken.mmi"

minimap2 -t "${THREADS}" -d "${DCS_INDEX}" "${DCS_REF}"
minimap2 -t "${THREADS}" -d "${HUMAN_INDEX}" "${HUMAN_REF}"
minimap2 -t "${THREADS}" -d "${CHICKEN_INDEX}" "${CHICKEN_REF}"

# ============================================================================
# STEP 2: Process each FASTQ file
# ============================================================================

# Summary stats file
SUMMARY_FILE="${OUTDIR}/filtering_summary.txt"
{
    echo "=========================================="
    echo "Contamination Filtering Summary"
    echo "=========================================="
    echo "Date: $(date)"
    echo "Input directory: ${INPUT_DIR}"
    echo ""
    printf "%-40s %12s %12s %10s %10s %10s %10s %10s %10s %10s\n" \
        "Sample" "Total_Reads" "Clean_Reads" "DCS_Rmvd" "DCS_Pct" "Human_Rmvd" "Human_Pct" "Chkn_Rmvd" "Chkn_Pct" "Tot_Pct"
    echo "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
} > "${SUMMARY_FILE}"

SAMPLE_COUNT=0

for INPUT_FASTQ in "${FASTQ_FILES[@]}"; do
    SAMPLE_COUNT=$((SAMPLE_COUNT + 1))
    
    # Get base name for output files
    BASENAME=$(basename "${INPUT_FASTQ}" .fastq.gz)
    BASENAME=$(basename "${BASENAME}" .fq.gz)
    
    log_message "=============================================="
    log_message "Processing sample ${SAMPLE_COUNT}/${NUM_FILES}: ${BASENAME}"
    log_message "=============================================="
    
    # Create sample-specific temp directory
    SAMPLE_TMPDIR="${TMPDIR}/${BASENAME}"
    mkdir -p "${SAMPLE_TMPDIR}"
    
    # ------------------------------------------------------------------
    # Helper: map reads to one reference, emit unmapped as FASTQ
    # Returns mapped read count via stdout of samtools flagstat
    # ------------------------------------------------------------------
    map_and_filter() {
        local REF_INDEX="$1"
        local IN_FASTQ="$2"
        local OUT_FASTQ="$3"
        local TMP_BAM="${SAMPLE_TMPDIR}/tmp_mapped.bam"

        minimap2 -t "${THREADS}" -x map-ont -a --secondary=no \
            "${REF_INDEX}" "${IN_FASTQ}" 2>/dev/null | \
            samtools view -@ "${THREADS}" -bS - > "${TMP_BAM}"

        # Count primary mapped reads only (exclude unmapped=4, secondary=256, supplementary=2048)
        local MAPPED
        MAPPED=$(samtools view -@ "${THREADS}" -c -F 2308 "${TMP_BAM}")

        # Write unmapped reads to next-stage FASTQ
        samtools view -@ "${THREADS}" -b -f 4 "${TMP_BAM}" | \
            samtools fastq -@ "${THREADS}" - | \
            gzip > "${OUT_FASTQ}"

        rm -f "${TMP_BAM}"
        echo "${MAPPED}"
    }

    # ------------------------------------------------------------------
    # Sequential filtering: DCS → Human → Chicken
    # ------------------------------------------------------------------
    TOTAL_READS=$(zcat "${INPUT_FASTQ}" | awk 'NR%4==1' | wc -l)

    log_message "[${BASENAME}] Step 1/3 — removing DCS reads"
    AFTER_DCS="${SAMPLE_TMPDIR}/${BASENAME}_after_dcs.fastq.gz"
    DCS_REMOVED=$(map_and_filter "${DCS_INDEX}" "${INPUT_FASTQ}" "${AFTER_DCS}")

    log_message "[${BASENAME}] Step 2/3 — removing human reads"
    AFTER_HUMAN="${SAMPLE_TMPDIR}/${BASENAME}_after_human.fastq.gz"
    HUMAN_REMOVED=$(map_and_filter "${HUMAN_INDEX}" "${AFTER_DCS}" "${AFTER_HUMAN}")

    log_message "[${BASENAME}] Step 3/3 — removing chicken reads"
    CHICKEN_REMOVED=$(map_and_filter "${CHICKEN_INDEX}" "${AFTER_HUMAN}" "${OUTDIR}/${BASENAME}_filtered.fastq.gz")

    CLEAN_READS=$(zcat "${OUTDIR}/${BASENAME}_filtered.fastq.gz" | awk 'NR%4==1' | wc -l)
    TOTAL_REMOVED=$((TOTAL_READS - CLEAN_READS))

    # ------------------------------------------------------------------
    # Compute percentages (relative to original total)
    # ------------------------------------------------------------------
    pct() { echo "scale=2; $1 * 100 / $2" | bc; }

    if [[ ${TOTAL_READS} -gt 0 ]]; then
        DCS_PCT=$(pct "${DCS_REMOVED}" "${TOTAL_READS}")
        HUMAN_PCT=$(pct "${HUMAN_REMOVED}" "${TOTAL_READS}")
        CHICKEN_PCT=$(pct "${CHICKEN_REMOVED}" "${TOTAL_READS}")
        TOTAL_PCT=$(pct "${TOTAL_REMOVED}" "${TOTAL_READS}")
    else
        DCS_PCT="0.00"; HUMAN_PCT="0.00"; CHICKEN_PCT="0.00"; TOTAL_PCT="0.00"
    fi

    # ------------------------------------------------------------------
    # Individual sample stats file
    # ------------------------------------------------------------------
    STATS_FILE="${OUTDIR}/${BASENAME}_filtering_stats.txt"
    {
        echo "=========================================="
        echo "Contamination Filtering Statistics"
        echo "Sample: ${BASENAME}"
        echo "=========================================="
        echo "Date: $(date)"
        echo "Input file: ${INPUT_FASTQ}"
        echo ""
        echo "Read counts:"
        echo "  Total input reads:          ${TOTAL_READS}"
        echo "  DCS reads removed:          ${DCS_REMOVED} (${DCS_PCT}%)"
        echo "  Human reads removed:        ${HUMAN_REMOVED} (${HUMAN_PCT}%)"
        echo "  Chicken reads removed:      ${CHICKEN_REMOVED} (${CHICKEN_PCT}%)"
        echo "  Total removed:              ${TOTAL_REMOVED} (${TOTAL_PCT}%)"
        echo "  Clean reads retained:       ${CLEAN_READS}"
        echo ""
        echo "Output file: ${OUTDIR}/${BASENAME}_filtered.fastq.gz"
        echo "=========================================="
    } > "${STATS_FILE}"

    # Add row to summary table
    printf "%-40s %12d %12d %10d %9s%% %10d %9s%% %10d %9s%% %9s%%\n" \
        "${BASENAME}" "${TOTAL_READS}" "${CLEAN_READS}" \
        "${DCS_REMOVED}" "${DCS_PCT}" \
        "${HUMAN_REMOVED}" "${HUMAN_PCT}" \
        "${CHICKEN_REMOVED}" "${CHICKEN_PCT}" \
        "${TOTAL_PCT}" >> "${SUMMARY_FILE}"

    # Cleanup sample temp files
    rm -rf "${SAMPLE_TMPDIR}"

    log_message "[${BASENAME}] Complete — ${CLEAN_READS} reads retained (DCS: ${DCS_PCT}%, Human: ${HUMAN_PCT}%, Chicken: ${CHICKEN_PCT}% removed)"
done

# ============================================================================
# FINAL SUMMARY
# ============================================================================

echo "" >> "${SUMMARY_FILE}"
echo "=========================================="  >> "${SUMMARY_FILE}"
echo "Processing complete: $(date)" >> "${SUMMARY_FILE}"
echo "==========================================" >> "${SUMMARY_FILE}"

# ============================================================================
# CLEANUP
# ============================================================================

log_message "Cleaning up temporary files"
rm -rf "${TMPDIR}"

# ============================================================================
# COMPLETION
# ============================================================================

log_message "=============================================="
log_message "Pipeline complete!"
log_message "Processed ${NUM_FILES} samples"
log_message "=============================================="
log_message ""
log_message "Output files in: ${OUTDIR}"
log_message "Summary statistics: ${SUMMARY_FILE}"
log_message ""

# Display summary
cat "${SUMMARY_FILE}"
