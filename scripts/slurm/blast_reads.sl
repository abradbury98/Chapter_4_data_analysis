#!/bin/bash -e

#SBATCH --job-name=blast_host_check
#SBATCH --account=massey04083
#SBATCH --time=06:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --array=0-5
#SBATCH --output=logs/%A_%a_blast.out
#SBATCH --error=logs/%A_%a_blast.err

# ============================================================================
# BLAST reads against chicken genome to check for residual host sequences
# Array job: one task per FASTQ file (6 files, indices 0-5)
# ============================================================================

module purge
module load BLAST/2.16.0-GCC-12.3.0

# ============================================================================
# USER CONFIGURATION
# ============================================================================

INPUT_DIR="/home/abradbur/00_nesi_projects/00_nesi_projects/massey04083_nobackup/poultry_processing_chapter2/ASC/host_removed_reads"
CHICKEN_REF="/home/abradbur/00_nesi_projects/00_nesi_projects/massey04083_nobackup/poultry_processing_chapter2/ASC_BFBB/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna"
OUTDIR="$(pwd)/blast_host_results"
DB_DIR="$(pwd)/chicken_blastdb"

# ============================================================================

THREADS=${SLURM_CPUS_PER_TASK:-8}
mkdir -p "${OUTDIR}" "${DB_DIR}" logs

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"; }

# ------------------------------------------------------------------
# Step 1: Build chicken BLAST database (only needs to run once —
#         subsequent array tasks will reuse it)
# ------------------------------------------------------------------
if [[ ! -f "${DB_DIR}/chicken.nhr" ]]; then
    log "Building chicken BLAST database"
    makeblastdb \
        -in "${CHICKEN_REF}" \
        -dbtype nucl \
        -out "${DB_DIR}/chicken" \
        -title "GalGal_broiler"
    log "Database built"
else
    log "Chicken BLAST database already exists, skipping build"
fi

# ------------------------------------------------------------------
# Step 2: Resolve FASTQ for this array task
# ------------------------------------------------------------------
FASTQ_FILES=(${INPUT_DIR}/*.fastq.gz)
FASTQ="${FASTQ_FILES[$SLURM_ARRAY_TASK_ID]}"
BASENAME=$(basename "${FASTQ}" .fastq.gz)

log "Array task ${SLURM_ARRAY_TASK_ID}: processing ${BASENAME}"

# ------------------------------------------------------------------
# Step 3: Convert FASTQ.GZ -> FASTA
# ------------------------------------------------------------------
FASTA="${OUTDIR}/${BASENAME}.fasta"
log "Converting FASTQ to FASTA"
zcat "${FASTQ}" | awk 'NR%4==1{ print ">" substr($0,2) } NR%4==2{ print }' > "${FASTA}"
READ_COUNT=$(grep -c "^>" "${FASTA}")
log "Total reads: ${READ_COUNT}"

# ------------------------------------------------------------------
# Step 4: BLAST against chicken genome — top 3 hits per read
# ------------------------------------------------------------------
BLAST_OUT="${OUTDIR}/${BASENAME}_blast_raw.tsv"
log "Running blastn against chicken genome"

blastn \
    -query "${FASTA}" \
    -db "${DB_DIR}/chicken" \
    -out "${BLAST_OUT}" \
    -outfmt "6 qseqid sseqid pident length evalue bitscore stitle" \
    -max_target_seqs 3 \
    -max_hsps 1 \
    -num_threads "${THREADS}" \
    -evalue 1e-5 \
    -perc_identity 80

log "BLAST complete — $(wc -l < "${BLAST_OUT}") chicken hits found"

# ------------------------------------------------------------------
# Step 5: Add rank and sample columns
# ------------------------------------------------------------------
RANKED_OUT="${OUTDIR}/${BASENAME}_chicken_top3.tsv"

awk -v sample="${BASENAME}" '
BEGIN { OFS="\t"; prev=""; rank=0 }
{
    if ($1 != prev) { prev=$1; rank=1 } else { rank++ }
    if (rank <= 3) { print sample, $0, rank }
}
' "${BLAST_OUT}" > "${RANKED_OUT}"

# Summary: how many reads had any chicken hit
READS_WITH_HITS=$(awk '{print $2}' "${RANKED_OUT}" | sort -u | wc -l)
log "Reads with chicken hits: ${READS_WITH_HITS} / ${READ_COUNT}"
log "Done: ${RANKED_OUT}"
log "=============================================="
