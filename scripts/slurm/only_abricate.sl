#!/bin/bash -e
#SBATCH --account=massey04083
#SBATCH --job-name=abricate_only
#SBATCH --time=12:00:00
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=12
#SBATCH --output=log/%x_%j.out
#SBATCH --error=log/%x_%j.err

set -euo pipefail
shopt -s nullglob

# ============================================================================
# ABRicate Only - Reads and Contigs
# Databases: ecoh, ecoli_vf, megares, vfdb, plasmidfinder
# ============================================================================

log_step() { printf '\033[1;32m[%s] %s\033[0m\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$1"; }

# -------------------------
# Configuration
# -------------------------
WD=$(pwd)

INPUT_DIR="/home/abradbur/00_nesi_projects/00_nesi_projects/massey04083_nobackup/poultry_processing_chapter2/EVI_RVS/host_removed_reads"
FLYE_DIR="${WD}/flye_assembly"
ABRICATE_READS_DIR="${WD}/abricate_reads"
THREADS=${SLURM_CPUS_PER_TASK:-8}

ABRICATE_DBS=("ecoh" "ecoli_vf" "megares" "vfdb" "plasmidfinder")

LOG_DIR="${WD}/logs"
mkdir -p "$LOG_DIR" log

# ============================================================================
# STEP 1: ABRicate on READS
# ============================================================================
log_step "Step 1: ABRicate on reads"

module purge
module load seqtk/1.5-GCC-12.3.0
module load ABRicate/1.0.0-GCC-11.3.0-Perl-5.34.1

abricate --version | tee -a "$LOG_DIR/module_versions.txt"

abricate --setupdb > "$LOG_DIR/abricate_setup.log" 2>&1 || log_step "abricate --setupdb may have failed or been run earlier"

log_step "Available ABRicate databases:"
abricate --list 2>&1 | tee -a "$LOG_DIR/abricate_databases.txt"

mkdir -p "$ABRICATE_READS_DIR"
FASTA_TMP="${ABRICATE_READS_DIR}/tmp_fasta"
mkdir -p "$FASTA_TMP"

for fq in "$INPUT_DIR"/*.fastq.gz "$INPUT_DIR"/*.fq.gz; do
    [[ -f "$fq" ]] || continue
    sample=$(basename "$fq" | sed -E 's/_filtered\.fastq\.gz$//; s/\.fastq\.gz$//; s/\.fq\.gz$//')

    log_step "Converting FASTQ to FASTA for ABRicate: $sample"
    fasta_file="${FASTA_TMP}/${sample}_reads.fasta"
    seqtk seq -A "$fq" > "$fasta_file"

    sample_abr_dir="${ABRICATE_READS_DIR}/${sample}"
    mkdir -p "$sample_abr_dir"

    for db in "${ABRICATE_DBS[@]}"; do
        log_step "Running ABRicate ($db) on reads: $sample"
        out_tsv="${sample_abr_dir}/${sample}_reads_${db}.tsv"

        abricate --db "$db" "$fasta_file" > "$out_tsv" 2> "$LOG_DIR/${sample}.abricate_reads_${db}.log" || {
            log_step "ABRicate ($db) failed for $sample reads (see log)"
            continue
        }

        hits=$(tail -n +2 "$out_tsv" | wc -l || echo 0)
        log_step "ABRicate ($db) found $hits hits in $sample reads"
    done

    rm -f "$fasta_file"
done

rmdir "$FASTA_TMP" 2>/dev/null || true

log_step "Step 1 complete: ABRicate on reads finished."

# ============================================================================
# STEP 2: ABRicate on CONTIGS
# ============================================================================
log_step "Step 2: ABRicate on contigs"

module purge
module load ABRicate/1.0.0-GCC-11.3.0-Perl-5.34.1

mkdir -p "${FLYE_DIR}/abricate_results_contigs"

for file in "$FLYE_DIR"/flye_output_*/assembly.fasta; do
    [[ -f "$file" ]] || continue
    sample_dir=$(dirname "$file")
    sample=$(basename "$sample_dir" | sed 's/^flye_output_//')

    sample_abr_dir="${FLYE_DIR}/abricate_results_contigs/${sample}"
    mkdir -p "$sample_abr_dir"

    for db in "${ABRICATE_DBS[@]}"; do
        log_step "Running ABRicate ($db) on contigs: $sample"
        out_tsv="${sample_abr_dir}/${sample}_contigs_${db}.tsv"

        abricate --db "$db" "$file" > "$out_tsv" 2> "$LOG_DIR/${sample}.abricate_contigs_${db}.log" || {
            log_step "ABRicate ($db) failed for $sample contigs (see log)"
            continue
        }

        hits=$(tail -n +2 "$out_tsv" | wc -l || echo 0)
        log_step "ABRicate ($db) found $hits hits in $sample contigs"
    done
done

log_step "Step 2 complete: ABRicate on contigs finished."

# ============================================================================
# FINAL SUMMARY
# ============================================================================
log_step "Pipeline complete."
log_step "  ABRicate (reads):   $ABRICATE_READS_DIR"
log_step "  ABRicate (contigs): ${FLYE_DIR}/abricate_results_contigs"
log_step "  Logs:               $LOG_DIR"
