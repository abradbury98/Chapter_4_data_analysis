#!/bin/bash -e
#SBATCH --account=massey04083
#SBATCH --job-name=KRAKEN2_nt
#SBATCH --time=08:00:00
#SBATCH --mem=240GB
#SBATCH --cpus-per-task=24
#SBATCH --output=log/%x_%j.out
#SBATCH --error=log/%x_%j.err

set -euo pipefail
shopt -s nullglob

# -------------------------
# Helper functions
# -------------------------
log_step() { printf '\033[1;32m[%s] %s\033[0m\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$1"; }
log_module_version() {
    local name="$1"; shift
    if "$@" &>/dev/null; then
        echo "[$name] $("$@" 2>&1 | head -n 1)"
    else
        echo "[$name] ERROR running version command"
    fi
}

# -------------------------
# Configuration
# -------------------------
WD=$(pwd)

# Input directory - host-removed reads
INPUT_DIR="/home/abradbur/00_nesi_projects/00_nesi_projects/massey04083_nobackup/poultry_processing_chapter2/ASC/host_removed_reads"

# Output directories
KRAKEN_READS_DIR="${WD}/kraken_reads_nt"
MERGED_RESULTS_DIR="${WD}/merged_amr_taxonomy"
THREADS=${SLURM_CPUS_PER_TASK:-16}
DB=/opt/nesi/db/Kraken2/nt

# Central log directory
LOG_DIR="${WD}/logs"

# ============================================================================
# STEP 1: Kraken2 classification on READS
# ============================================================================
log_step "🧫 Step 1: Kraken2 classification on reads"

module purge
module load Kraken2/2.1.6-GCC-12.3.0

log_module_version "Kraken2" kraken2 --version | tee -a "$LOG_DIR/module_versions.txt"

mkdir -p "$KRAKEN_READS_DIR"

for fq in "$INPUT_DIR"/*.fastq.gz "$INPUT_DIR"/*.fq.gz; do
    [[ -f "$fq" ]] || continue
    sample=$(basename "$fq" | sed -E 's/_filtered\.fastq\.gz$//; s/\.fastq\.gz$//; s/\.fq\.gz$//')
    out_prefix="${KRAKEN_READS_DIR}/${sample}"
    
    log_step "Running Kraken2 on reads for sample: $sample"
    kraken2 \
	--db $DB \
        --threads "$THREADS" \
 	--confidence 0.1 \
        --use-names \
        --output "${out_prefix}_reads_kraken2.out" \
        --report "${out_prefix}_reads_kraken2.report" \
        "$fq" 2> "$LOG_DIR/${sample}.kraken2_reads.log" || {
            log_step "⚠️ Kraken2 failed for $sample reads (see log)"
            continue
        }
    log_step "Kraken2 complete for $sample reads"
done

log_step "✅ Step 1 complete: Kraken2 on reads finished."