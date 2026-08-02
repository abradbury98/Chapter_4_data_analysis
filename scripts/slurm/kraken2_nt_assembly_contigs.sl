#!/bin/bash -e
#SBATCH --account=massey04083
#SBATCH --job-name=KRAKEN2_nt_assembly
#SBATCH --time=24:00:00
#SBATCH --mem=500GB
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
ASSEMBLY_BASE="/home/abradbur/00_nesi_projects/00_nesi_projects/massey04083_nobackup/poultry_processing_chapter2/EVI_RVS/flye_assembly"
KRAKEN_OUT_DIR="${ASSEMBLY_BASE}/kraken_nt"
LOG_DIR="${KRAKEN_OUT_DIR}/logs"
THREADS=${SLURM_CPUS_PER_TASK:-16}
DB=/opt/nesi/db/Kraken2/nt

mkdir -p "$KRAKEN_OUT_DIR" "$LOG_DIR"

# ============================================================================
# STEP 1: Kraken2 classification on assemblies
# ============================================================================
log_step "Step 1: Kraken2 classification on assemblies (nt database)"

module purge
module load Kraken2/2.1.6-GCC-12.3.0

log_module_version "Kraken2" kraken2 --version | tee -a "$LOG_DIR/module_versions.txt"

for sample_dir in "$ASSEMBLY_BASE"/flye_output_*/; do
    [[ -d "$sample_dir" ]] || continue

    # Extract sample name by stripping flye_output_ prefix
    sample=$(basename "$sample_dir" | sed 's/^flye_output_//')

    # Expected FASTA files
    fasta_m1000="${sample_dir}${sample}_assembly.m1000.fasta"
    fasta_all="${sample_dir}assembly.fasta"

    out_dir="${KRAKEN_OUT_DIR}/${sample}"
    mkdir -p "$out_dir"

    # --- Run on >1000bp filtered assembly ---
    if [[ -f "$fasta_m1000" ]]; then
        log_step "Running Kraken2 on m1000 assembly: $sample"
        kraken2 \
            --db "$DB" \
            --threads "$THREADS" \
            --confidence 0.1 \
            --memory-mapping \
            --use-names \
            --output  "${out_dir}/${sample}_m1000_kraken2.out" \
            --report  "${out_dir}/${sample}_m1000_kraken2.report" \
            "$fasta_m1000" \
            2> "$LOG_DIR/${sample}_m1000.kraken2.log" || {
                log_step "WARNING: Kraken2 failed for ${sample} m1000 assembly (see log)"
            }
        log_step "Kraken2 complete for ${sample} m1000 assembly"
    else
        log_step "WARNING: m1000 fasta not found, skipping: $fasta_m1000"
    fi

    # --- Run on all-contig assembly ---
    if [[ -f "$fasta_all" ]]; then
        log_step "Running Kraken2 on full assembly: $sample"
        kraken2 \
            --db "$DB" \
            --threads "$THREADS" \
            --confidence 0.1 \
            --memory-mapping \
            --use-names \
            --output  "${out_dir}/${sample}_all_kraken2.out" \
            --report  "${out_dir}/${sample}_all_kraken2.report" \
            "$fasta_all" \
            2> "$LOG_DIR/${sample}_all.kraken2.log" || {
                log_step "WARNING: Kraken2 failed for ${sample} full assembly (see log)"
            }
        log_step "Kraken2 complete for ${sample} full assembly"
    else
        log_step "WARNING: assembly.fasta not found, skipping: $fasta_all"
    fi

done

log_step "Step 1 complete: Kraken2 on assemblies finished."
