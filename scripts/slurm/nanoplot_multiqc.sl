#!/bin/bash -e

#SBATCH --job-name=nanoplot_multiqc
#SBATCH --account=massey04083
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/%j_nanoplot_multiqc.out
#SBATCH --error=logs/%j_nanoplot_multiqc.err

# ============================================================================
# NanoPlot + MultiQC QC Pipeline
# ============================================================================
# Runs NanoPlot on every *.fastq.gz in INPUT_DIR, then aggregates with MultiQC.
# Change INPUT_DIR to reuse this script for different sample folders.
# ============================================================================

# ============================================================================
# USER CONFIGURATION — edit these two lines only
# ============================================================================

INPUT_DIR="/path/to/your/fastq/folder"
WD=$(pwd)

# ============================================================================
# DERIVED PATHS
# ============================================================================

NANOPLOT_DIR="$WD/Nanoplot"
MULTIQC_DIR="$WD/Nanoplot/MultiQC"
LOG_DIR="$WD/logs"

THREADS=${SLURM_CPUS_PER_TASK:-16}
MAX_JOBS=6

# ============================================================================
# SETUP
# ============================================================================

mkdir -p "$NANOPLOT_DIR" "$MULTIQC_DIR" "$LOG_DIR"

log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# ============================================================================
# VALIDATION
# ============================================================================

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "ERROR: INPUT_DIR not found: $INPUT_DIR"
    exit 1
fi

FASTQ_FILES=( "$INPUT_DIR"/*.fastq.gz )
if [[ ! -e "${FASTQ_FILES[0]}" ]]; then
    echo "ERROR: No *.fastq.gz files found in $INPUT_DIR"
    exit 1
fi

log_message "Input directory : $INPUT_DIR"
log_message "NanoPlot output : $NANOPLOT_DIR"
log_message "MultiQC output  : $MULTIQC_DIR"
log_message "Files found     : ${#FASTQ_FILES[@]}"
log_message "Threads         : $THREADS"

# ============================================================================
# STEP 1: NanoPlot (parallelised up to MAX_JOBS)
# ============================================================================

log_message "STEP 1: Running NanoPlot"

module purge
module load NanoPlot/1.43.0-foss-2023a-Python-3.11.6
export PYTHONNOUSERSITE=1   # prevent ~/.local packages from overriding module's plotly/narwhals

job_count=0

for fq in "${FASTQ_FILES[@]}"; do
    sample=$(basename "$fq" .fastq.gz)
    out_dir="$NANOPLOT_DIR/$sample"

    if [[ -f "$out_dir/NanoPlot-report.html" ]]; then
        log_message "  Skipping $sample (already done)"
        continue
    fi

    mkdir -p "$out_dir"
    log_message "  Starting NanoPlot: $sample"

    NanoPlot \
        --fastq "$fq" \
        --threads "$THREADS" \
        --plots dot \
        --maxlength 1000000 \
        --outdir "$out_dir" \
        >"$LOG_DIR/${sample}_NanoPlot.log" 2>&1 &

    job_count=$((job_count + 1))
    if [[ "$job_count" -ge "$MAX_JOBS" ]]; then
        wait -n
        job_count=$((job_count - 1))
    fi
done

wait
log_message "STEP 1 complete: NanoPlot finished for all samples"

# ============================================================================
# STEP 2: MultiQC — per-sample then global aggregate
# ============================================================================

log_message "STEP 2: Running MultiQC"

module purge
module load MultiQC/1.15-gimkl-2022a-Python-3.10.5

# Build a simple HTML index page linking all per-sample reports
INDEX="$MULTIQC_DIR/index.html"
cat > "$INDEX" <<'HTML'
<html>
<head>
  <title>QC Dashboard</title>
  <style>
    body { font-family: sans-serif; margin: 20px; }
    .links { display: flex; flex-wrap: wrap; gap: 8px; margin-bottom: 16px; }
    a { padding: 6px 10px; background: #eef; border-radius: 5px; text-decoration: none; }
    iframe { width: 100%; height: 88vh; border: none; }
  </style>
</head>
<body>
<h1>QC Dashboard</h1>
<div class="links">
HTML

for fq in "${FASTQ_FILES[@]}"; do
    sample=$(basename "$fq" .fastq.gz)
    out_sub="$MULTIQC_DIR/$sample"
    mkdir -p "$out_sub"

    log_message "  MultiQC: $sample"

    multiqc "$NANOPLOT_DIR/$sample" \
        --force \
        --filename "multiqc_${sample}.html" \
        -o "$out_sub" \
        >"$LOG_DIR/${sample}_multiqc.log" 2>&1 || \
        log_message "  WARNING: MultiQC failed for $sample — check $LOG_DIR/${sample}_multiqc.log"

    echo "  <a href=\"${sample}/multiqc_${sample}.html\" target=\"content\">${sample}</a>" >> "$INDEX"
done

# Global report across all samples
log_message "  MultiQC: global aggregate report"
multiqc "$NANOPLOT_DIR" \
    --force \
    --filename "multiqc_all_samples.html" \
    -o "$MULTIQC_DIR" \
    >"$LOG_DIR/multiqc_all.log" 2>&1 || \
    log_message "  WARNING: Global MultiQC failed — check $LOG_DIR/multiqc_all.log"

cat >> "$INDEX" <<'HTML'
</div>
<hr>
<iframe name="content" src="multiqc_all_samples.html"></iframe>
</body>
</html>
HTML

# ============================================================================
# DONE
# ============================================================================

log_message "Pipeline complete"
log_message "  NanoPlot reports : $NANOPLOT_DIR/<sample>/NanoPlot-report.html"
log_message "  MultiQC aggregate: $MULTIQC_DIR/multiqc_all_samples.html"
log_message "  Dashboard index  : $MULTIQC_DIR/index.html"
