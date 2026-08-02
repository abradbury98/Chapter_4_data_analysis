#!/bin/bash
#SBATCH --account=massey04083
#SBATCH --job-name=CheckM2_GTDBTk
#SBATCH --time=12:00:00
#SBATCH --mem=180GB
#SBATCH --cpus-per-task=16
#SBATCH --output=log/%x_%j.out
#SBATCH --error=log/%x_%j.err

THREADS=${SLURM_CPUS_PER_TASK:-16}
set -euo pipefail

# --------------------------------------------------
# Set up directories and variables
# --------------------------------------------------
WD=$(pwd)
OUTPUT_DIR="${OUTPUT_DIR:-${WD}/dastool_results}"
LOG_DIR="${LOG_DIR:-${WD}/log}"

mkdir -p "$LOG_DIR"

VER_FILE="$LOG_DIR/checkm2_gtdbtk_$(date '+%Y%m%d_%H%M%S').log"
echo "CheckM2 + GTDB-Tk Pipeline Log - $(date)" > "$VER_FILE"
log_step() {
  printf "[%s] %s\n" "$(date '+%Y-%m-%d %H:%M:%S')" "$1" | tee -a "$VER_FILE"
}

log_step "Working directory: $WD"
log_step "Output directory: $OUTPUT_DIR"

# --------------------------------------------------
# Detect samples from dastool_results directory
# --------------------------------------------------
samples=()
for dir in "$OUTPUT_DIR"/*/; do
    sample=$(basename "$dir")
    if [ -d "$dir" ]; then
        samples+=("$sample")
    fi
done

if [ ${#samples[@]} -eq 0 ]; then
    echo "No sample directories found in $OUTPUT_DIR"
    exit 1
fi

log_step "Found ${#samples[@]} samples: ${samples[*]}"

# ==================================================
# PART 1: Run CheckM2 on DAS Tool bins
# ==================================================
log_step ""
log_step "=========================================="
log_step "PART 1: Running CheckM2 on DAS Tool bins"
log_step "=========================================="

module purge
module load CheckM2/1.0.1-Miniconda3

for sample in "${samples[@]}"; do
    log_step ""

    dastool_bins="$OUTPUT_DIR/${sample}/dastool_output/${sample}_DASTool_bins"

    if [ ! -d "$dastool_bins" ] || [ -z "$(ls -A "$dastool_bins"/*.fa 2>/dev/null)" ]; then
        log_step "⚠️ No DAS Tool bins found for $sample - skipping CheckM2"
        continue
    fi

    bin_count=$(ls "$dastool_bins"/*.fa 2>/dev/null | wc -l)
    log_step "Running CheckM2 on $bin_count bins for $sample"

    checkm2_output="$OUTPUT_DIR/${sample}/checkm2_dastool"
    mkdir -p "$checkm2_output"

    if checkm2 predict \
        --input "$dastool_bins" \
        --output-directory "$checkm2_output" \
        --extension fa \
        --threads $THREADS \
        > "$checkm2_output/checkm2.log" 2>&1; then
        if [ -f "$checkm2_output/quality_report.tsv" ]; then
            log_step "✅ CheckM2 complete for $sample"
            log_step "CheckM2 Summary:"
            head -n 11 "$checkm2_output/quality_report.tsv" | while IFS=$'\t' read -r line; do
                log_step "  $line"
            done
        else
            log_step "⚠️ CheckM2 finished but quality_report.tsv not found for $sample"
        fi
    else
        log_step "⚠️ CheckM2 failed for $sample - check $checkm2_output/checkm2.log"
    fi
done

# ==================================================
# PART 2: Run GTDB-Tk on DAS Tool bins
# ==================================================
log_step ""
log_step "=========================================="
log_step "PART 2: Running GTDB-Tk on DAS Tool bins"
log_step "=========================================="

module purge
module load GTDB-Tk/2.4.0-foss-2023a-Python-3.11.6

for sample in "${samples[@]}"; do
    log_step ""

    dastool_bins="$OUTPUT_DIR/${sample}/dastool_output/${sample}_DASTool_bins"

    if [ ! -d "$dastool_bins" ] || [ -z "$(ls -A "$dastool_bins"/*.fa 2>/dev/null)" ]; then
        log_step "⚠️ No DAS Tool bins found for $sample - skipping GTDB-Tk"
        continue
    fi

    bin_count=$(ls "$dastool_bins"/*.fa 2>/dev/null | wc -l)
    log_step "Running GTDB-Tk on $bin_count bins for $sample"

    gtdbtk_output="$OUTPUT_DIR/${sample}/gtdbtk_dastool"
    mkdir -p "$gtdbtk_output"

    if gtdbtk classify_wf \
        --genome_dir "$dastool_bins" \
        --out_dir "$gtdbtk_output" \
        --extension fa \
        --cpus $THREADS \
        --skip_ani_screen \
        > "$gtdbtk_output/gtdbtk.log" 2>&1; then
        if [ -f "$gtdbtk_output/gtdbtk.bac120.summary.tsv" ] || [ -f "$gtdbtk_output/gtdbtk.ar53.summary.tsv" ]; then
            log_step "✅ GTDB-Tk complete for $sample"
            if [ -f "$gtdbtk_output/gtdbtk.bac120.summary.tsv" ]; then
                bac_count=$(tail -n +2 "$gtdbtk_output/gtdbtk.bac120.summary.tsv" | wc -l)
                log_step "  Bacterial genomes classified: $bac_count"
            fi
            if [ -f "$gtdbtk_output/gtdbtk.ar53.summary.tsv" ]; then
                arc_count=$(tail -n +2 "$gtdbtk_output/gtdbtk.ar53.summary.tsv" | wc -l)
                log_step "  Archaeal genomes classified: $arc_count"
            fi
        else
            log_step "⚠️ GTDB-Tk finished but no summary files found for $sample"
        fi
    else
        log_step "⚠️ GTDB-Tk failed for $sample - check $gtdbtk_output/gtdbtk.log"
    fi
done

# ==================================================
# PART 3: Generate summary report
# ==================================================
log_step ""
log_step "=========================================="
log_step "PART 3: Generating summary report"
log_step "=========================================="

summary_file="$OUTPUT_DIR/SUMMARY_REPORT.txt"

cat > "$summary_file" << 'EOF'
===============================================
DAS Tool Bin Refinement Summary
===============================================
EOF

echo "Generated: $(date)" >> "$summary_file"
echo "" >> "$summary_file"

for sample in "${samples[@]}"; do
    echo "Sample: $sample" >> "$summary_file"
    echo "-------------------" >> "$summary_file"

    # Count DAS Tool bins
    dastool_count=0
    if [ -d "$OUTPUT_DIR/${sample}/dastool_output/${sample}_DASTool_bins" ]; then
        dastool_count=$(ls "$OUTPUT_DIR/${sample}/dastool_output/${sample}_DASTool_bins"/*.fa 2>/dev/null | wc -l)
    fi
    echo "  DAS Tool refined bins: $dastool_count" >> "$summary_file"

    # CheckM2 results
    if [ -f "$OUTPUT_DIR/${sample}/checkm2_dastool/quality_report.tsv" ]; then
        echo "" >> "$summary_file"
        echo "  CheckM2 Quality Assessment:" >> "$summary_file"
        awk -F'\t' 'NR>1 {
            comp=$2; cont=$3;
            sum_comp+=comp; sum_cont+=cont; count++;
            if (comp >= 90 && cont <= 5) hq++;
            if (comp >= 50 && cont <= 10) mq++;
        } END {
            if (count>0) {
                printf "    Total bins assessed: %d\n", count;
                printf "    Average Completeness: %.2f%%\n", sum_comp/count;
                printf "    Average Contamination: %.2f%%\n", sum_cont/count;
                printf "    High Quality (>=90%% comp, <=5%% cont): %d bins\n", hq;
                printf "    Medium Quality (>=50%% comp, <=10%% cont): %d bins\n", mq;
            }
        }' "$OUTPUT_DIR/${sample}/checkm2_dastool/quality_report.tsv" >> "$summary_file"
    else
        echo "  CheckM2: no results (bins absent or run failed)" >> "$summary_file"
    fi

    # GTDB-Tk results
    if [ -f "$OUTPUT_DIR/${sample}/gtdbtk_dastool/gtdbtk.bac120.summary.tsv" ] || \
       [ -f "$OUTPUT_DIR/${sample}/gtdbtk_dastool/gtdbtk.ar53.summary.tsv" ]; then
        echo "" >> "$summary_file"
        echo "  GTDB-Tk Taxonomy:" >> "$summary_file"
        if [ -f "$OUTPUT_DIR/${sample}/gtdbtk_dastool/gtdbtk.bac120.summary.tsv" ]; then
            bac_count=$(tail -n +2 "$OUTPUT_DIR/${sample}/gtdbtk_dastool/gtdbtk.bac120.summary.tsv" | wc -l)
            echo "    Bacterial genomes: $bac_count" >> "$summary_file"
        fi
        if [ -f "$OUTPUT_DIR/${sample}/gtdbtk_dastool/gtdbtk.ar53.summary.tsv" ]; then
            arc_count=$(tail -n +2 "$OUTPUT_DIR/${sample}/gtdbtk_dastool/gtdbtk.ar53.summary.tsv" | wc -l)
            echo "    Archaeal genomes: $arc_count" >> "$summary_file"
        fi
    else
        echo "  GTDB-Tk: no results (bins absent or run failed)" >> "$summary_file"
    fi

    echo "" >> "$summary_file"
done

echo "Output Locations:" >> "$summary_file"
echo "  CheckM2 results: $OUTPUT_DIR/\${sample}/checkm2_dastool/quality_report.tsv" >> "$summary_file"
echo "  GTDB-Tk results: $OUTPUT_DIR/\${sample}/gtdbtk_dastool/" >> "$summary_file"
echo "===============================================" >> "$summary_file"

log_step "Summary report saved to: $summary_file"
cat "$summary_file" | tee -a "$VER_FILE"

log_step ""
log_step "Pipeline complete!"
log_step "Results directory: $OUTPUT_DIR"
