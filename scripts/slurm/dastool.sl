#!/bin/bash
#SBATCH --account=massey04083
#SBATCH --job-name=DASTool_Fixed
#SBATCH --time=12:00:00
#SBATCH --mem=180GB
#SBATCH --cpus-per-task=16
#SBATCH --output=log/%x_%j.out
#SBATCH --error=log/%x_%j.err

THREADS=${SLURM_CPUS_PER_TASK:-16}
set -euo pipefail
shopt -s nullglob

# --------------------------------------------------
# Set up directories and variables
# --------------------------------------------------
WD=$(pwd)
RUN_DATE=$(date '+%Y%m%d')
ASSEMBLY_DIR="${ASSEMBLY_DIR:-${WD}/flye_assembly}"
OUTPUT_DIR="${OUTPUT_DIR:-${WD}/dastool_${RUN_DATE}}"
LOG_DIR="${LOG_DIR:-${WD}/log}"
CHECKM_DIR="checkm_${RUN_DATE}"
CHECKM2_DIR="checkm2_${RUN_DATE}"
GTDBTK_DIR="gtdbtk_${RUN_DATE}"

mkdir -p "$LOG_DIR" "$OUTPUT_DIR"

VER_FILE="$LOG_DIR/dastool_fixed_$(date '+%Y%m%d_%H%M%S').log"
echo "📜 DAS Tool Fixed Pipeline Log - $(date)" > "$VER_FILE"
log_step() {
  printf "[%s] %s\n" "$(date '+%Y-%m-%d %H:%M:%S')" "$1" | tee -a "$VER_FILE"
}

log_step "Working directory: $WD"
log_step "Assembly directory: $ASSEMBLY_DIR"
log_step "Output directory: $OUTPUT_DIR"

# --------------------------------------------------
# Manual function to create scaffolds2bin files
# --------------------------------------------------
create_scaffolds2bin() {
    local bin_dir="$1"
    local extension="$2"
    local output_file="$3"
    
    log_step "  Creating scaffolds2bin from $bin_dir with extension .$extension"
    
    > "$output_file"  # Create empty file
    
    for bin_file in "$bin_dir"/*."$extension"; do
        if [ ! -f "$bin_file" ]; then
            continue
        fi
        
        bin_name=$(basename "$bin_file" ."$extension")
        
        # Extract contig names and assign to this bin
        grep "^>" "$bin_file" | sed 's/^>//' | sed 's/ .*//' | while read contig; do
            echo -e "${contig}\t${bin_name}"
        done >> "$output_file"
    done
    
    line_count=$(wc -l < "$output_file")
    log_step "  Created scaffolds2bin with $line_count entries"
}

# --------------------------------------------------
# Detect samples from existing bin directories
# --------------------------------------------------
samples=()
for dir in bins_metabat2_*; do
    if [ -d "$dir" ]; then
        sample=$(echo "$dir" | sed 's/^bins_metabat2_//')
        samples+=("$sample")
    fi
done

if [ ${#samples[@]} -eq 0 ]; then
    echo "❌ No binning directories found! Expected bins_metabat2_* directories"
    exit 1
fi

log_step "Found ${#samples[@]} samples to process: ${samples[*]}"

# ==================================================
# PART 1: Create scaffolds2bin files manually
# ==================================================
log_step ""
log_step "=========================================="
log_step "PART 1: Creating scaffolds2bin files manually"
log_step "=========================================="

for sample in "${samples[@]}"; do
    log_step ""
    log_step "Processing sample: $sample"
    
    mkdir -p "$OUTPUT_DIR/${sample}/dastool_input"
    
    # MetaBAT2 scaffolds2bin
    if [ -d "bins_metabat2_${sample}" ] && [ -n "$(ls -A bins_metabat2_${sample}/*.fa 2>/dev/null)" ]; then
        log_step "Creating MetaBAT2 scaffolds2bin file"
        create_scaffolds2bin "bins_metabat2_${sample}" "fa" "$OUTPUT_DIR/${sample}/dastool_input/metabat2_scaffolds2bin.tsv"
    else
        log_step "⚠️ No MetaBAT2 bins found for $sample"
    fi
    
    # MaxBin2 scaffolds2bin
    if [ -d "bins_maxbin2_${sample}" ] && [ -n "$(ls -A bins_maxbin2_${sample}/*.fasta 2>/dev/null)" ]; then
        log_step "Creating MaxBin2 scaffolds2bin file"
        create_scaffolds2bin "bins_maxbin2_${sample}" "fasta" "$OUTPUT_DIR/${sample}/dastool_input/maxbin2_scaffolds2bin.tsv"
    else
        log_step "⚠️ No MaxBin2 bins found for $sample"
    fi
    
    # CONCOCT scaffolds2bin
    if [ -d "bins_concoct_${sample}/fasta_bins" ] && [ -n "$(ls -A bins_concoct_${sample}/fasta_bins/*.fa 2>/dev/null)" ]; then
        log_step "Creating CONCOCT scaffolds2bin file"
        create_scaffolds2bin "bins_concoct_${sample}/fasta_bins" "fa" "$OUTPUT_DIR/${sample}/dastool_input/concoct_scaffolds2bin.tsv"
    else
        log_step "⚠️ No CONCOCT bins found for $sample"
    fi
done

# ==================================================
# PART 2: Run DAS Tool
# ==================================================
log_step ""
log_step "=========================================="
log_step "PART 2: Running DAS Tool"
log_step "=========================================="

module purge
module load DAS_Tool/1.1.5-gimkl-2022a-R-4.2.1

for sample in "${samples[@]}"; do
    log_step ""
    log_step "Running DAS Tool for sample: $sample"
    
    # Find the assembly file - try multiple patterns
    assembly=""
    for pattern in \
        "${ASSEMBLY_DIR}/flye_output_${sample}/${sample}_assembly.m1000.fasta" \
        "${ASSEMBLY_DIR}/flye_output_${sample}/assembly.m1000.fasta" \
        "${ASSEMBLY_DIR}/flye_output_${sample}/assembly.fasta"; do
        if [ -f "$pattern" ]; then
            assembly="$pattern"
            break
        fi
    done
    
    if [ -z "$assembly" ]; then
        log_step "⚠️ Assembly file not found for $sample, skipping DAS Tool"
        continue
    fi
    
    log_step "Using assembly: $assembly"
    
    # Build DAS Tool input list
    bin_files=""
    bin_labels=""
    binner_count=0

    if [ -f "$OUTPUT_DIR/${sample}/dastool_input/metabat2_scaffolds2bin.tsv" ] && [ -s "$OUTPUT_DIR/${sample}/dastool_input/metabat2_scaffolds2bin.tsv" ]; then
        bin_files="$OUTPUT_DIR/${sample}/dastool_input/metabat2_scaffolds2bin.tsv"
        bin_labels="MetaBAT2"
        binner_count=$((binner_count + 1))
        log_step "  Including MetaBAT2 bins"
    fi

    if [ -f "$OUTPUT_DIR/${sample}/dastool_input/maxbin2_scaffolds2bin.tsv" ] && [ -s "$OUTPUT_DIR/${sample}/dastool_input/maxbin2_scaffolds2bin.tsv" ]; then
        if [ -n "$bin_files" ]; then
            bin_files="${bin_files},$OUTPUT_DIR/${sample}/dastool_input/maxbin2_scaffolds2bin.tsv"
            bin_labels="${bin_labels},MaxBin2"
        else
            bin_files="$OUTPUT_DIR/${sample}/dastool_input/maxbin2_scaffolds2bin.tsv"
            bin_labels="MaxBin2"
        fi
        binner_count=$((binner_count + 1))
        log_step "  Including MaxBin2 bins"
    fi

    if [ -f "$OUTPUT_DIR/${sample}/dastool_input/concoct_scaffolds2bin.tsv" ] && [ -s "$OUTPUT_DIR/${sample}/dastool_input/concoct_scaffolds2bin.tsv" ]; then
        if [ -n "$bin_files" ]; then
            bin_files="${bin_files},$OUTPUT_DIR/${sample}/dastool_input/concoct_scaffolds2bin.tsv"
            bin_labels="${bin_labels},CONCOCT"
        else
            bin_files="$OUTPUT_DIR/${sample}/dastool_input/concoct_scaffolds2bin.tsv"
            bin_labels="CONCOCT"
        fi
        binner_count=$((binner_count + 1))
        log_step "  Including CONCOCT bins"
    fi

    if [ -z "$bin_files" ]; then
        log_step "⚠️ No binning results found for $sample - skipping DAS Tool and downstream steps"
        continue
    fi

    if [ "$binner_count" -lt 3 ]; then
        log_step "⚠️ Only $binner_count/3 binners available for $sample - running DAS Tool on: $bin_labels"
    else
        log_step "Running DAS Tool with all 3 binners: $bin_labels"
    fi

    mkdir -p "$OUTPUT_DIR/${sample}/dastool_output"

    if DAS_Tool \
        -i "$bin_files" \
        -l "$bin_labels" \
        -c "$assembly" \
        -o "$OUTPUT_DIR/${sample}/dastool_output/${sample}" \
        --write_bins \
        --threads $THREADS \
        --search_engine diamond \
        > "$OUTPUT_DIR/${sample}/dastool.log" 2>&1; then
        if [ -d "$OUTPUT_DIR/${sample}/dastool_output/${sample}_DASTool_bins" ]; then
            bin_count=$(ls "$OUTPUT_DIR/${sample}/dastool_output/${sample}_DASTool_bins"/*.fa 2>/dev/null | wc -l)
            log_step "✅ DAS Tool complete for $sample - produced $bin_count bins"
        else
            log_step "⚠️ DAS Tool ran but produced no bins for $sample"
        fi
    else
        log_step "⚠️ DAS Tool failed for $sample - check $OUTPUT_DIR/${sample}/dastool.log"
    fi
done

# ==================================================
# PART 3: Run CheckM on DAS Tool bins
# ==================================================
log_step ""
log_step "=========================================="
log_step "PART 3: Running CheckM on DAS Tool bins"
log_step "=========================================="

module purge
module load CheckM/1.2.3-foss-2023a-Python-3.11.6

for sample in "${samples[@]}"; do
    log_step ""
    
    dastool_bins="$OUTPUT_DIR/${sample}/dastool_output/${sample}_DASTool_bins"
    
    if [ ! -d "$dastool_bins" ] || [ -z "$(ls -A $dastool_bins/*.fa 2>/dev/null)" ]; then
        log_step "⚠️ No DAS Tool bins found for $sample, skipping CheckM"
        continue
    fi
    
    bin_count=$(ls $dastool_bins/*.fa 2>/dev/null | wc -l)
    log_step "Running CheckM on $bin_count DAS Tool bins for $sample"
    
    checkm_output="$OUTPUT_DIR/${sample}/${CHECKM_DIR}"
    mkdir -p "$checkm_output"
    
    if checkm lineage_wf \
        -t $THREADS \
        -x fa \
        --tab_table \
        --file "$checkm_output/checkm_results.tsv" \
        "$dastool_bins" \
        "$checkm_output" \
        > "$checkm_output/checkm.log" 2>&1; then
        if [ -f "$checkm_output/checkm_results.tsv" ]; then
            log_step "✅ CheckM complete for $sample"
            log_step "CheckM Summary:"
            head -n 11 "$checkm_output/checkm_results.tsv" | while IFS=$'\t' read -r line; do
                log_step "  $line"
            done
        else
            log_step "⚠️ CheckM finished but checkm_results.tsv not found for $sample"
        fi
    else
        log_step "⚠️ CheckM failed for $sample - check $checkm_output/checkm.log"
    fi
done

# ==================================================
# PART 4: Run CheckM2 on DAS Tool bins
# ==================================================
log_step ""
log_step "=========================================="
log_step "PART 4: Running CheckM2 on DAS Tool bins"
log_step "=========================================="

module purge
module load CheckM2/1.0.1-Miniconda3

for sample in "${samples[@]}"; do
    log_step ""

    dastool_bins="$OUTPUT_DIR/${sample}/dastool_output/${sample}_DASTool_bins"

    if [ ! -d "$dastool_bins" ] || [ -z "$(ls -A $dastool_bins/*.fa 2>/dev/null)" ]; then
        log_step "⚠️ No DAS Tool bins found for $sample, skipping CheckM2"
        continue
    fi

    bin_count=$(ls $dastool_bins/*.fa 2>/dev/null | wc -l)
    log_step "Running CheckM2 on $bin_count DAS Tool bins for $sample"

    checkm2_output="$OUTPUT_DIR/${sample}/${CHECKM2_DIR}"
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
# PART 5: Run GTDB-Tk on DAS Tool bins
# ==================================================
log_step ""
log_step "=========================================="
log_step "PART 5: Running GTDB-Tk on DAS Tool bins"
log_step "=========================================="

module purge
module load GTDB-Tk/2.4.0-foss-2023a-Python-3.11.6

for sample in "${samples[@]}"; do
    log_step ""
    
    dastool_bins="$OUTPUT_DIR/${sample}/dastool_output/${sample}_DASTool_bins"
    
    if [ ! -d "$dastool_bins" ] || [ -z "$(ls -A $dastool_bins/*.fa 2>/dev/null)" ]; then
        log_step "⚠️ No DAS Tool bins found for $sample, skipping GTDB-Tk"
        continue
    fi
    
    bin_count=$(ls $dastool_bins/*.fa 2>/dev/null | wc -l)
    log_step "Running GTDB-Tk on $bin_count DAS Tool bins for $sample"
    
    gtdbtk_output="$OUTPUT_DIR/${sample}/${GTDBTK_DIR}"
    mkdir -p "$gtdbtk_output"
    
    gtdbtk classify_wf \
        --genome_dir "$dastool_bins" \
        --out_dir "$gtdbtk_output" \
        --extension fa \
        --cpus $THREADS \
        --skip_ani_screen \
        > "$gtdbtk_output/gtdbtk.log" 2>&1
    
    if [ $? -eq 0 ] && { [ -f "$gtdbtk_output/gtdbtk.bac120.summary.tsv" ] || [ -f "$gtdbtk_output/gtdbtk.ar53.summary.tsv" ]; }; then
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
        log_step "⚠️ GTDB-Tk failed for $sample - check $gtdbtk_output/gtdbtk.log"
    fi
done

# ==================================================
# PART 6: Generate summary report
# ==================================================
log_step ""
log_step "=========================================="
log_step "PART 6: Generating summary report"
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
    
    # Count individual bins
    metabat2_count=$(ls bins_metabat2_${sample}/*.fa 2>/dev/null | wc -l)
    maxbin2_count=$(ls bins_maxbin2_${sample}/*.fasta 2>/dev/null | wc -l)
    concoct_count=$(ls bins_concoct_${sample}/fasta_bins/*.fa 2>/dev/null | wc -l)
    
    echo "  Input Binning Results:" >> "$summary_file"
    echo "    MetaBAT2: $metabat2_count bins" >> "$summary_file"
    echo "    MaxBin2:  $maxbin2_count bins" >> "$summary_file"
    echo "    CONCOCT:  $concoct_count bins" >> "$summary_file"
    
    # Count DAS Tool bins
    dastool_count=0
    if [ -d "$OUTPUT_DIR/${sample}/dastool_output/${sample}_DASTool_bins" ]; then
        dastool_count=$(ls "$OUTPUT_DIR/${sample}/dastool_output/${sample}_DASTool_bins"/*.fa 2>/dev/null | wc -l)
    fi
    echo "    DAS Tool: $dastool_count refined bins" >> "$summary_file"
    
    # CheckM results
    if [ -f "$OUTPUT_DIR/${sample}/${CHECKM_DIR}/checkm_results.tsv" ]; then
        echo "" >> "$summary_file"
        echo "  CheckM Quality Assessment:" >> "$summary_file"
        
        awk -F'\t' 'NR>1 {
            comp=$12; cont=$13;
            sum_comp+=comp; sum_cont+=cont; count++;
            if (comp >= 90 && cont <= 5) hq++;
            if (comp >= 50 && cont <= 10) mq++;
        } END {
            if (count>0) {
                printf "    Total bins assessed: %d\n", count;
                printf "    Average Completeness: %.2f%%\n", sum_comp/count;
                printf "    Average Contamination: %.2f%%\n", sum_cont/count;
                printf "    High Quality (≥90%% comp, ≤5%% cont): %d bins\n", hq;
                printf "    Medium Quality (≥50%% comp, ≤10%% cont): %d bins\n", mq;
            }
        }' "$OUTPUT_DIR/${sample}/${CHECKM_DIR}/checkm_results.tsv" >> "$summary_file"
    fi
    
    # CheckM2 results
    if [ -f "$OUTPUT_DIR/${sample}/${CHECKM2_DIR}/quality_report.tsv" ]; then
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
        }' "$OUTPUT_DIR/${sample}/${CHECKM2_DIR}/quality_report.tsv" >> "$summary_file"
    fi

    # GTDB-Tk results
    if [ -f "$OUTPUT_DIR/${sample}/${GTDBTK_DIR}/gtdbtk.bac120.summary.tsv" ] || \
       [ -f "$OUTPUT_DIR/${sample}/${GTDBTK_DIR}/gtdbtk.ar53.summary.tsv" ]; then
        echo "" >> "$summary_file"
        echo "  GTDB-Tk Taxonomy:" >> "$summary_file"
        
        if [ -f "$OUTPUT_DIR/${sample}/${GTDBTK_DIR}/gtdbtk.bac120.summary.tsv" ]; then
            bac_count=$(tail -n +2 "$OUTPUT_DIR/${sample}/${GTDBTK_DIR}/gtdbtk.bac120.summary.tsv" | wc -l)
            echo "    Bacterial genomes: $bac_count" >> "$summary_file"
        fi
        
        if [ -f "$OUTPUT_DIR/${sample}/${GTDBTK_DIR}/gtdbtk.ar53.summary.tsv" ]; then
            arc_count=$(tail -n +2 "$OUTPUT_DIR/${sample}/${GTDBTK_DIR}/gtdbtk.ar53.summary.tsv" | wc -l)
            echo "    Archaeal genomes: $arc_count" >> "$summary_file"
        fi
    fi
    
    echo "" >> "$summary_file"
done

echo "" >> "$summary_file"
echo "Output Locations:" >> "$summary_file"
echo "  DAS Tool bins: $OUTPUT_DIR/\${sample}/dastool_output/\${sample}_DASTool_bins/" >> "$summary_file"
echo "  CheckM results: $OUTPUT_DIR/\${sample}/${CHECKM_DIR}/checkm_results.tsv" >> "$summary_file"
echo "  CheckM2 results: $OUTPUT_DIR/\${sample}/${CHECKM2_DIR}/quality_report.tsv" >> "$summary_file"
echo "  GTDB-Tk results: $OUTPUT_DIR/\${sample}/${GTDBTK_DIR}/" >> "$summary_file"
echo "===============================================" >> "$summary_file"

log_step "Summary report saved to: $summary_file"
cat "$summary_file" | tee -a "$VER_FILE"

log_step ""
log_step "🎉 DAS Tool refinement pipeline complete!"
log_step "Results directory: $OUTPUT_DIR"