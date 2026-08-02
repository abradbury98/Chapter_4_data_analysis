#!/bin/bash

#SBATCH --job-name=host_dep_benchmark
#SBATCH --account=massey04083
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/%j_benchmark.out
#SBATCH --error=logs/%j_benchmark.err

set -eo pipefail

# ============================================================================
# Host Depletion Parameter Benchmark — Nanopore Metagenomics
# ============================================================================
# Tests multiple parameter sets for DCS / human / chicken host read removal
# using two aligners:
#   • minimap2  — designed for long reads; the production aligner
#   • Bowtie2   — designed for short reads; included for comparison.
#                 NOTE: Bowtie2 is NOT optimised for ONT reads. It uses
#                 seed-extension on short k-mers and may miss long-read
#                 alignments that minimap2 catches. Expect lower host removal
#                 rates with Bowtie2 and treat those results as a lower bound.
#                 --local mode is used throughout as it performs better on
#                 long reads than end-to-end.
#
# Residual contamination assessed with:
#   1. Kraken2 — taxonomic breakdown (% unclassified, human, chicken)
#   2. BLAST vs chicken genome — direct alignment check for residual chicken
#
# Key tradeoff:
#   More sensitive removal (lower k / lower score) → fewer chicken reads
#   survive, but higher risk of removing genuine pathogen reads that share
#   short regions of similarity with the host genome.
#
# Output: benchmark_results/benchmark_summary.tsv — one row per sample × param set
# ============================================================================

module purge
module load minimap2/2.30-GCC-12.3.0
module load Bowtie2/2.5.4-GCC-12.3.0
module load HTSlib/1.23.1-GCC-12.3.0
module load SAMtools/1.23.1-GCC-12.3.0
module load seqtk/1.5-GCC-12.3.0
module load BLAST/2.16.0-GCC-12.3.0
module load Kraken2/2.1.6-GCC-12.3.0   

# ============================================================================
# USER CONFIGURATION
# ============================================================================

INPUT_DIR="/home/abradbur/00_nesi_projects/00_nesi_projects/massey04083_nobackup/poultry_processing_chapter2/host_removal_assess/fasta_files"
WD=$(pwd)
OUTDIR="${WD}/benchmark_results"

DCS_REF="/home/abradbur/00_nesi_projects/00_nesi_projects/massey04083_nobackup/poultry_processing_chapter2/ASC_BFBB/DCS.fasta"
HUMAN_REF="/home/abradbur/00_nesi_projects/00_nesi_projects/massey04083_nobackup/poultry_processing_chapter2/ASC_BFBB/GCF_000001405.40_GRCh38.p14_genomic.fna"
CHICKEN_REF="/home/abradbur/00_nesi_projects/00_nesi_projects/massey04083_nobackup/poultry_processing_chapter2/ASC_BFBB/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna"

# Kraken2 database — update to your database on NeSI.
# To find available databases run: find /opt/nesi/db/ -name "*.k2d" 2>/dev/null | head
# IMPORTANT: The standard Kraken2 database (RefSeq bacteria/viral/human) does
# NOT include Gallus gallus. kraken_pct_chicken will be 0 unless your database
# was built with vertebrate genomes (custom build including chicken RefSeq).
# Check with: kraken2-inspect --db $KRAKEN2_DB | grep -i "gallus"
# If chicken is absent, rely solely on blast_pct_chicken for host assessment.
KRAKEN2_DB="/opt/nesi/db/kraken2/standard"   # EDIT THIS — or set RUN_KRAKEN2=false

# Set to false to skip Kraken2 (useful if database path is unknown or unavailable).
RUN_KRAKEN2=false

# Subsample each sample to this many reads before host removal.
# Keeps benchmark runtime manageable; use 0 for all reads (very slow).
SUBSAMPLE_READS=500000

# Number of reads to BLAST against chicken genome for assessment.
# BLAST is slow on long ONT reads — keep this small (5k–20k is sufficient).
BLAST_SUBSAMPLE=10000

# Set to false to skip BLAST assessment (Kraken2 only — much faster).
RUN_BLAST=true

THREADS=${SLURM_CPUS_PER_TASK:-16}

# ============================================================================
# PARAMETER SETS
# Format: "label|tool|extra_flags"
#
# tool = minimap2 or bowtie2
#
# --- minimap2 flags (applied on top of -x map-ont base preset) ---
#   -k INT   k-mer size (map-ont default: 15)
#            Lower k → more sensitive, catches more divergent host reads,
#            but slower and more likely to tag pathogen reads as host.
#   -s INT   minimum DP alignment score (map-ont default: 100)
#            Lower s → accepts weaker alignments as host-mapped.
#
# --- bowtie2 flags ---
#   --local               local alignment mode (better for long reads than end-to-end)
#   --sensitive-local     preset: -D 15 -R 2 -N 0 -L 20 -i S,1,0.75
#   --very-sensitive-local preset: -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
#   --score-min L,0,1.0   custom score minimum (lower = accepts weaker hits)
# ============================================================================

PARAM_SETS=(
    # minimap2 — baseline and sensitivity sweep
    "mm2_default|minimap2|"
    "mm2_k12|minimap2|-k 12"
    "mm2_k10|minimap2|-k 10"
    "mm2_k12_s60|minimap2|-k 12 -s 60"
    "mm2_k10_s40|minimap2|-k 10 -s 40"
    # bowtie2 — local alignment, three sensitivity presets
    "bt2_local|bowtie2|--local"
    "bt2_sensitive_local|bowtie2|--sensitive-local"
    "bt2_very_sensitive_local|bowtie2|--very-sensitive-local"
)

# ============================================================================
# SETUP
# ============================================================================

mkdir -p "${OUTDIR}" logs

log()        { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"; }
check_file() { [[ -f "$1" ]] || { echo "ERROR: File not found: $1"; exit 1; }; }
check_dir()  { [[ -d "$1" ]] || { echo "ERROR: Directory not found: $1"; exit 1; }; }

check_dir  "${INPUT_DIR}"
check_file "${DCS_REF}"
check_file "${HUMAN_REF}"
check_file "${CHICKEN_REF}"
[[ "${RUN_KRAKEN2}" == "true" ]] && check_dir "${KRAKEN2_DB}"

FASTQ_FILES=( $(ls "${INPUT_DIR}"/*.fastq.gz "${INPUT_DIR}"/*.fq.gz 2>/dev/null || true) )
NUM_FILES=${#FASTQ_FILES[@]}
[[ ${NUM_FILES} -gt 0 ]] || { echo "ERROR: No FASTQ files in ${INPUT_DIR}"; exit 1; }
log "Found ${NUM_FILES} FASTQ files"
log "Testing ${#PARAM_SETS[@]} parameter sets"
[[ ${SUBSAMPLE_READS} -gt 0 ]] && log "Subsampling to ${SUBSAMPLE_READS} reads per sample"

# ============================================================================
# BUILD BOWTIE2 INDICES (once — param-independent, reused for all bt2 sets)
# Building for large genomes (human ~3.1 GB) requires --large-index and can
# take 1–2 hours; indices are cached and skipped if already present.
# ============================================================================

BT2_INDEX_DIR="${OUTDIR}/bowtie2_indices"
mkdir -p "${BT2_INDEX_DIR}"

BT2_DCS_IDX="${BT2_INDEX_DIR}/dcs"
BT2_HUMAN_IDX="${BT2_INDEX_DIR}/human"
BT2_CHICKEN_IDX="${BT2_INDEX_DIR}/chicken"

build_bt2_index() {
    local REF="$1"
    local IDX="$2"
    local LABEL="$3"
    if [[ ! -f "${IDX}.1.bt2" && ! -f "${IDX}.1.bt2l" ]]; then
        log "Building Bowtie2 index: ${LABEL}"
        bowtie2-build --threads "${THREADS}" --large-index "${REF}" "${IDX}"
        log "Bowtie2 index ready: ${LABEL}"
    else
        log "Bowtie2 index already exists, skipping: ${LABEL}"
    fi
}

# Only build if any bowtie2 param sets are defined
if printf '%s\n' "${PARAM_SETS[@]}" | grep -q '|bowtie2|'; then
    build_bt2_index "${DCS_REF}"     "${BT2_DCS_IDX}"     "DCS"
    build_bt2_index "${HUMAN_REF}"   "${BT2_HUMAN_IDX}"   "human"
    build_bt2_index "${CHICKEN_REF}" "${BT2_CHICKEN_IDX}" "chicken"
fi

# ============================================================================
# BUILD CHICKEN BLAST DATABASE (once, reused across all parameter sets)
# ============================================================================

if [[ "${RUN_BLAST}" == "true" ]]; then
    CHICKEN_BLASTDB="${OUTDIR}/chicken_blastdb/chicken"
    mkdir -p "${OUTDIR}/chicken_blastdb"
    if [[ ! -f "${CHICKEN_BLASTDB}.nhr" ]]; then
        log "Building BLAST database from chicken reference (one-time setup)"
        makeblastdb \
            -in "${CHICKEN_REF}" \
            -dbtype nucl \
            -out "${CHICKEN_BLASTDB}" \
            -parse_seqids \
            -title "GalGal_broiler"
        log "BLAST database ready"
    fi
fi

# ============================================================================
# FUNCTIONS
# ============================================================================

# minimap2: map reads to reference; write unmapped reads to OUT_FASTQ.
# Echoes count of primary-mapped (host) reads removed.
filter_minimap2() {
    local EXTRA_FLAGS="$1"
    local REF="$2"
    local IN_FASTQ="$3"
    local OUT_FASTQ="$4"
    local TMP_BAM="${WORK_TMPDIR}/filter_$$.bam"

    # shellcheck disable=SC2086
    minimap2 -t "${THREADS}" -x map-ont --secondary=no "${EXTRA_FLAGS}" -a \
        "${REF}" "${IN_FASTQ}" 2>/dev/null \
        | samtools view -@ "${THREADS}" -bS - > "${TMP_BAM}"

    # Count primary mapped reads (exclude unmapped=4, secondary=256, supplementary=2048)
    local MAPPED
    MAPPED=$(samtools view -@ "${THREADS}" -c -F 2308 "${TMP_BAM}")

    # Keep only unmapped reads for the next stage
    samtools view -@ "${THREADS}" -b -f 4 "${TMP_BAM}" \
        | samtools fastq -@ "${THREADS}" - \
        | gzip > "${OUT_FASTQ}"

    rm -f "${TMP_BAM}"
    echo "${MAPPED}"
}

# Bowtie2: same logic, different aligner.
# Uses local alignment mode throughout (better for long ONT reads than end-to-end).
# Bowtie2 outputs both mapped and unmapped reads to SAM by default; we then
# extract unmapped with samtools flag 4, identical to the minimap2 approach.
filter_bowtie2() {
    local EXTRA_FLAGS="$1"
    local IDX="$2"
    local IN_FASTQ="$3"
    local OUT_FASTQ="$4"
    local TMP_BAM="${WORK_TMPDIR}/filter_$$.bam"

    # shellcheck disable=SC2086
    bowtie2 -x "${IDX}" \
        -U "${IN_FASTQ}" \
        -p "${THREADS}" \
        "${EXTRA_FLAGS}" \
        2>/dev/null \
        | samtools view -@ "${THREADS}" -bS - > "${TMP_BAM}"

    local MAPPED
    MAPPED=$(samtools view -@ "${THREADS}" -c -F 2308 "${TMP_BAM}")

    samtools view -@ "${THREADS}" -b -f 4 "${TMP_BAM}" \
        | samtools fastq -@ "${THREADS}" - \
        | gzip > "${OUT_FASTQ}"

    rm -f "${TMP_BAM}"
    echo "${MAPPED}"
}

# Dispatcher: calls the correct filter function based on tool name.
# Returns DCS/human/chicken refs or bt2 index prefixes depending on tool.
run_depletion_step() {
    local TOOL="$1"
    local EXTRA_FLAGS="$2"
    local REF_OR_IDX="$3"
    local IN_FASTQ="$4"
    local OUT_FASTQ="$5"

    if [[ "${TOOL}" == "minimap2" ]]; then
        filter_minimap2 "${EXTRA_FLAGS}" "${REF_OR_IDX}" "${IN_FASTQ}" "${OUT_FASTQ}"
    elif [[ "${TOOL}" == "bowtie2" ]]; then
        filter_bowtie2  "${EXTRA_FLAGS}" "${REF_OR_IDX}" "${IN_FASTQ}" "${OUT_FASTQ}"
    else
        echo "ERROR: Unknown tool '${TOOL}'" >&2; exit 1
    fi
}

# Run Kraken2 and parse the report for key taxa.
# Echoes three tab-separated values: pct_unclassified, pct_human, pct_chicken
run_kraken2() {
    local IN_FASTQ="$1"
    local OUT_PREFIX="$2"

    kraken2 \
        --db "${KRAKEN2_DB}" \
        --threads "${THREADS}" \
        --gzip-compressed \
        --report "${OUT_PREFIX}.report" \
        --output "${OUT_PREFIX}.kraken" \
        "${IN_FASTQ}" > /dev/null 2>&1

    # Kraken2 report columns: pct  clade_reads  taxon_reads  rank_code  taxid  name
    local PCT_UNCLASSIFIED PCT_HUMAN PCT_CHICKEN
    PCT_UNCLASSIFIED=$(awk '$4 == "U"  { printf "%.2f", $1 }' "${OUT_PREFIX}.report")
    PCT_HUMAN=$(       awk '$5 == 9606 { printf "%.2f", $1 }' "${OUT_PREFIX}.report")
    PCT_CHICKEN=$(     awk '$5 == 9031 { printf "%.2f", $1 }' "${OUT_PREFIX}.report")

    echo "${PCT_UNCLASSIFIED:-0.00}	${PCT_HUMAN:-0.00}	${PCT_CHICKEN:-0.00}"
}

# BLAST a subsample of clean reads against the chicken BLAST database.
# Echoes % of reads with at least one hit (evalue <= 1e-5).
blast_vs_chicken() {
    local IN_FASTQ="$1"
    local OUT_PREFIX="$2"
    local N_READS="$3"

    local FASTA="${OUT_PREFIX}_blast_input.fasta"
    local BLAST_OUT="${OUT_PREFIX}_blast_chicken.tsv"

    # Convert to FASTA with fixed seed — same reads tested across all param sets
    seqtk sample -s 42 "${IN_FASTQ}" "${N_READS}" 2>/dev/null \
        | awk 'NR%4==1{ print ">" substr($0,2) } NR%4==2{ print }' > "${FASTA}"

    local QUERY_COUNT
    QUERY_COUNT=$(grep -c "^>" "${FASTA}" || echo 0)

    if [[ ${QUERY_COUNT} -eq 0 ]]; then
        rm -f "${FASTA}"
        echo "0.00"
        return
    fi

    blastn \
        -query "${FASTA}" \
        -db "${CHICKEN_BLASTDB}" \
        -out "${BLAST_OUT}" \
        -outfmt "6 qseqid sseqid pident length evalue bitscore" \
        -max_target_seqs 1 \
        -max_hsps 1 \
        -num_threads "${THREADS}" \
        -evalue 1e-5 \
        2>/dev/null

    local HIT_COUNT
    HIT_COUNT=$(awk '!seen[$1]++{ c++ } END{ print c+0 }' "${BLAST_OUT}")

    awk -v hits="${HIT_COUNT}" -v total="${QUERY_COUNT}" \
        'BEGIN { printf "%.2f\n", (hits * 100) / total }' /dev/null

    rm -f "${FASTA}"
}

# ============================================================================
# SUMMARY TABLE HEADER
# ============================================================================

SUMMARY_TSV="${OUTDIR}/benchmark_summary.tsv"
{
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" \
        "param_set" "tool" "sample" \
        "total_reads" "clean_reads" "pct_retained" \
        "dcs_removed" "human_removed" "chicken_removed"
    [[ "${RUN_KRAKEN2}" == "true" ]] && printf "\t%s\t%s\t%s" \
        "kraken_pct_unclassified" "kraken_pct_human" "kraken_pct_chicken"
    [[ "${RUN_BLAST}" == "true" ]] && printf "\t%s" "blast_pct_chicken"
    printf "\n"
} > "${SUMMARY_TSV}"

# ============================================================================
# MAIN BENCHMARK LOOP
# ============================================================================

for PARAM_ENTRY in "${PARAM_SETS[@]}"; do
    IFS='|' read -r LABEL TOOL EXTRA_FLAGS <<< "${PARAM_ENTRY}"

    log "=============================================="
    log "Parameter set: [${LABEL}]  tool: ${TOOL}  flags: '${EXTRA_FLAGS:-<defaults>}'"
    log "=============================================="

    PARAM_OUTDIR="${OUTDIR}/${LABEL}"
    WORK_TMPDIR="${PARAM_OUTDIR}/tmp"
    mkdir -p "${PARAM_OUTDIR}" "${WORK_TMPDIR}"

    # Set reference paths based on tool
    if [[ "${TOOL}" == "minimap2" ]]; then
        DCS_REF_OR_IDX="${DCS_REF}"
        HUMAN_REF_OR_IDX="${HUMAN_REF}"
        CHICKEN_REF_OR_IDX="${CHICKEN_REF}"
    else
        DCS_REF_OR_IDX="${BT2_DCS_IDX}"
        HUMAN_REF_OR_IDX="${BT2_HUMAN_IDX}"
        CHICKEN_REF_OR_IDX="${BT2_CHICKEN_IDX}"
    fi

    for INPUT_FASTQ in "${FASTQ_FILES[@]}"; do
        BASENAME=$(basename "${INPUT_FASTQ}" .fastq.gz)
        BASENAME=$(basename "${BASENAME}" .fq.gz)

        log "[${LABEL}/${BASENAME}] ---- starting ----"

        SAMPLE_TMPDIR="${WORK_TMPDIR}/${BASENAME}"
        mkdir -p "${SAMPLE_TMPDIR}"

        # ------------------------------------------------------------------
        # Subsample input (fixed seed = identical reads across all param sets)
        # ------------------------------------------------------------------
        if [[ ${SUBSAMPLE_READS} -gt 0 ]]; then
            WORKING_FASTQ="${SAMPLE_TMPDIR}/${BASENAME}_sub.fastq.gz"
            log "[${LABEL}/${BASENAME}] Subsampling to ${SUBSAMPLE_READS} reads"
            seqtk sample -s 42 "${INPUT_FASTQ}" "${SUBSAMPLE_READS}" 2>/dev/null \
                | gzip > "${WORKING_FASTQ}"
        else
            WORKING_FASTQ="${INPUT_FASTQ}"
        fi

        TOTAL_READS=$(zcat "${WORKING_FASTQ}" | awk 'NR%4==1' | wc -l)
        log "[${LABEL}/${BASENAME}] Input reads: ${TOTAL_READS}"

        # ------------------------------------------------------------------
        # Sequential host removal: DCS → human → chicken
        # Each step checks the output is non-empty before continuing so a
        # silent failure (e.g. OOM-killed aligner) causes a loud abort
        # rather than cascading zeros through the rest of the pipeline.
        # ------------------------------------------------------------------
        check_reads() {
            local FILE="$1" LABEL="$2"
            local N
            N=$(zcat "${FILE}" | awk 'NR%4==1' | wc -l)
            [[ ${N} -gt 0 ]] || { echo "ERROR: ${LABEL} produced 0 reads — aligner may have been OOM-killed. Increase --mem."; exit 1; }
            echo "${N}"
        }

        AFTER_DCS="${SAMPLE_TMPDIR}/${BASENAME}_after_dcs.fastq.gz"
        DCS_REMOVED=$(run_depletion_step "${TOOL}" "${EXTRA_FLAGS}" \
            "${DCS_REF_OR_IDX}" "${WORKING_FASTQ}" "${AFTER_DCS}")
        check_reads "${AFTER_DCS}" "DCS step" > /dev/null
        log "[${LABEL}/${BASENAME}] DCS removed: ${DCS_REMOVED}"

        AFTER_HUMAN="${SAMPLE_TMPDIR}/${BASENAME}_after_human.fastq.gz"
        HUMAN_REMOVED=$(run_depletion_step "${TOOL}" "${EXTRA_FLAGS}" \
            "${HUMAN_REF_OR_IDX}" "${AFTER_DCS}" "${AFTER_HUMAN}")
        check_reads "${AFTER_HUMAN}" "human step" > /dev/null
        log "[${LABEL}/${BASENAME}] Human removed: ${HUMAN_REMOVED}"

        CLEAN_FASTQ="${PARAM_OUTDIR}/${BASENAME}_filtered.fastq.gz"
        CHICKEN_REMOVED=$(run_depletion_step "${TOOL}" "${EXTRA_FLAGS}" \
            "${CHICKEN_REF_OR_IDX}" "${AFTER_HUMAN}" "${CLEAN_FASTQ}")
        log "[${LABEL}/${BASENAME}] Chicken removed: ${CHICKEN_REMOVED}"

        CLEAN_READS=$(zcat "${CLEAN_FASTQ}" | awk 'NR%4==1' | wc -l)
        PCT_RETAINED=$(echo "scale=2; ${CLEAN_READS} * 100 / ${TOTAL_READS}" | bc)
        log "[${LABEL}/${BASENAME}] Retained: ${CLEAN_READS}/${TOTAL_READS} (${PCT_RETAINED}%)"

        # ------------------------------------------------------------------
        # Assessment 1: Kraken2 taxonomic classification (optional)
        # ------------------------------------------------------------------
        K_UNCLASS="N/A"; K_HUMAN="N/A"; K_CHICKEN="N/A"
        if [[ "${RUN_KRAKEN2}" == "true" ]]; then
            log "[${LABEL}/${BASENAME}] Running Kraken2"
            KRAKEN_PREFIX="${PARAM_OUTDIR}/${BASENAME}_kraken2"
            KRAKEN_STATS=$(run_kraken2 "${CLEAN_FASTQ}" "${KRAKEN_PREFIX}")
            IFS=$'\t' read -r K_UNCLASS K_HUMAN K_CHICKEN <<< "${KRAKEN_STATS}"
            log "[${LABEL}/${BASENAME}] Kraken2: unclassified=${K_UNCLASS}%  human=${K_HUMAN}%  chicken=${K_CHICKEN}%"
        else
            log "[${LABEL}/${BASENAME}] Kraken2 skipped (RUN_KRAKEN2=false)"
        fi

        # ------------------------------------------------------------------
        # Assessment 2: BLAST vs chicken genome (optional)
        # ------------------------------------------------------------------
        BLAST_CHICKEN="N/A"
        if [[ "${RUN_BLAST}" == "true" ]]; then
            log "[${LABEL}/${BASENAME}] Running BLAST vs chicken (n=${BLAST_SUBSAMPLE} reads)"
            BLAST_CHICKEN=$(blast_vs_chicken "${CLEAN_FASTQ}" \
                "${PARAM_OUTDIR}/${BASENAME}" "${BLAST_SUBSAMPLE}")
            log "[${LABEL}/${BASENAME}] BLAST chicken hits: ${BLAST_CHICKEN}%"
        fi

        # ------------------------------------------------------------------
        # Append row to summary table
        # ------------------------------------------------------------------
        {
            printf "%s\t%s\t%s\t%d\t%d\t%s\t%d\t%d\t%d" \
                "${LABEL}" "${TOOL}" "${BASENAME}" \
                "${TOTAL_READS}" "${CLEAN_READS}" "${PCT_RETAINED}" \
                "${DCS_REMOVED}" "${HUMAN_REMOVED}" "${CHICKEN_REMOVED}"
            [[ "${RUN_KRAKEN2}" == "true" ]] && printf "\t%s\t%s\t%s" \
                "${K_UNCLASS}" "${K_HUMAN}" "${K_CHICKEN}"
            [[ "${RUN_BLAST}" == "true" ]] && printf "\t%s" "${BLAST_CHICKEN}"
            printf "\n"
        } >> "${SUMMARY_TSV}"

        rm -rf "${SAMPLE_TMPDIR}"
        log "[${LABEL}/${BASENAME}] ---- done ----"
    done

    rm -rf "${WORK_TMPDIR}"
done

# ============================================================================
# FINAL REPORT
# ============================================================================

log "=============================================="
log "Benchmark complete — ${#PARAM_SETS[@]} parameter sets × ${NUM_FILES} samples"
log "Summary: ${SUMMARY_TSV}"
log "=============================================="
echo ""
echo "================================================================================================"
echo "BENCHMARK SUMMARY"
echo "================================================================================================"
column -t -s $'\t' "${SUMMARY_TSV}"
echo ""
echo "Interpretation guide:"
echo "  pct_retained              — lower = more aggressive removal (risk of losing pathogen reads)"
echo "  chicken_removed           — reads removed at the chicken alignment step"
echo "  kraken_pct_chicken        — residual chicken by taxonomy (0 if chicken not in your db)"
echo "  blast_pct_chicken         — residual chicken by direct alignment; primary metric for host removal"
echo "  kraken_pct_unclassified   — reads with no taxonomic assignment; high is expected/normal"
echo "  kraken_pct_human          — residual human; should approach 0 after depletion"
echo ""
echo "Goal: minimise blast_pct_chicken while keeping pct_retained close to mm2_default."
echo "Bowtie2 results will typically show higher residual chicken — this is expected"
echo "because Bowtie2 is not designed for long Nanopore reads."
