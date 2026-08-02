#!/usr/bin/env bash
# =============================================================================
# rRNA SNP Resistance Pipeline — Nanopore Metagenomics
# =============================================================================
# Input:  pre-basecalled, QC-filtered FASTQ reads
# Output: per-sample VCF annotated with resistance SNPs, coverage stats,
#         and a summary TSV of resistance calls
#
# Dependencies (must be on PATH):
#   minimap2, samtools, medaka (or clair3), bgzip, tabix, bcftools
#   Optional: bedtools (for coverage QC), snpEff/snpSift (for annotation)
#
# Usage:
#   bash rrna_snp_pipeline.sh [options]
#
# Options:
#   -i  Input FASTQ (or directory of FASTQs, will be merged)   [required]
#   -r  Reference FASTA (16S + 23S rRNA sequences)             [required]
#   -o  Output directory                                        [default: ./results]
#   -s  Sample name                                             [default: sample]
#   -b  Resistance SNP BED file (chrom pos ref alt drug)       [optional]
#   -t  Threads                                                 [default: 4]
#   -m  Medaka model (e.g. r941_min_hac_g507)                  [default: r1041_e82_400bps_sup_v5.0.0]
#   -q  Minimum mapping quality                                 [default: 20]
#   -d  Minimum read depth for SNP calling                     [default: 20]
#   -f  Minimum allele frequency for variant reporting          [default: 0.05]
#   -c  Variant caller: medaka | clair3 | longshot              [default: medaka]
#   -h  Show this help
# =============================================================================

set -euo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────────
INPUT=""
REFERENCE=""
OUTDIR="./results"
SAMPLE="sample"
SNP_BED=""
THREADS=4
MEDAKA_MODEL="r1041_e82_400bps_sup_v5.0.0"
MIN_MAPQ=20
MIN_DEPTH=20
MIN_AF=0.05
CALLER="medaka"

# ── Colours for logging ───────────────────────────────────────────────────────
RED='\033[0;31m'; YELLOW='\033[1;33m'; GREEN='\033[0;32m'
CYAN='\033[0;36m'; NC='\033[0m'

log()  { echo -e "${CYAN}[$(date '+%H:%M:%S')] $*${NC}"; }
warn() { echo -e "${YELLOW}[WARN] $*${NC}"; }
err()  { echo -e "${RED}[ERROR] $*${NC}"; exit 1; }
ok()   { echo -e "${GREEN}[OK]   $*${NC}"; }

# ── Argument parsing ──────────────────────────────────────────────────────────
usage() {
  grep '^#' "$0" | grep -v '#!/' | sed 's/^# \?//'
  exit 0
}

while getopts ":i:r:o:s:b:t:m:q:d:f:c:h" opt; do
  case $opt in
    i) INPUT="$OPTARG" ;;
    r) REFERENCE="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    s) SAMPLE="$OPTARG" ;;
    b) SNP_BED="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    m) MEDAKA_MODEL="$OPTARG" ;;
    q) MIN_MAPQ="$OPTARG" ;;
    d) MIN_DEPTH="$OPTARG" ;;
    f) MIN_AF="$OPTARG" ;;
    c) CALLER="$OPTARG" ;;
    h) usage ;;
    :) err "Option -$OPTARG requires an argument." ;;
    \?) err "Unknown option: -$OPTARG" ;;
  esac
done

# ── Validate required args ────────────────────────────────────────────────────
[[ -z "$INPUT"     ]] && err "Input FASTQ/directory required (-i)"
[[ -z "$REFERENCE" ]] && err "Reference FASTA required (-r)"
[[ "$CALLER" =~ ^(medaka|clair3|longshot)$ ]] || err "Caller must be medaka, clair3, or longshot"

# ── Dependency checks ─────────────────────────────────────────────────────────
check_dep() {
  command -v "$1" &>/dev/null || err "Dependency not found: $1. Please install or add to PATH."
}

check_dep minimap2
check_dep samtools
check_dep bgzip
check_dep tabix
check_dep bcftools

case "$CALLER" in
  medaka)   check_dep medaka ;;
  clair3)   check_dep run_clair3.sh ;;
  longshot) check_dep longshot ;;
esac

# ── Setup output directories ──────────────────────────────────────────────────
LOGDIR="${OUTDIR}/logs"
ALIGNDIR="${OUTDIR}/alignments"
VCFDIR="${OUTDIR}/variants"
REPORTDIR="${OUTDIR}/report"

mkdir -p "$LOGDIR" "$ALIGNDIR" "$VCFDIR" "$REPORTDIR"

LOGFILE="${LOGDIR}/${SAMPLE}_pipeline.log"
exec > >(tee -a "$LOGFILE") 2>&1

log "=================================================="
log "  rRNA SNP Resistance Pipeline"
log "  Sample:    $SAMPLE"
log "  Input:     $INPUT"
log "  Reference: $REFERENCE"
log "  Caller:    $CALLER"
log "  Threads:   $THREADS"
log "  MinMapQ:   $MIN_MAPQ | MinDepth: $MIN_DEPTH | MinAF: $MIN_AF"
log "=================================================="

# =============================================================================
# STEP 1: Merge input FASTQs if a directory is given
# =============================================================================
log "STEP 1: Preparing input reads"

MERGED_FASTQ="${OUTDIR}/${SAMPLE}_merged.fastq.gz"

if [[ -d "$INPUT" ]]; then
  log "  Merging FASTQs from directory: $INPUT"
  FQCOUNT=$(find "$INPUT" -name "*.fastq.gz" -o -name "*.fq.gz" -o -name "*.fastq" -o -name "*.fq" | wc -l)
  [[ "$FQCOUNT" -eq 0 ]] && err "No FASTQ files found in $INPUT"
  log "  Found $FQCOUNT FASTQ file(s)"

  # Merge: handle both gzipped and plain
  find "$INPUT" -name "*.fastq.gz" -o -name "*.fq.gz" | sort | xargs -I{} zcat {} 2>/dev/null | \
    gzip > "$MERGED_FASTQ" || true
  find "$INPUT" -name "*.fastq" -o -name "*.fq" | sort | xargs -I{} cat {} 2>/dev/null | \
    gzip >> "$MERGED_FASTQ" || true

elif [[ -f "$INPUT" ]]; then
  if [[ "$INPUT" == *.gz ]]; then
    ln -sf "$(realpath "$INPUT")" "$MERGED_FASTQ"
  else
    gzip -c "$INPUT" > "$MERGED_FASTQ"
  fi
else
  err "Input not found: $INPUT"
fi

TOTAL_READS=$(zcat "$MERGED_FASTQ" | awk 'NR%4==1' | wc -l)
log "  Total reads to process: $TOTAL_READS"
ok "STEP 1 complete"

# =============================================================================
# STEP 2: Index reference (if not already indexed)
# =============================================================================
log "STEP 2: Indexing reference"

[[ -f "$REFERENCE" ]] || err "Reference FASTA not found: $REFERENCE"

if [[ ! -f "${REFERENCE}.fai" ]]; then
  log "  Indexing with samtools faidx"
  samtools faidx "$REFERENCE"
fi

# Minimap2 index (speeds up repeated mapping runs)
MM2_INDEX="${OUTDIR}/ref.mmi"
if [[ ! -f "$MM2_INDEX" ]]; then
  log "  Building Minimap2 index"
  minimap2 -d "$MM2_INDEX" "$REFERENCE" -t "$THREADS" \
    2>> "${LOGDIR}/minimap2_index.log"
fi

ok "STEP 2 complete"

# =============================================================================
# STEP 3: Align reads to rRNA reference
# =============================================================================
log "STEP 3: Aligning reads with Minimap2 (map-ont)"

RAW_BAM="${ALIGNDIR}/${SAMPLE}_raw.bam"
FILTERED_BAM="${ALIGNDIR}/${SAMPLE}_filtered.bam"
FINAL_BAM="${ALIGNDIR}/${SAMPLE}.bam"

minimap2 \
  -ax map-ont \
  -t "$THREADS" \
  --secondary=no \
  -Y \
  "$MM2_INDEX" \
  "$MERGED_FASTQ" \
  2>> "${LOGDIR}/minimap2_align.log" \
| samtools sort \
  -@ "$((THREADS > 2 ? THREADS - 2 : 1))" \
  -o "$RAW_BAM"

samtools index -@ "$THREADS" "$RAW_BAM"

MAPPED=$(samtools view -c -F 4 "$RAW_BAM")
UNMAPPED=$(samtools view -c -f 4 "$RAW_BAM")
log "  Mapped reads:   $MAPPED"
log "  Unmapped reads: $UNMAPPED"

if [[ "$MAPPED" -eq 0 ]]; then
  err "No reads mapped. Check that your reference contains the correct rRNA sequences."
fi

# ── Filter: mapQ, primary alignments only ─────────────────────────────────
log "  Filtering: mapQ≥${MIN_MAPQ}, primary alignments, removing supplementary"

samtools view \
  -@ "$THREADS" \
  -b \
  -q "$MIN_MAPQ" \
  -F 2308 \
  "$RAW_BAM" \
| samtools sort -@ "$THREADS" -o "$FILTERED_BAM"

samtools index -@ "$THREADS" "$FILTERED_BAM"

FILTERED=$(samtools view -c -F 4 "$FILTERED_BAM")
log "  Reads after quality filter: $FILTERED"

# ── Mark duplicates (flag only, don't remove — rRNA is multicopy) ──────────
# We flag but keep duplicates since rRNA paralogs are legitimate multicopy genes.
# For strict analysis, set -r to also remove them, but this will reduce depth.
samtools markdup \
  -@ "$THREADS" \
  --no-PG \
  "$FILTERED_BAM" \
  "$FINAL_BAM" \
  2>> "${LOGDIR}/markdup.log"

samtools index -@ "$THREADS" "$FINAL_BAM"

ok "STEP 3 complete — final BAM: $FINAL_BAM"

# =============================================================================
# STEP 4: Coverage QC
# =============================================================================
log "STEP 4: Coverage statistics"

COV_FILE="${REPORTDIR}/${SAMPLE}_coverage.tsv"

samtools coverage "$FINAL_BAM" > "$COV_FILE"

# Per-reference coverage summary
log "  Reference coverage summary:"
awk 'NR>1 {
  printf "    %-40s  depth=%.1f  breadth=%.1f%%  len=%d\n", $1, $7, $6, $3
}' "$COV_FILE"

# Warn if any reference has low depth
LOW_DEPTH=$(awk -v d="$MIN_DEPTH" 'NR>1 && $7 < d {count++} END {print count+0}' "$COV_FILE")
if [[ "$LOW_DEPTH" -gt 0 ]]; then
  warn "$LOW_DEPTH reference(s) have mean depth < ${MIN_DEPTH}×. SNP calls may be unreliable."
fi

# Per-position depth (useful for downstream filtering)
samtools depth -a "$FINAL_BAM" > "${REPORTDIR}/${SAMPLE}_depth.tsv"

ok "STEP 4 complete — coverage stats: $COV_FILE"

# =============================================================================
# STEP 5: Variant calling
# =============================================================================
log "STEP 5: Variant calling with $CALLER"

VCF_RAW="${VCFDIR}/${SAMPLE}_raw.vcf.gz"
VCF_FILTERED="${VCFDIR}/${SAMPLE}_filtered.vcf.gz"

case "$CALLER" in

  # ── Medaka ────────────────────────────────────────────────────────────────
  medaka)
    log "  Running medaka haploid variant"
    MEDAKA_OUTDIR="${VCFDIR}/medaka_${SAMPLE}"
    mkdir -p "$MEDAKA_OUTDIR"

    medaka haploid_variant \
      -i "$FINAL_BAM" \
      -r "$REFERENCE" \
      -o "$MEDAKA_OUTDIR" \
      -m "$MEDAKA_MODEL" \
      -t "$THREADS" \
      2>> "${LOGDIR}/medaka.log"

    # Medaka outputs medaka.annotated.vcf — copy and compress
    MEDAKA_VCF=$(find "$MEDAKA_OUTDIR" -name "*.annotated.vcf" | head -1)
    [[ -z "$MEDAKA_VCF" ]] && MEDAKA_VCF=$(find "$MEDAKA_OUTDIR" -name "*.vcf" | head -1)
    [[ -z "$MEDAKA_VCF" ]] && err "Medaka did not produce a VCF. Check ${LOGDIR}/medaka.log"

    bgzip -c "$MEDAKA_VCF" > "$VCF_RAW"
    tabix -p vcf "$VCF_RAW"
    ;;

  # ── Clair3 ────────────────────────────────────────────────────────────────
  clair3)
    log "  Running Clair3"
    # Clair3 requires a model directory. Set CLAIR3_MODEL or pass via env.
    CLAIR3_MODEL="${CLAIR3_MODEL:-/opt/models/ont}"
    [[ -d "$CLAIR3_MODEL" ]] || err "Clair3 model directory not found: $CLAIR3_MODEL (set CLAIR3_MODEL env var)"

    CLAIR3_OUTDIR="${VCFDIR}/clair3_${SAMPLE}"
    mkdir -p "$CLAIR3_OUTDIR"

    run_clair3.sh \
      --bam_fn="$FINAL_BAM" \
      --ref_fn="$REFERENCE" \
      --threads="$THREADS" \
      --platform="ont" \
      --model_path="$CLAIR3_MODEL" \
      --output="$CLAIR3_OUTDIR" \
      --haploid_sensitive \
      --min_coverage="$MIN_DEPTH" \
      2>> "${LOGDIR}/clair3.log"

    CLAIR3_VCF="${CLAIR3_OUTDIR}/merge_output.vcf.gz"
    [[ -f "$CLAIR3_VCF" ]] || err "Clair3 did not produce output. Check ${LOGDIR}/clair3.log"
    cp "$CLAIR3_VCF" "$VCF_RAW"
    cp "${CLAIR3_VCF}.tbi" "${VCF_RAW}.tbi"
    ;;

  # ── Longshot ──────────────────────────────────────────────────────────────
  longshot)
    log "  Running Longshot"
    VCF_LONGSHOT="${VCFDIR}/${SAMPLE}_longshot.vcf"

    longshot \
      --bam "$FINAL_BAM" \
      --ref "$REFERENCE" \
      --out "$VCF_LONGSHOT" \
      --min_alt_frac "$MIN_AF" \
      --min_cov "$MIN_DEPTH" \
      --strand_bias_pvalue_cutoff 0.01 \
      2>> "${LOGDIR}/longshot.log"

    bgzip -c "$VCF_LONGSHOT" > "$VCF_RAW"
    tabix -p vcf "$VCF_RAW"
    ;;
esac

TOTAL_VARIANTS=$(bcftools view -H "$VCF_RAW" | wc -l)
log "  Raw variants called: $TOTAL_VARIANTS"

# =============================================================================
# STEP 6: VCF filtering — depth, allele frequency, SNPs only
# =============================================================================
log "STEP 6: Filtering variants"

# Keep only SNPs (no indels — rRNA resistance is SNP-mediated)
# Apply minimum depth and allele frequency filters
bcftools view \
  --type snps \
  "$VCF_RAW" \
| bcftools filter \
  --include "INFO/DP >= ${MIN_DEPTH}" \
  --soft-filter "LowDepth" \
  - \
| bcftools filter \
  --include "MIN(FORMAT/AF) >= ${MIN_AF} || MIN(FORMAT/VAF) >= ${MIN_AF} || (AC/AN) >= ${MIN_AF}" \
  --soft-filter "LowAF" \
  - \
| bcftools view \
  --apply-filters "PASS,." \
  -O z -o "$VCF_FILTERED"

tabix -p vcf "$VCF_FILTERED"

FILTERED_VARIANTS=$(bcftools view -H "$VCF_FILTERED" | wc -l)
log "  Variants after filtering: $FILTERED_VARIANTS"

ok "STEP 6 complete — filtered VCF: $VCF_FILTERED"

# =============================================================================
# STEP 7: Resistance SNP annotation (if BED file provided)
# =============================================================================
log "STEP 7: Resistance SNP annotation"

RESISTANCE_TSV="${REPORTDIR}/${SAMPLE}_resistance_calls.tsv"

if [[ -n "$SNP_BED" && -f "$SNP_BED" ]]; then
  log "  Annotating against resistance BED: $SNP_BED"

  # Expected BED format (tab-separated):
  # CHROM  POS-1  POS  REF  ALT  DRUG_CLASS  GENE  DESCRIPTION
  # (POS-1 = 0-based start, POS = 1-based end for BED convention)
  # Example:
  #   23S_rRNA  2057  2058  A  G  Macrolides  23S  A2058G_clarithromycin_resistance

  # Intersect filtered VCF positions with resistance BED using bcftools
  {
    echo -e "SAMPLE\tCHROM\tPOS\tREF\tALT\tDEPTH\tAF\tDRUG_CLASS\tGENE\tDESCRIPTION\tFILTER"

    # Extract variant positions from VCF
    bcftools query \
      -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\t[%AF]\t%FILTER\n' \
      "$VCF_FILTERED" 2>/dev/null | \
    while IFS=$'\t' read -r chrom pos ref alt dp af filter; do
      # Look up this position in the resistance BED (BED POS = VCF POS for SNPs)
      match=$(awk -v c="$chrom" -v p="$pos" -v r="$ref" -v a="$alt" \
        'BEGIN{OFS="\t"} $1==c && $3==p && $4==r && $5==a {print $6, $7, $8}' \
        "$SNP_BED")

      if [[ -n "$match" ]]; then
        drug=$(echo "$match" | cut -f1)
        gene=$(echo "$match" | cut -f2)
        desc=$(echo "$match" | cut -f3)
        echo -e "${SAMPLE}\t${chrom}\t${pos}\t${ref}\t${alt}\t${dp:-NA}\t${af:-NA}\t${drug}\t${gene}\t${desc}\t${filter}"
      fi
    done
  } > "$RESISTANCE_TSV"

  HITS=$(tail -n +2 "$RESISTANCE_TSV" | wc -l)
  log "  Resistance SNPs detected: $HITS"

  if [[ "$HITS" -gt 0 ]]; then
    warn "RESISTANCE CALLS FOUND — see $RESISTANCE_TSV"
    echo ""
    echo "  ┌─ Resistance summary ─────────────────────────────────────────"
    tail -n +2 "$RESISTANCE_TSV" | awk -F'\t' \
      '{printf "  │  %-8s  %-5s→%-5s  AF=%-6s  %-20s  %s\n", $2":"$3, $4, $5, $7, $8, $10}'
    echo "  └───────────────────────────────────────────────────────────────"
    echo ""
  else
    log "  No known resistance SNPs detected at the applied thresholds."
  fi

else
  warn "No resistance BED file provided (-b). Skipping annotation."
  warn "Provide a BED file with columns: CHROM POS-1 POS REF ALT DRUG GENE DESC"
  warn "See README section below for known rRNA resistance positions."

  # Still produce an unannotated variant table for manual review
  {
    echo -e "SAMPLE\tCHROM\tPOS\tREF\tALT\tDEPTH\tAF\tFILTER"
    bcftools query \
      -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\t[%AF]\t%FILTER\n' \
      "$VCF_FILTERED" | \
    awk -v s="$SAMPLE" 'BEGIN{OFS="\t"} {print s, $0}'
  } > "$RESISTANCE_TSV"
fi

ok "STEP 7 complete — results: $RESISTANCE_TSV"

# =============================================================================
# STEP 8: Summary report
# =============================================================================
log "STEP 8: Writing summary report"

SUMMARY="${REPORTDIR}/${SAMPLE}_summary.txt"

{
  echo "========================================================"
  echo "  rRNA SNP Resistance Pipeline — Run Summary"
  echo "  Sample:  $SAMPLE"
  echo "  Date:    $(date)"
  echo "========================================================"
  echo ""
  echo "--- Input ---"
  echo "  Total reads:          $TOTAL_READS"
  echo "  Mapped reads:         $MAPPED"
  echo "  Post-filter reads:    $FILTERED"
  echo ""
  echo "--- Alignment ---"
  samtools flagstat "$FINAL_BAM" | sed 's/^/  /'
  echo ""
  echo "--- Coverage (per reference) ---"
  awk 'NR>1 {printf "  %-40s  mean_depth=%.1f  breadth=%.1f%%\n", $1, $7, $6}' "$COV_FILE"
  echo ""
  echo "--- Variant Calling ---"
  echo "  Caller:               $CALLER"
  echo "  Raw variants:         $TOTAL_VARIANTS"
  echo "  Filtered variants:    $FILTERED_VARIANTS"
  echo "  Min depth:            ${MIN_DEPTH}×"
  echo "  Min allele freq:      $MIN_AF"
  echo ""
  echo "--- Resistance Calls ---"
  if [[ -n "$SNP_BED" && -f "$SNP_BED" ]]; then
    HITS=$(tail -n +2 "$RESISTANCE_TSV" | wc -l)
    echo "  Resistance SNPs:      $HITS"
    if [[ "$HITS" -gt 0 ]]; then
      echo ""
      cat "$RESISTANCE_TSV"
    fi
  else
    echo "  (No resistance BED provided — see $RESISTANCE_TSV for all variants)"
  fi
  echo ""
  echo "--- Output Files ---"
  echo "  Alignments:   $FINAL_BAM"
  echo "  VCF (raw):    $VCF_RAW"
  echo "  VCF (filt):   $VCF_FILTERED"
  echo "  Coverage:     $COV_FILE"
  echo "  Results:      $RESISTANCE_TSV"
  echo "  Log:          $LOGFILE"
  echo ""
  echo "========================================================"
} > "$SUMMARY"

cat "$SUMMARY"

ok "Pipeline complete — summary: $SUMMARY"

# =============================================================================
# README — Known rRNA resistance SNPs (E. coli K-12 numbering)
# =============================================================================
# To use the annotation step, create a BED file (-b) like this:
#
# 23S_rRNA  2057  2058  A  G  Macrolides        23S  A2058G_clarithromycin_azithromycin
# 23S_rRNA  2058  2059  A  G  Macrolides        23S  A2059G_erythromycin
# 23S_rRNA  2575  2576  G  T  Oxazolidinones    23S  G2576U_linezolid
# 23S_rRNA  2450  2451  A  G  Oxazolidinones    23S  A2451G_linezolid
# 23S_rRNA  2503  2504  A  G  Amphenicols       23S  A2503G_chloramphenicol
# 23S_rRNA  2057  2058  A  T  Macrolides        23S  A2058T_macrolide
# 16S_rRNA  1407  1408  A  G  Aminoglycosides   16S  A1408G_gentamicin_tobramycin
# 16S_rRNA  1490  1491  G  T  Aminoglycosides   16S  G1491T_aminoglycosides
# 16S_rRNA  516   517   C  T  Tetracyclines     16S  C517T_tetracycline (some orgs)
#
# Adjust CHROM names to match your reference FASTA headers.
# Positions are 1-based in VCF but the BED file uses 0-based start / 1-based end.
# =============================================================================