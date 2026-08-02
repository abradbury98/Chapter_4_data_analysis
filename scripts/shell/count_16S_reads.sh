#!/bin/bash
# =============================================================================
# count_16S_reads.sh
#
# Map host-depleted ONT reads to a 16S rRNA reference database using minimap2.
# Outputs one BLAST-format TSV per sample — these feed directly into
# normalise_megares_reads.py as --sixteens_dir input.
#
# Why minimap2 and not BLASTN?
#   ONT reads are long (~1–10 kb) with higher error rates (~5–15 % raw).
#   BLASTN is optimised for short, accurate reads and is extremely slow on ONT.
#   minimap2 --preset map-ont handles variable-length, error-prone reads natively
#   and is orders of magnitude faster.
#
# 16S reference database
#   Recommended: SILVA SSU NR99 representative sequences (clustered at 99 %)
#   Download: https://www.arb-silva.de/download/arb-files/
#     → SILVA_138.2_SSURef_NR99_tax_silva_trunc.fasta.gz  (~180 MB)
#   Unzip and point REF_FASTA to the .fasta file.
#   On NeSI, check /nesi/project/nesi00000/databases/ or ask your support team.
#
# Output TSV columns (BLAST -outfmt 6 compatible + qlen + slen):
#   qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
#   pident  = (residue_matches / alignment_block_length) × 100
#   bitscore = residue_matches (higher = better; used for best-hit selection)
#   evalue   = 0  (not computed by minimap2; placeholder only)
#   mismatch = alignment_block_length − residue_matches (approximate)
#   gapopen  = 0  (not directly available from PAF; conservative placeholder)
#
# Filters applied (adjust below):
#   --min_identity  80 %  (lenient for 16S; ONT error rate ~5–15 %)
#   --min_coverage  80 %  of the 16S reference gene covered by the alignment
#
# Usage (interactive or SLURM):
#   bash count_16S_reads.sh \
#       --fastq_dir  /nesi/nobackup/.../host_removed_reads \
#       --ref        /nesi/nobackup/.../SILVA_138.2_SSURef_NR99.fasta \
#       --out_dir    /nesi/nobackup/.../16S_blast_tsvs \
#       --threads    16
#
# SLURM submission:
#   sbatch count_16S_reads.sh --fastq_dir ... --ref ... --out_dir ...
# =============================================================================

# ── SLURM directives (uncomment to submit as a job) ──────────────────────────
#SBATCH --job-name=count_16S
#SBATCH --account=YOUR_ACCOUNT        # replace with your NeSI project code
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=logs/count_16S_%j.out
#SBATCH --error=logs/count_16S_%j.err
# -----------------------------------------------------------------------------

set -euo pipefail

# ── Defaults ─────────────────────────────────────────────────────────────────
FASTQ_DIR=""
REF_FASTA=""
OUT_DIR="16S_blast_tsvs"
THREADS=8
MIN_IDENTITY=80      # minimum % identity to keep a 16S hit
MIN_COVERAGE=0.80    # minimum fraction of reference gene covered
SAMPLE_MAP=""        # optional TSV: barcode_stem <tab> sample_name

# ── Parse arguments ───────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --fastq_dir)    FASTQ_DIR="$2";    shift 2 ;;
        --ref)          REF_FASTA="$2";    shift 2 ;;
        --out_dir)      OUT_DIR="$2";      shift 2 ;;
        --threads)      THREADS="$2";      shift 2 ;;
        --min_identity) MIN_IDENTITY="$2"; shift 2 ;;
        --min_coverage) MIN_COVERAGE="$2"; shift 2 ;;
        --sample_map)   SAMPLE_MAP="$2";   shift 2 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

# ── Validate inputs ────────────────────────────────────────────────────────────
if [[ -z "$FASTQ_DIR" || -z "$REF_FASTA" ]]; then
    echo "ERROR: --fastq_dir and --ref are required."
    echo ""
    echo "Usage:"
    echo "  bash count_16S_reads.sh \\"
    echo "      --fastq_dir  /path/to/host_removed_reads \\"
    echo "      --ref        /path/to/SILVA_138.2_SSURef_NR99.fasta \\"
    echo "      --out_dir    /path/to/16S_blast_tsvs \\"
    echo "      --sample_map /path/to/host_removed_reads/sample_map.tsv \\"
    echo "      --threads    8"
    exit 1
fi

if [[ ! -d "$FASTQ_DIR" ]]; then
    echo "ERROR: FASTQ directory not found: $FASTQ_DIR"; exit 1
fi
if [[ ! -f "$REF_FASTA" ]]; then
    echo "ERROR: 16S reference FASTA not found: $REF_FASTA"
    echo ""
    echo "Download SILVA SSU NR99 from:"
    echo "  https://www.arb-silva.de/download/arb-files/"
    echo "  File: SILVA_138.2_SSURef_NR99_tax_silva_trunc.fasta.gz"
    echo "  Then: gunzip SILVA_138.2_SSURef_NR99_tax_silva_trunc.fasta.gz"
    exit 1
fi

# ── Load modules (NeSI) ────────────────────────────────────────────────────────
# Uncomment the lines appropriate to your NeSI environment.
# module purge
# module load minimap2/2.28-GCC-12.3.0
# module load SAMtools/1.19-GCC-12.3.0   # optional; only needed if using SAM output

# Verify minimap2 is available
if ! command -v minimap2 &>/dev/null; then
    echo "ERROR: minimap2 not found. Load it with 'module load minimap2' or install it."
    exit 1
fi

echo "minimap2 version: $(minimap2 --version)"
echo "Reference:  $REF_FASTA"
echo "FASTQ dir:  $FASTQ_DIR"
echo "Output dir: $OUT_DIR"
echo "Threads:    $THREADS"
echo "Min identity: ${MIN_IDENTITY}%  |  Min coverage: ${MIN_COVERAGE}"
[[ -n "$SAMPLE_MAP" ]] && echo "Sample map: $SAMPLE_MAP" || echo "Sample map: (none — using FASTQ stem as sample name)"
echo ""

# ── Sample map lookup (bash 3.2-compatible; uses awk instead of declare -A) ───
# lookup_sample STEM → prints sample name, or empty string if not found
lookup_sample() {
    local stem="$1"
    if [[ -z "$SAMPLE_MAP" ]]; then
        echo ""
        return
    fi
    awk -F'\t' -v stem="$stem" '
        NR == 1 { next }          # skip header
        $1 == stem { print $2; exit }
    ' "$SAMPLE_MAP"
}

if [[ -n "$SAMPLE_MAP" ]]; then
    if [[ ! -f "$SAMPLE_MAP" ]]; then
        echo "ERROR: sample map file not found: $SAMPLE_MAP"; exit 1
    fi
    N_MAP=$(awk 'NR>1 && NF>=2' "$SAMPLE_MAP" | wc -l | tr -d ' ')
    echo "Loaded $N_MAP barcode→sample mappings from $SAMPLE_MAP"
    echo ""
fi

mkdir -p "$OUT_DIR" logs

# ── Index the reference (once; reused across samples) ────────────────────────
REF_INDEX="${OUT_DIR}/silva_16S_index.mmi"
if [[ ! -f "$REF_INDEX" ]]; then
    echo "Indexing 16S reference with minimap2 ..."
    minimap2 -x map-ont -d "$REF_INDEX" "$REF_FASTA" -t "$THREADS"
    echo "  Index written: $REF_INDEX"
else
    echo "Using existing minimap2 index: $REF_INDEX"
fi
echo ""

# ── Process each FASTQ ────────────────────────────────────────────────────────
PROCESSED=0
SKIPPED=0

# Collect FASTQ files (plain and gzipped)
shopt -s nullglob
FASTQ_FILES=("$FASTQ_DIR"/*.fastq "$FASTQ_DIR"/*.fq \
              "$FASTQ_DIR"/*.fastq.gz "$FASTQ_DIR"/*.fq.gz)
shopt -u nullglob

if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No FASTQ files found in $FASTQ_DIR"
    exit 1
fi

echo "Found ${#FASTQ_FILES[@]} FASTQ files to process."
echo ""

for FASTQ in "${FASTQ_FILES[@]}"; do
    # ── Derive sample name ────────────────────────────────────────────────────
    BASENAME=$(basename "$FASTQ")
    # Strip extension(s): .fastq.gz, .fq.gz, .fastq, .fq
    STEM="${BASENAME%.fastq.gz}"
    STEM="${STEM%.fq.gz}"
    STEM="${STEM%.fastq}"
    STEM="${STEM%.fq}"

    # Look up sample name in map; fall back to stripping common suffixes from stem
    SAMPLE=$(lookup_sample "$STEM")
    if [[ -z "$SAMPLE" ]]; then
        SAMPLE="$STEM"
        for SUFFIX in _host_removed _host_remove _dehost _filtered _clean; do
            SAMPLE="${SAMPLE%${SUFFIX}}"
        done
        if [[ -n "$SAMPLE_MAP" ]]; then
            echo "  [WARN] $STEM not found in sample map — using stem as sample name: $SAMPLE"
        fi
    fi

    OUT_TSV="${OUT_DIR}/${SAMPLE}_16S.tsv"

    if [[ -f "$OUT_TSV" ]]; then
        echo "  [SKIP] $SAMPLE — output already exists"
        ((SKIPPED++))
        continue
    fi

    echo "  Processing: $SAMPLE"
    PAF_TMP="${OUT_DIR}/${SAMPLE}_16S_tmp.paf"

    # ── minimap2: map ONT reads to 16S reference ─────────────────────────────
    # --secondary=no   : keep only primary alignments (one best hit per read per ref)
    # -N 1             : report at most 1 secondary alignment (combined with --secondary=no
    #                    this gives one alignment per read per reference sequence)
    # -c               : output CIGAR strings (needed for accurate mismatch counts)
    # Note: minimap2 may output multiple alignments per read if the read spans
    # multiple 16S genes; these are all valid 16S hits and are kept below.
    minimap2 \
        -x map-ont \
        --secondary=no \
        -c \
        -t "$THREADS" \
        "$REF_INDEX" \
        "$FASTQ" \
        > "$PAF_TMP" 2>>"logs/count_16S_${SAMPLE}.err"

    # ── Convert PAF → BLAST-format TSV ───────────────────────────────────────
    # PAF columns (0-indexed):
    #  0  qname       1  qlen       2  qstart(0-based)  3  qend
    #  4  strand      5  tname      6  tlen
    #  7  tstart(0-based)  8  tend
    #  9  residue_matches  10 block_len  11 mapq
    #
    # Output 14 columns to match BLAST -outfmt 6 with qlen + slen:
    #  qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
    awk -v min_id="$MIN_IDENTITY" -v min_cov="$MIN_COVERAGE" '
    BEGIN { OFS="\t" }
    /^[^#]/ {
        qname  = $1
        qlen   = $2
        qstart = $3 + 1          # convert 0-based to 1-based
        qend   = $4
        tname  = $6
        tlen   = $7
        tstart = $8 + 1
        tend   = $9
        res_match = $10
        block_len = $11
        mapq   = $12

        # Skip unmapped reads (tname == "*")
        if (tname == "*") next

        # Compute identity and coverage
        if (block_len <= 0) next
        pident = (res_match / block_len) * 100
        cov    = block_len / tlen          # fraction of reference covered

        # Apply filters
        if (pident < min_id) next
        if (cov    < min_cov) next

        mismatch = block_len - res_match   # approximate (ignores insertions)
        gapopen  = 0
        evalue   = 0
        bitscore = res_match               # residue matches as quality proxy

        print qname, tname, pident, block_len, mismatch, gapopen, \
              qstart, qend, tstart, tend, evalue, bitscore, qlen, tlen
    }
    ' "$PAF_TMP" > "$OUT_TSV"

    # Clean up temporary PAF
    rm -f "$PAF_TMP"

    N_HITS=$(wc -l < "$OUT_TSV")
    N_READS=$(awk '{print $1}' "$OUT_TSV" | sort -u | wc -l)
    echo "    → $N_HITS alignments across $N_READS unique reads — saved to $(basename "$OUT_TSV")"
    ((PROCESSED++))
done

echo ""
echo "Done. Processed $PROCESSED sample(s); skipped $SKIPPED (already existed)."
echo "Output directory: $OUT_DIR"

# ── Summary table ──────────────────────────────────────────────────────────────
SUMMARY="${OUT_DIR}/16S_count_summary.tsv"
echo -e "sample\tn_alignments\tn_unique_16S_reads" > "$SUMMARY"
for TSV in "$OUT_DIR"/*_16S.tsv; do
    [[ -f "$TSV" ]] || continue
    SNAME=$(basename "$TSV" _16S.tsv)
    N_ALN=$(wc -l < "$TSV")
    N_UNQ=$(awk '{print $1}' "$TSV" | sort -u | wc -l)
    echo -e "${SNAME}\t${N_ALN}\t${N_UNQ}" >> "$SUMMARY"
done
echo ""
echo "Summary written: $SUMMARY"
cat "$SUMMARY"
