#!/bin/bash -e
#SBATCH --job-name=mag_phylo
#SBATCH --account=massey04083
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=phylo_allMAGs_%j.log
#SBATCH --error=phylo_allMAGs_%j.err

# ── Paths ─────────────────────────────────────────────────────────────────────
MAG_DIRS=(
    "/home/abradbur/00_nesi_projects/00_nesi_projects/massey04083_nobackup/poultry_processing_chapter2/MAGs/ASC_CAT"
    "/home/abradbur/00_nesi_projects/00_nesi_projects/massey04083_nobackup/poultry_processing_chapter2/MAGs/EVI_CAT"
)
OUT_DIR="/home/abradbur/00_nesi_projects/00_nesi_projects/massey04083_nobackup/poultry_processing_chapter2/phylo_output_allMAGs"
THREADS=$SLURM_CPUS_PER_TASK

mkdir -p "$OUT_DIR/prokka" "$OUT_DIR/panaroo" "$OUT_DIR/iqtree"

# ── Step 1: Prokka annotation ─────────────────────────────────────────────────
echo "=== Step 1: Prokka annotation ==="
module purge
module load prokka/1.14.6-apptainer

while IFS= read -r fasta; do
    sample=$(basename "$fasta" .fa)
    if [ -f "$OUT_DIR/prokka/$sample/$sample.gff" ]; then
        echo "  Skipping $sample (already done)"
        continue
    fi
    echo "  Annotating $sample..."
    prokka \
        --outdir "$OUT_DIR/prokka/$sample" \
        --prefix "$sample" \
        --kingdom Bacteria \
        --cpus "$THREADS" \
        --force \
        "$fasta"
done < <(find "${MAG_DIRS[@]}" -name "*.fa")

# ── Step 2: Panaroo core genome analysis ─────────────────────────────────────
echo ""
echo "=== Step 2: Panaroo core genome analysis ==="
module purge
module load panaroo/1.3.0-gimkl-2022a-Python-3.10.5 2>/dev/null || {
    echo "  Panaroo module not found, installing via pip..."
    module load Python/3.11.3-gimkl-2022a
    pip install --user panaroo 2>/dev/null || true
    export PATH="$HOME/.local/bin:$PATH"
    module load MAFFT/7.490-GCC-11.3.0
    module load CD-HIT/4.8.1-GCC-11.3.0
}

panaroo \
    -i "$OUT_DIR/prokka/"*/*.gff \
    -o "$OUT_DIR/panaroo" \
    --clean-mode strict \
    -a core \
    --aligner mafft \
    --core_threshold 0.98 \
    -t "$THREADS"

# ── Step 3: IQ-TREE3 phylogenetic tree ───────────────────────────────────────
echo ""
echo "=== Step 3: IQ-TREE3 ==="
module purge
module load IQ-TREE/3.1.1-foss-2023a

ALIGNMENT="$OUT_DIR/panaroo/core_gene_alignment.aln"

if [ ! -f "$ALIGNMENT" ]; then
    echo "ERROR: Core alignment not found at $ALIGNMENT"
    echo "Files in panaroo output:"
    ls "$OUT_DIR/panaroo/"*.aln 2>/dev/null || echo "  (no .aln files)"
    exit 1
fi

iqtree3 \
    -s "$ALIGNMENT" \
    -m GTR+G \
    -B 1000 \
    -T "$THREADS" \
    --prefix "$OUT_DIR/iqtree/all_mags_tree"

echo ""
echo "=== Done! ==="
echo "Tree file: $OUT_DIR/iqtree/all_mags_tree.treefile"
