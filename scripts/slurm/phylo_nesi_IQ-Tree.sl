#!/bin/bash -e
#SBATCH --job-name=campy_phylo
#SBATCH --account=massey04083
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=phylo_IQTree_%j.log
#SBATCH --error=phylo_IQTree_%j.err

# ── Paths ─────────────────────────────────────────────────────────────────────
MAG_DIRS=(
    "/home/abradbur/00_nesi_projects/00_nesi_projects/massey04083_nobackup/poultry_processing_chapter2/MAGs/ASC_CAT"
    "/home/abradbur/00_nesi_projects/00_nesi_projects/massey04083_nobackup/poultry_processing_chapter2/MAGs/EVI_CAT"
)
OUT_DIR="/home/abradbur/00_nesi_projects/00_nesi_projects/massey04083_nobackup/poultry_processing_chapter2/phylo_output"
THREADS=$SLURM_CPUS_PER_TASK


# ── Step 3: IQ-TREE2 phylogenetic tree ───────────────────────────────────────
echo ""
echo "=== Step 3: IQ-TREE2 ==="
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
    --prefix "$OUT_DIR/iqtree/campylobacter_coli_tree"

echo ""
echo "=== Done! ==="
echo "Tree file: $OUT_DIR/iqtree/campylobacter_coli_tree.treefile"
