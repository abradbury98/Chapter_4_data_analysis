#!/usr/bin/env bash
# Phylogenetic tree pipeline: Prokka -> Panaroo -> IQ-TREE2
# Input: Campylobacter coli MAGs in FASTA format
# Output: Maximum-likelihood phylogenetic tree
#
# Prerequisites (run once before this script):
#   CONDA_SUBDIR=osx-64 mamba create -n prokka_env -c bioconda -c conda-forge prokka -y
#   mamba create -n panaroo_env -c conda-forge python=3.9 mafft cd-hit -y
#   conda run -n panaroo_env pip install panaroo
#   curl -L https://github.com/iqtree/iqtree2/releases/download/v2.3.6/iqtree-2.3.6-MacOSX-arm.tar.gz \
#       | tar -xz -C ~/

set -euo pipefail

MAG_DIR="/Users/alicebradbury/Desktop/Campylobacter MAGs"
OUT_DIR="/Users/alicebradbury/Desktop/phylo_output"
IQTREE2="$HOME/iqtree-2.3.6-MacOSX-arm/bin/iqtree2"
THREADS=4

mkdir -p "$OUT_DIR/prokka" "$OUT_DIR/panaroo" "$OUT_DIR/iqtree"

# ── Step 1: Annotate each MAG with Prokka ────────────────────────────────────
echo "=== Step 1: Prokka annotation ==="
for fasta in "$MAG_DIR"/*.fasta; do
    sample=$(basename "$fasta" .fasta)
    if [ -f "$OUT_DIR/prokka/$sample/$sample.gff" ]; then
        echo "  Skipping $sample (already done)"
        continue
    fi
    echo "  Annotating $sample..."
    conda run -n prokka_env prokka \
        --outdir "$OUT_DIR/prokka/$sample" \
        --prefix "$sample" \
        --genus Campylobacter \
        --species coli \
        --kingdom Bacteria \
        --cpus "$THREADS" \
        --force \
        "$fasta"
done

# ── Step 2: Core genome alignment with Panaroo ───────────────────────────────
echo ""
echo "=== Step 2: Panaroo core genome analysis ==="

conda run -n panaroo_env panaroo \
    -i "$OUT_DIR/prokka/"*/*.gff \
    -o "$OUT_DIR/panaroo" \
    --clean-mode strict \
    -a core \
    --aligner mafft \
    --core_threshold 0.98 \
    -t "$THREADS"

# ── Step 3: Build phylogenetic tree with IQ-TREE2 ────────────────────────────
echo ""
echo "=== Step 3: IQ-TREE2 phylogenetic tree ==="
ALIGNMENT="$OUT_DIR/panaroo/core_gene_alignment.aln"

if [ ! -f "$ALIGNMENT" ]; then
    echo "ERROR: Core alignment not found at $ALIGNMENT"
    echo "Files present in panaroo output:"
    ls "$OUT_DIR/panaroo/"*.aln 2>/dev/null || echo "  (no .aln files found)"
    exit 1
fi

"$IQTREE2" \
    -s "$ALIGNMENT" \
    -m GTR+G \
    -B 1000 \
    -T "$THREADS" \
    --prefix "$OUT_DIR/iqtree/campylobacter_coli_tree"

echo ""
echo "=== Done! ==="
echo "Tree file: $OUT_DIR/iqtree/campylobacter_coli_tree.treefile"
echo "Visualise with FigTree or upload to iTOL: https://itol.embl.de"
