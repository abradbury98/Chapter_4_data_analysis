#!/bin/bash -e
#SBATCH --job-name=mag_phylo_perspecies
#SBATCH --account=massey04083
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=phylo_perspecies_%j.log
#SBATCH --error=phylo_perspecies_%j.err

# ── Paths ─────────────────────────────────────────────────────────────────────
PROKKA_DIR="/home/abradbur/00_nesi_projects/00_nesi_projects/massey04083_nobackup/poultry_processing_chapter2/phylo_output_allMAGs/prokka"
OUT_DIR="/home/abradbur/00_nesi_projects/00_nesi_projects/massey04083_nobackup/poultry_processing_chapter2/phylo_output_perspecies"
THREADS=$SLURM_CPUS_PER_TASK
MIN_MAGS=3   # minimum MAGs required to build a tree for a species

mkdir -p "$OUT_DIR"

# ── Step 1: Discover species with enough MAGs ─────────────────────────────────
echo "=== Discovering species ==="

declare -A SPECIES_COUNT

for dir in "$PROKKA_DIR"/*/; do
    dirname=$(basename "$dir")
    # Extract genus_species (last two underscore-delimited fields)
    species=$(echo "$dirname" | rev | cut -d'_' -f1,2 | rev)
    SPECIES_COUNT[$species]=$(( ${SPECIES_COUNT[$species]:-0} + 1 ))
done

echo "Species found:"
for species in "${!SPECIES_COUNT[@]}"; do
    echo "  ${SPECIES_COUNT[$species]}x  $species"
done

# ── Load modules ──────────────────────────────────────────────────────────────
module purge
module load panaroo/1.3.0-gimkl-2022a-Python-3.10.5 2>/dev/null || {
    echo "  Panaroo module not found, installing via pip..."
    module load Python/3.11.3-gimkl-2022a
    pip install --user panaroo 2>/dev/null || true
    export PATH="$HOME/.local/bin:$PATH"
    module load MAFFT/7.490-GCC-11.3.0
    module load CD-HIT/4.8.1-GCC-11.3.0
}

# ── Step 2 & 3: Panaroo + IQ-TREE per species ────────────────────────────────
for species in "${!SPECIES_COUNT[@]}"; do

    count=${SPECIES_COUNT[$species]}
    if [ "$count" -lt "$MIN_MAGS" ]; then
        echo ""
        echo "Skipping $species ($count MAGs — below minimum of $MIN_MAGS)"
        continue
    fi

    echo ""
    echo "======================================================"
    echo "Processing: $species ($count MAGs)"
    echo "======================================================"

    SPECIES_OUT="$OUT_DIR/$species"
    mkdir -p "$SPECIES_OUT/panaroo" "$SPECIES_OUT/iqtree"

    # Collect GFF files for this species
    GFF_FILES=()
    for dir in "$PROKKA_DIR"/*_${species}/; do
        sample=$(basename "$dir")
        gff="$dir/${sample}.gff"
        if [ -f "$gff" ]; then
            GFF_FILES+=("$gff")
        fi
    done

    echo "  Found ${#GFF_FILES[@]} GFF files"

    if [ "${#GFF_FILES[@]}" -lt "$MIN_MAGS" ]; then
        echo "  Skipping — not enough GFF files found"
        continue
    fi

    # ── Panaroo ───────────────────────────────────────────────────────────────
    if [ -f "$SPECIES_OUT/panaroo/core_gene_alignment.aln" ]; then
        echo "  Panaroo already done — skipping"
    else
        echo "  Running Panaroo..."
        panaroo \
            -i "${GFF_FILES[@]}" \
            -o "$SPECIES_OUT/panaroo" \
            --clean-mode strict \
            -a core \
            --aligner mafft \
            --core_threshold 0.98 \
            -t "$THREADS"
    fi

    # ── IQ-TREE ───────────────────────────────────────────────────────────────
    ALIGNMENT="$SPECIES_OUT/panaroo/core_gene_alignment.aln"

    if [ ! -f "$ALIGNMENT" ] || [ ! -s "$ALIGNMENT" ]; then
        echo "  WARNING: No core alignment produced for $species — skipping IQ-TREE"
        echo "  (core genome may be too small at 98% threshold)"
        continue
    fi

    if [ -f "$SPECIES_OUT/iqtree/${species}_tree.treefile" ]; then
        echo "  IQ-TREE already done — skipping"
        continue
    fi

    echo "  Running IQ-TREE3..."
    module purge
    module load IQ-TREE/3.1.1-foss-2023a

    iqtree3 \
        -s "$ALIGNMENT" \
        -m GTR+G \
        -B 1000 \
        -T "$THREADS" \
        --prefix "$SPECIES_OUT/iqtree/${species}_tree"

    # Reload panaroo module for next species
    module purge
    module load panaroo/1.3.0-gimkl-2022a-Python-3.10.5 2>/dev/null || {
        module load Python/3.11.3-gimkl-2022a
        export PATH="$HOME/.local/bin:$PATH"
        module load MAFFT/7.490-GCC-11.3.0
        module load CD-HIT/4.8.1-GCC-11.3.0
    }

done

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "======================================================"
echo "=== All done! Tree files produced: ==="
find "$OUT_DIR" -name "*.treefile" | sort
echo "======================================================"
