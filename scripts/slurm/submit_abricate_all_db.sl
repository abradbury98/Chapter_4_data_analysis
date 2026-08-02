#!/bin/bash
# submit_abricate_all_db.sl — run abricate across multiple databases
# Usage: bash submit_abricate_all_db.sl <input_dir>
#
# Example:
#   bash submit_abricate_all_db.sl \
#     /home/abradbur/.../CH01-06_CHINA_RVS/host_removed_reads
#
# Produces one output folder per database, e.g.:
#   <parent>/abricate_ecoli_vf/
#   <parent>/abricate_ecoh/
#   <parent>/abricate_vfdb/
#   <parent>/abricate_megares/
#   <parent>/abricate_plasmidfinder/

SCRIPTS="/home/abradbur/00_nesi_projects/00_nesi_projects/massey04083_nobackup/retail_poultry_chapter3/China/CH01-06_CHINA_RVS"
SUBMIT="${SCRIPTS}/submit_abricate.sl"

# ── Databases to run ──────────────────────────────────────────────────────────
DATABASES=(
    ecoli_vf
    ecoh
    vfdb
    megares
    plasmidfinder
)

# ── Arguments ─────────────────────────────────────────────────────────────────
INDIR="${1}"

if [[ -z "$INDIR" ]]; then
    echo "Usage: bash submit_abricate_all_db.sl <input_dir>"
    exit 1
fi

if [[ ! -d "$INDIR" ]]; then
    echo "ERROR: Input directory not found: ${INDIR}"
    exit 1
fi

if [[ ! -f "$SUBMIT" ]]; then
    echo "ERROR: submit_abricate.sl not found at ${SUBMIT}"
    exit 1
fi

# ── Submit one array job per database ────────────────────────────────────────
echo "Input dir : ${INDIR}"
echo "Databases : ${DATABASES[*]}"
echo ""

for DB in "${DATABASES[@]}"; do
    echo "--- Submitting: ${DB} ---"
    bash "${SUBMIT}" "${INDIR}" "${DB}"
    echo ""
done

echo "All databases submitted."
echo "Monitor with: squeue -u abradbur"
