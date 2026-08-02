#!/bin/bash
# submit_abricate.sl — generic abricate submission script
# Usage: bash submit_abricate.sl <input_dir> [database]
#
# Examples:
#   bash submit_abricate.sl /path/to/host_removed_reads ecoli_vf
#   bash submit_abricate.sl /path/to/host_removed_reads vfdb
#   bash submit_abricate.sl /path/to/host_removed_reads ncbi
#
# Output is written to a sibling folder of the input dir:
#   <parent_of_input_dir>/abricate_<database>/
#
# Both scripts (submit_abricate.sl and abricate_array.sl) must be in:
#   /home/abradbur/00_nesi_projects/00_nesi_projects/massey04083_nobackup/scripts/

SCRIPTS="/home/abradbur/00_nesi_projects/00_nesi_projects/massey04083_nobackup/retail_poultry_chapter3/China/CH01-06_CHINA_RVS"
TEMPLATE="${SCRIPTS}/abricate_array.sl"

# ── Arguments ─────────────────────────────────────────────────────────────────
INDIR="${1}"
DB="${2:-ecoli_vf}"   # default to ecoli_vf if not specified

if [[ -z "$INDIR" ]]; then
    echo "Usage: bash submit_abricate.sl <input_dir> [database]"
    echo "  database defaults to ecoli_vf if not specified"
    exit 1
fi

if [[ ! -d "$INDIR" ]]; then
    echo "ERROR: Input directory not found: ${INDIR}"
    exit 1
fi

if [[ ! -f "$TEMPLATE" ]]; then
    echo "ERROR: Template script not found: ${TEMPLATE}"
    echo "Make sure abricate_array.sl is in ${SCRIPTS}"
    exit 1
fi

# ── Derived paths ─────────────────────────────────────────────────────────────
PARENT=$(dirname "${INDIR}")
OUTDIR="${PARENT}/abricate_${DB}"
LOGDIR="${PARENT}/logs_abricate_${DB}"

mkdir -p "${OUTDIR}" "${LOGDIR}"

# ── Count files ───────────────────────────────────────────────────────────────
N=$(find "${INDIR}" -maxdepth 1 -name "*.fastq.gz" | wc -l)

if [[ "$N" -eq 0 ]]; then
    echo "ERROR: No .fastq.gz files found in ${INDIR}"
    exit 1
fi

echo "Input dir : ${INDIR}"
echo "Database  : ${DB}"
echo "Output    : ${OUTDIR}"
echo "Files     : ${N} fastq.gz"

# ── Patch template ────────────────────────────────────────────────────────────
LAST=$(( N - 1 ))
SUBMIT="${PARENT}/abricate_${DB}_submit.sl"

sed \
    -e "s|%NJOBS%|${LAST}|g" \
    -e "s|%INDIR%|${INDIR}|g" \
    -e "s|%OUTDIR%|${OUTDIR}|g" \
    -e "s|%LOGDIR%|${LOGDIR}|g" \
    -e "s|%DB%|${DB}|g" \
    "${TEMPLATE}" > "${SUBMIT}"

# ── Submit array ──────────────────────────────────────────────────────────────
JOB_ID=$(sbatch --parsable "${SUBMIT}")
echo "Submitted array job ${JOB_ID} (${N} tasks)"

# ── Submit summary after array completes ──────────────────────────────────────
sbatch --dependency=afterok:${JOB_ID} \
       --job-name=abricate_summary \
       --account=massey04083 \
       --time=00:10:00 \
       --cpus-per-task=1 \
       --mem=2G \
       --output="${LOGDIR}/summary.out" \
       --wrap="module purge; module load ABRicate/1.0.0-GCC-11.3.0-Perl-5.34.1; \
               abricate --summary ${OUTDIR}/*_${DB}.tab > ${OUTDIR}/summary_${DB}.tab; \
               echo 'Done: ${OUTDIR}/summary_${DB}.tab'"

echo ""
echo "Results will be in : ${OUTDIR}/"
echo "Logs in            : ${LOGDIR}/"
echo "Monitor with       : squeue -u abradbur"
