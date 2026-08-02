#!/bin/bash
#SBATCH --job-name=abricate
#SBATCH --account=massey04083
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --array=0-%NJOBS%
#SBATCH --output=%LOGDIR%/abricate_%A_%a.out
#SBATCH --error=%LOGDIR%/abricate_%A_%a.err

# ── Arguments injected by submit_abricate.sl ──────────────────────────────────
INDIR="%INDIR%"
OUTDIR="%OUTDIR%"
DB="%DB%"

# ── Load abricate ─────────────────────────────────────────────────────────────
module purge
module load ABRicate/1.0.0-GCC-11.3.0-Perl-5.34.1

# ── Get file for this array task ──────────────────────────────────────────────
mapfile -t FILES < <(find "${INDIR}" -maxdepth 1 -name "*.fastq.gz" | sort)
FILE="${FILES[$SLURM_ARRAY_TASK_ID]}"

if [[ -z "$FILE" ]]; then
    echo "No file for array index ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

SAMPLE=$(basename "${FILE}" .fastq.gz)

echo "=============================="
echo "Sample  : ${SAMPLE}"
echo "File    : ${FILE}"
echo "DB      : ${DB}"
echo "Started : $(date)"
echo "=============================="

abricate \
    --datadir /home/abradbur/00_nesi_projects/00_nesi_projects/massey04083_nobackup/abricate_db \
    --db "${DB}" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --minid 80 \
    --mincov 80 \
    "${FILE}" \
    > "${OUTDIR}/${SAMPLE}_${DB}.tab"

echo "Finished : $(date)"
echo "Output   : ${OUTDIR}/${SAMPLE}_${DB}.tab"
