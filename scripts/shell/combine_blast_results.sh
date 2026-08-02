#!/bin/bash
# Run after all array tasks complete to merge results into one table.
# Usage: bash combine_blast_results.sh

OUTDIR="blast_results"
FINAL="${OUTDIR}/all_samples_blast_top3.tsv"

echo -e "sample\tread_id\tsubject_id\tpct_identity\talign_length\tevalue\tbitscore\tsubject_title\tscientific_name\thit_rank" \
    > "${FINAL}"

cat "${OUTDIR}"/*_blast_top3.tsv >> "${FINAL}"

echo "Combined table: ${FINAL}"
echo "Total hits: $(tail -n +2 "${FINAL}" | wc -l)"
