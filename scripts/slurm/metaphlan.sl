#!/bin/bash -e
#SBATCH --account       massey04083
#SBATCH --job-name      metaphlan
#SBATCH --time          08:00:00
#SBATCH --cpus-per-task 16
#SBATCH --mem           120G
#SBATCH --output        metaphlan_%j.out
#SBATCH --error         metaphlan_%j.err

module purge
module load MetaPhlAn/4.1.0-gimkl-2022a-Python-3.10.5

INDIR=/home/abradbur/00_nesi_projects/00_nesi_projects/massey04083_nobackup/poultry_processing_chapter2/ASC_BFBB/host_removed_reads
OUTDIR=/home/abradbur/00_nesi_projects/00_nesi_projects/massey04083_nobackup/poultry_processing_chapter2/ASC_BFBB/metaphlan_output
DBDIR=/home/abradbur/00_nesi_projects/00_nesi_projects/massey04083_nobackup/metaphlan_db

mkdir -p $OUTDIR

for READS in $INDIR/*.fastq.gz; do
    SAMPLE=$(basename $READS | grep -o 'barcode[0-9]*')

    if [[ -s $OUTDIR/${SAMPLE}_profile.txt ]]; then
        echo "Skipping $SAMPLE (already done)"
        continue
    fi

    echo "Processing $SAMPLE..."

    metaphlan $READS \
        --input_type fastq \
        --bowtie2db $DBDIR \
        --nproc $SLURM_CPUS_PER_TASK \
        --read_min_len 300 \
        --min_alignment_len 300 \
        -o $OUTDIR/${SAMPLE}_profile.txt \
        --bowtie2out $OUTDIR/${SAMPLE}.bowtie2.bz2

done

merge_metaphlan_tables.py $OUTDIR/*_profile.txt > $OUTDIR/merged_abundance_table.txt

echo "Done."