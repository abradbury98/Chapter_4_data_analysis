# Chapter 5 Data Analysis

> Analysis scripts for Chapter 4 - Poultry processing plant.  
> Workflows cover taxonomic classification, AMR profiling, host depletion benchmarking, MAG assembly and refinement, and phylogenetic reconstruction.

---

## Repository Structure

```
Chapter_4_data_analysis/
└── scripts/
    ├── R/        Standalone R scripts — sequencing QC and rarefaction plots
    ├── python/   Python utilities — AMR processing and Kraken2 normalisation
    ├── rmd/      R Markdown notebooks — chapter figures and statistical analyses
    ├── shell/    Shell pipelines — local phylogenetics and resistance co-location
    └── slurm/    SLURM job scripts — HPC execution on NeSI
```

---

## Script Directories

| Directory | What it contains | README |
|-----------|-----------------|--------|
| [`scripts/R/`](scripts/R/) | Pore activity plots, rarefaction curves | [→ README](scripts/R/README.md) |
| [`scripts/python/`](scripts/python/) | ABRicate processing, Kraken2 BPM normalisation, MetaPhlAn visualisations | [→ README](scripts/python/README.md) |
| [`scripts/rmd/`](scripts/rmd/) | Diversity, AMR, MAG quality, base tracking, qPCR analyses | [→ README](scripts/rmd/README.md) |
| [`scripts/shell/`](scripts/shell/) | IQ-TREE2 phylogenetics, BLAST merging, resistance co-location | [→ README](scripts/shell/README.md) |
| [`scripts/slurm/`](scripts/slurm/) | Kraken2, MetaPhlAn, ABRicate, DAS Tool, GTDB-Tk on NeSI | [→ README](scripts/slurm/README.md) |

---

## Quick-start Commands

| Script type | Command |
|-------------|---------|
| R script | `Rscript scripts/R/<script>.R` |
| R Markdown | `Rscript -e "rmarkdown::render('scripts/rmd/<notebook>.Rmd')"` |
| Python | `python scripts/python/<script>.py` |
| Shell pipeline | `bash scripts/shell/<script>.sh` |
| SLURM job | `sbatch scripts/slurm/<script>.sl` |

---

## Before Running Any Script

1. **Update paths** — most scripts have hard-coded Desktop or NeSI paths at the top. Change these to match your directory layout.
2. **SLURM account** — all `.sl` files use `--account=massey04083`. Change this for a different NeSI project or cluster.
3. **Conda environments** — shell scripts activate named environments (`prokka`, `panaroo`). Ensure these exist or install the tools into the active environment.
4. **Modules on NeSI** — verify module names with `module avail <tool>` before submitting a job.
5. **Databases** — Kraken2 (nt), ABRicate (CARD, ResFinder, MEGARes, NCBI), and BLAST databases must be accessible at the paths referenced in each script.

---

## Tool Dependencies

| Tool | Version | Purpose |
|------|---------|---------|
| Kraken2 | 2.1.6 | Taxonomic classification of reads and contigs |
| MetaPhlAn | 4.1.0 | Species-level community profiling |
| ABRicate | — | AMR gene detection |
| hAMRonization | — | AMR output standardisation and deduplication |
| NanoPlot | 1.43.0 | Nanopore read quality visualisation |
| MultiQC | 1.15 | QC report aggregation |
| minimap2 | — | Long-read host depletion alignment |
| Bowtie2 | — | Short-read host depletion alignment |
| DAS Tool | 1.1.5 | Metagenomic bin refinement |
| CheckM | 1.2.3 | MAG quality (marker-gene based) |
| CheckM2 | 1.0.1 | MAG quality (ML-based) |
| GTDB-Tk | 2.4.0 | MAG taxonomic classification |
| Prokka | — | Genome annotation |
| Panaroo | — | Core genome alignment |
| IQ-TREE2 | 2.3.6 | Maximum-likelihood phylogenetic trees |
| BLAST+ | — | Sequence alignment and contamination checks |
| R / tidyverse / ggplot2 | — | Visualisation and statistics |
| Python / pandas / openpyxl | — | Data processing |

---

## Data

Raw sequencing data and processed outputs are not stored in this repository.  
Scripts reference external storage on NeSI (`/nesi/nobackup/massey04083/`) and local Desktop paths.
