# Chapter_4_data_analysis

This repository contains analysis scripts grouped by language and execution environment under `scripts/`.

## Repository structure

- `scripts/R` – R and R Markdown scripts for plotting and rarefaction analyses.
- `scripts/python` – Python scripts for processing and normalizing AMR and taxonomy outputs.
- `scripts/rmd` – R Markdown notebooks for chapter-level analysis and figure generation.
- `scripts/shell` – Shell pipelines and helper scripts for AMR, BLAST, phylogeny, and resistance workflows.
- `scripts/slurm` – SLURM job submission scripts for HPC execution of key workflows.

## Script documentation

Each scripts subfolder has its own README with file-level descriptions:

- [`scripts/R/README.md`](scripts/R/README.md)
- [`scripts/python/README.md`](scripts/python/README.md)
- [`scripts/rmd/README.md`](scripts/rmd/README.md)
- [`scripts/shell/README.md`](scripts/shell/README.md)
- [`scripts/slurm/README.md`](scripts/slurm/README.md)

## Notes

- Most workflows appear to be tailored to metagenomics, AMR profiling, host depletion benchmarking, and phylogenetic analysis.
- Many scripts expect local/HPC-specific paths and software environments; review each script before execution.
