# scripts/shell

Shell pipelines and helpers that orchestrate AMR profiling, BLAST processing, phylogeny, and resistance co-location analysis.

## Files

- `combine_blast_results.sh`  
  Combines BLAST result files into a single dataset.

- `phylo_pipeline.sh`  
  End-to-end or multi-step phylogenetic workflow orchestration.

- `rRNA_SNP_check_resistance.sh`  
  Checks rRNA SNP-associated resistance patterns.

- `resistance_colocation_pipeline_ASC.sh`  
  Pipeline for resistance co-location analysis (ASC workflow variant).

- `resistance_colocation_pipeline_ASC 2.sh`  
  Alternate/copy variant of the ASC resistance co-location pipeline.

- `run_abricate_MAGs.sh`  
  Runs ABRicate against MAG inputs.

- `run_abricate_chunks.sh`  
  Runs ABRicate in chunks/batches for scalability.

- `run_checkm2_gtdbtk.sl`  
  Script for CheckM2/GTDB-Tk execution (SLURM-style extension despite location).

## Typical usage

```bash
bash scripts/shell/<script_name>.sh
```

Some scripts may require editing paths, database locations, and conda/module setup before execution.
