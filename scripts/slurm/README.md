# scripts/slurm

SLURM job scripts for running computational workflows on HPC infrastructure.

## Files

- `abricate_array.sl` – ABRicate array job submission.
- `blast_reads.sl` – BLAST-based read processing jobs.
- `checkm2_gtdbtk.sl` – CheckM2 and GTDB-Tk workflow jobs.
- `dastool.sl` – DAS Tool execution workflow.
- `host_dep_breakdown.sl` – Host depletion breakdown analyses.
- `host_depletion_benchmark.sl` – Host depletion benchmarking workflow.
- `kraken2_nt.sl` – Kraken2 classification against nt.
- `kraken2_nt_assembly_contigs.sl` – Kraken2 contig-level classification jobs.
- `metaphlan.sl` – MetaPhlAn profiling jobs.
- `nanoplot_multiqc.sl` – NanoPlot/MultiQC reporting jobs.
- `only_abricate.sl` – ABRicate-only pipeline execution.
- `phylo_nesi_IQ-Tree.sl` – IQ-TREE phylogenetics jobs on NeSI.
- `phylo_nesi_allMAGs.sl` – Phylogenetics workflow across all MAGs.
- `phylo_nesi_perspecies.sl` – Per-species phylogenetic job submission.
- `run_abricate_kraken2.sl` – Combined ABRicate + Kraken2 job script.
- `run_checkm2_gtdbtk.sl` – Wrapper/runner for CheckM2 + GTDB-Tk jobs.
- `submit_abricate.sl` – ABRicate submission helper.
- `submit_abricate_all_db.sl` – ABRicate submission across all databases.

## Typical usage

Submit jobs with:

```bash
sbatch scripts/slurm/<script_name>.sl
```

## HPC notes

- Update account/partition/time/memory directives as needed for your cluster.
- Verify required modules/conda environments and reference databases are available.
