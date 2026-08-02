# scripts/slurm

SLURM job submission scripts for HPC execution on NeSI (project `massey04083`).  
Submit with `sbatch scripts/slurm/<script>.sl`. Update account, paths, and resource directives before submitting.

---

## Contents

| Script | Category | CPUs | Memory | Time |
|--------|----------|------|--------|------|
| [`kraken2_nt.sl`](#kraken2_ntsl) | Taxonomic classification | 24 | 240 GB | 8 h |
| [`kraken2_nt_assembly_contigs.sl`](#kraken2_nt_assembly_contigssl) | Taxonomic classification | — | — | — |
| [`metaphlan.sl`](#metaphlansl) | Taxonomic classification | 16 | 120 GB | 8 h |
| [`nanoplot_multiqc.sl`](#nanoplot_multiqcsl) | Quality control | 16 | 32 GB | 12 h |
| [`host_depletion_benchmark.sl`](#host_depletion_benchmarksl) | Quality control | — | — | — |
| [`host_dep_breakdown.sl`](#host_dep_breakdownsl) | Quality control | — | — | — |
| [`abricate_array.sl`](#abricate_arraysl) | AMR profiling | — | — | — |
| [`only_abricate.sl`](#only_abricatesl) | AMR profiling | — | — | — |
| [`run_abricate_kraken2.sl`](#run_abricate_kraken2sl) | AMR profiling | — | — | — |
| [`submit_abricate.sl`](#submit_abricatesl) | AMR profiling | — | — | — |
| [`submit_abricate_all_db.sl`](#submit_abricate_all_dbsl) | AMR profiling | — | — | — |
| [`dastool.sl`](#dastoolsl) | MAG refinement | 16 | 180 GB | 12 h |
| [`checkm2_gtdbtk.sl`](#checkm2_gtdbtksl) | MAG refinement | — | — | — |
| [`run_checkm2_gtdbtk.sl`](#run_checkm2_gtdbtksl) | MAG refinement | — | — | — |
| [`phylo_nesi_IQ-Tree.sl`](#phylo_nesi_iq-treesl) | Phylogenetics | — | — | — |
| [`phylo_nesi_allMAGs.sl`](#phylo_nesi_allmagssal) | Phylogenetics | — | — | — |
| [`phylo_nesi_perspecies.sl`](#phylo_nesi_perspeciessl) | Phylogenetics | — | — | — |
| [`blast_reads.sl`](#blast_readssl) | Other | — | — | — |

---

## Before Submitting Any Job

1. **Account** — all scripts use `--account=massey04083`. Change this for a different NeSI project or cluster.
2. **Input/output paths** — update the hardcoded directory variables at the top of each script.
3. **Modules** — verify module availability on your cluster:
   ```bash
   module avail <tool>
   ```
4. **Resources** — adjust `--time`, `--mem`, and `--cpus-per-task` to match queue policies and your dataset size.
5. **Databases** — confirm the following are accessible before submitting:
   - Kraken2 nt database
   - ABRicate databases (CARD, ResFinder, MEGARes, NCBI)
   - MetaPhlAn database
   - GTDB-Tk reference database
   - Reference genomes for host depletion (DCS, human GRCh38, chicken broiler)

---

## Taxonomic Classification

### `kraken2_nt.sl`

Classifies host-depleted reads against the NCBI nucleotide (nt) database using Kraken2.

#### SLURM Directives

| Directive | Value |
|-----------|-------|
| `--account` | massey04083 |
| `--cpus-per-task` | 24 |
| `--mem` | 240G |
| `--time` | 08:00:00 |

#### Module

```
module load Kraken2/2.1.6-GCC-11.3.0
```

#### Inputs

Host-removed FASTQ files (`.fastq.gz`) from the poultry processing study directory. Update `$INPUT_DIR` at the top of the script.

#### Outputs

| File | Description |
|------|-------------|
| `<sample>.out` | Per-read Kraken2 classification assignments |
| `<sample>.report` | Taxonomy summary report |
| `logs/<sample>.err` | Per-sample error log |

#### Settings

Uses a 0.1 confidence threshold (`--confidence 0.1`). Iterates over all compressed FASTQ files; continues processing if individual samples fail.

#### Running

```bash
sbatch scripts/slurm/kraken2_nt.sl
```

---

### `kraken2_nt_assembly_contigs.sl`

Runs Kraken2 classification on assembled contigs rather than raw reads.

#### Inputs

Assembled contig FASTA files — output from a Flye or equivalent long-read assembler.

#### Outputs

Contig-level Kraken2 `.out` and `.report` files per sample.

#### Running

```bash
sbatch scripts/slurm/kraken2_nt_assembly_contigs.sl
```

---

### `metaphlan.sl`

Runs MetaPhlAn4 for species-level community profiling and merges all sample profiles into a single abundance table.

#### SLURM Directives

| Directive | Value |
|-----------|-------|
| `--account` | massey04083 |
| `--cpus-per-task` | 16 |
| `--mem` | 120G |
| `--time` | 08:00:00 |

#### Module

```
module load MetaPhlAn/4.1.0-gimkl-2022a-Python-3.10.5
```

#### Inputs

| Variable | Description |
|----------|-------------|
| `$INPUT_DIR` | Directory containing host-removed `.fastq.gz` files |
| `$DB_DIR` | MetaPhlAn database directory (`metaphlan_db/`) |

#### Outputs

| File | Description |
|------|-------------|
| `metaphlan_output/<sample>_profile.txt` | Per-sample taxonomic profile (clade-relative abundances) |
| `metaphlan_output/<sample>.bowtie2.bz2` | Bowtie2 alignment file |
| `metaphlan_output/merged_abundance_table.txt` | All samples combined — input for MetaPhlAn visualisations |

#### Settings

Minimum read length: 300 bp. Minimum alignment length: 300 bp. Skips already-completed samples.

#### Running

```bash
sbatch scripts/slurm/metaphlan.sl
```

---

## Quality Control

### `nanoplot_multiqc.sl`

Generates per-sample Nanopore QC reports with NanoPlot and aggregates them into a MultiQC dashboard.

#### SLURM Directives

| Directive | Value |
|-----------|-------|
| `--account` | massey04083 |
| `--cpus-per-task` | 16 |
| `--mem` | 32G |
| `--time` | 12:00:00 |

#### Modules

```
module load NanoPlot/1.43.0-foss-2023a-Python-3.11.6
module load MultiQC/1.15-gimkl-2022a-Python-3.10.5
```

#### Inputs

`*.fastq.gz` files in a user-specified working directory (`$WD`).

#### Outputs

| Location | Description |
|----------|-------------|
| `$WD/Nanoplot/<sample>/` | Per-sample NanoPlot HTML and stats files |
| `$WD/Nanoplot/MultiQC/multiqc_<sample>.html` | Per-sample MultiQC report |
| `$WD/Nanoplot/MultiQC/multiqc_all_samples.html` | Aggregated MultiQC dashboard |

Parallelises NanoPlot across up to 6 concurrent jobs. Skips already-completed samples. Generates an HTML index linking all per-sample reports.

#### Running

```bash
sbatch scripts/slurm/nanoplot_multiqc.sl
```

---

### `host_depletion_benchmark.sl`

Benchmarks eight parameter configurations for host sequence removal (DCS, human GRCh38, chicken broiler) across minimap2 and Bowtie2 aligners.

#### Inputs

| Input | Description |
|-------|-------------|
| FASTQ files | Raw Nanopore reads in `$INPUT_DIR` |
| DCS reference | FASTA for DNA CS spike-in |
| Human reference | GRCh38 FASTA |
| Chicken reference | Broiler chicken genome FASTA |

Reads are subsampled to 500,000 with a fixed seed for reproducibility before benchmarking.

#### Outputs

| File | Description |
|------|-------------|
| `benchmark_summary.tsv` | Reads retained and contamination % per configuration |
| Filtered FASTQ files | Host-depleted reads for each parameter set |
| Kraken2 reports | Residual contamination assessment |
| BLAST results | Per-sample-configuration alignment outputs |

#### Running

```bash
sbatch scripts/slurm/host_depletion_benchmark.sl
```

---

### `host_dep_breakdown.sl`

Reports the number of reads and bases removed at each host depletion step for per-sample pipeline tracking.

#### Inputs

Host depletion output logs or FASTQ files from each processing stage (DCS → human → chicken).

#### Outputs

Per-sample breakdown table: reads removed at each depletion step.

#### Running

```bash
sbatch scripts/slurm/host_dep_breakdown.sl
```

---

## AMR Profiling

### `abricate_array.sl`

SLURM array job that runs ABRicate in parallel across multiple input files.

#### Inputs

FASTA files (reads, contigs, or MAGs) and an ABRicate database (e.g., `card`, `resfinder`, `megares`).

#### Outputs

Per-file ABRicate TSV result files.

#### Running

```bash
sbatch scripts/slurm/abricate_array.sl
```

---

### `only_abricate.sl`

Single-sample ABRicate run — use for testing or when only a small number of files need processing.

#### Inputs / Outputs

Same format as `abricate_array.sl`.

#### Running

```bash
sbatch scripts/slurm/only_abricate.sl
```

---

### `run_abricate_kraken2.sl`

Integrated pipeline that runs Kraken2 taxonomic classification followed immediately by ABRicate AMR screening on the same input files.

#### Inputs

FASTQ or FASTA files.

#### Outputs

Kraken2 `.out`/`.report` files and ABRicate TSV files per sample.

#### Running

```bash
sbatch scripts/slurm/run_abricate_kraken2.sl
```

---

### `submit_abricate.sl`

Wrapper that submits ABRicate jobs for a predefined set of samples. Edit the sample list inside the script before submitting.

#### Running

```bash
sbatch scripts/slurm/submit_abricate.sl
```

---

### `submit_abricate_all_db.sl`

Submits ABRicate jobs across all configured databases (CARD, ResFinder, MEGARes, NCBI) for comprehensive AMR screening in one submission.

#### Running

```bash
sbatch scripts/slurm/submit_abricate_all_db.sl
```

---

## MAG Refinement

### `dastool.sl`

Consolidates bins from three binning tools using DAS Tool, then assesses quality and assigns taxonomy.

#### SLURM Directives

| Directive | Value |
|-----------|-------|
| `--account` | massey04083 |
| `--cpus-per-task` | 16 |
| `--mem` | 180G |
| `--time` | 12:00:00 |

#### Modules

```
module load DAS_Tool/1.1.5-gimkl-2022a-R-4.2.1
module load CheckM/1.2.3-gimkl-2022a-Python-3.10.5
module load CheckM2/1.0.1-gimkl-2022a-Python-3.10.5
module load GTDB-Tk/2.4.0-gimkl-2022a-Python-3.10.5
```

#### Inputs

| Directory | Description |
|-----------|-------------|
| `bins_metabat2_*/` | MetaBat2 bins (`.fa`) |
| `bins_maxbin2_*/` | MaxBin2 bins (`.fasta`) |
| `bins_concoct_*/fasta_bins/` | CONCOCT bins (`.fa`) |
| `flye_assembly/` | Flye assembly contigs |

#### Workflow

```
MetaBat2 + MaxBin2 + CONCOCT bins
    ↓  DAS Tool (bin consolidation)
    Refined bins
    ↓  CheckM (marker-gene quality)
    ↓  CheckM2 (ML-based quality)
    ↓  GTDB-Tk (taxonomic classification)
    Summary report
```

#### Outputs

| Output | Description |
|--------|-------------|
| `dastool_output/` | Refined bin FASTA files |
| `checkm_results.tsv` | Completeness and contamination (marker-gene) |
| `checkm2_results.tsv` | Completeness and contamination (ML-based) |
| `gtdbtk_summary.tsv` | Taxonomic classifications |
| `bin_summary_report.txt` | Combined statistics per bin |

#### Running

```bash
sbatch scripts/slurm/dastool.sl
```

---

### `checkm2_gtdbtk.sl`

Runs CheckM2 quality assessment and GTDB-Tk taxonomic classification on MAG FASTA files, without the DAS Tool refinement step.

#### Inputs

MAG FASTA files in a specified directory.

#### Outputs

CheckM2 quality report and GTDB-Tk taxonomy assignments.

#### Running

```bash
sbatch scripts/slurm/checkm2_gtdbtk.sl
```

---

### `run_checkm2_gtdbtk.sl`

Alternative / updated version of `checkm2_gtdbtk.sl`. Compare both before submitting to confirm which is current.

#### Running

```bash
sbatch scripts/slurm/run_checkm2_gtdbtk.sl
```

---

## Phylogenetics

### `phylo_nesi_IQ-Tree.sl`

Constructs a maximum-likelihood phylogenetic tree using IQ-TREE2 from a pre-computed core genome alignment.

#### Inputs

Core genome FASTA alignment — produced by Panaroo (from `phylo_nesi_allMAGs.sl` or the local `phylo_pipeline.sh`).

#### Outputs

| File | Description |
|------|-------------|
| `<prefix>.treefile` | Maximum-likelihood phylogenetic tree |
| `<prefix>.log` | IQ-TREE2 run log |
| Bootstrap files | Support values for each node |

#### Running

```bash
sbatch scripts/slurm/phylo_nesi_IQ-Tree.sl
```

---

### `phylo_nesi_allMAGs.sl`

Runs the full phylogenetic pipeline (annotation → core genome → tree) across all MAGs on NeSI.

#### Inputs

MAG FASTA files in a specified directory.

#### Outputs

GFF annotation files (Prokka), core genome FASTA alignment (Panaroo), and a maximum-likelihood phylogenetic tree (IQ-TREE2).

#### Running

```bash
sbatch scripts/slurm/phylo_nesi_allMAGs.sl
```

---

### `phylo_nesi_perspecies.sl`

Runs the phylogenetic pipeline independently for each species group identified in the MAG set.

#### Inputs

MAG FASTA files grouped by species (directory structure or a sample sheet defining groupings).

#### Outputs

Per-species core genome alignments and phylogenetic trees.

#### Running

```bash
sbatch scripts/slurm/phylo_nesi_perspecies.sl
```

---

## Other

### `blast_reads.sl`

BLASTs reads or contigs against a reference nucleotide or protein database for sequence identity verification or contamination checks.

#### Inputs

| Input | Description |
|-------|-------------|
| FASTA / FASTQ files | Reads or contigs to query |
| BLAST database | Nucleotide (nt) or protein (nr) database |

#### Outputs

Tabular BLAST output files (format 6) per sample.

#### Running

```bash
sbatch scripts/slurm/blast_reads.sl
```
