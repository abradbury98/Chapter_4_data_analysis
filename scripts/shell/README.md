# scripts/shell

Shell pipelines for local execution (Mac/Linux), orchestrating phylogenetic tree construction, BLAST result processing, AMR profiling, and resistance co-location analysis.

---

## Contents

| Script | Purpose |
|--------|---------|
| [`phylo_pipeline.sh`](#phylo_pipelinesh) | End-to-end *C. coli* MAG phylogenetics (Prokka → Panaroo → IQ-TREE2) |
| [`combine_blast_results.sh`](#combine_blast_resultssh) | Merge per-sample BLAST outputs into one file |
| [`rRNA_SNP_check_resistance.sh`](#rrna_snp_check_resistancesh) | Identify resistance-conferring SNPs in rRNA genes |
| [`resistance_colocation_pipeline_ASC.sh`](#resistance_colocation_pipeline_ascsh) | Resistance gene co-location analysis (ASC samples) |
| [`resistance_colocation_pipeline_ASC 2.sh`](#resistance_colocation_pipeline_asc-2sh) | Alternative version of the co-location pipeline |
| [`run_abricate_MAGs.sh`](#run_abricate_magssh) | Run ABRicate on all MAG FASTA files |
| [`run_abricate_chunks.sh`](#run_abricate_chunkssh) | Run ABRicate in batches on large file sets |
| [`run_checkm2_gtdbtk.sl`](#run_checkm2_gtdbtkcsl) | Local wrapper for CheckM2 and GTDB-Tk |

---

## Before Running

1. **Update hard-coded paths** — directory paths, binary locations, and database paths are set inside each script. Edit these before running.
2. **Conda environments** — scripts activate named environments. Ensure they exist:
   ```bash
   conda env list
   ```
3. **Tools on PATH** — verify tools are accessible:
   ```bash
   which prokka panaroo iqtree2 abricate blastn
   ```
4. **ABRicate databases** — confirm databases are downloaded:
   ```bash
   abricate --list
   ```

---

## Scripts

### `phylo_pipeline.sh`

Automates construction of a maximum-likelihood phylogenetic tree from *Campylobacter coli* MAGs.

#### Workflow

```
MAG FASTA files
    ↓  Prokka (annotation)
    GFF files
    ↓  Panaroo strict mode (core genome alignment, 98% presence threshold)
    core_gene_alignment.aln
    ↓  IQ-TREE2 GTR+G, 1000 bootstrap replicates
    phylogenetic tree (.treefile)
```

#### Inputs

| Input | Description |
|-------|-------------|
| MAG FASTA files | *C. coli* MAG sequences in a specified directory |
| IQ-TREE2 binary | v2.3.6 — path hardcoded in script (update for your system) |
| Conda environments | Named environments containing Prokka and Panaroo |

#### Outputs

| Output | Description |
|--------|-------------|
| `prokka/<sample>/` | GFF annotation files per MAG |
| `core_gene_alignment.aln` | Core genome FASTA alignment |
| `<prefix>.treefile` | Maximum-likelihood phylogenetic tree |
| Bootstrap files | Ready for FigTree or iTOL |

#### Dependencies

```bash
conda install -c bioconda prokka
conda install -c bioconda panaroo
# IQ-TREE2 v2.3.6 — download binary from https://github.com/Cibiv/IQ-TREE
```

#### Running

```bash
bash scripts/shell/phylo_pipeline.sh
```

> Update conda environment names and the IQ-TREE2 binary path at the top of the script before running.

---

### `combine_blast_results.sh`

Merges BLAST tabular output files from multiple samples into a single consolidated file.

#### Inputs

Individual BLAST output files (tabular format, e.g., `-outfmt 6`).

#### Outputs

Single merged BLAST results file.

#### Dependencies

Standard Unix tools (`cat`, `awk`) — no additional installs required.

#### Running

```bash
bash scripts/shell/combine_blast_results.sh
```

---

### `rRNA_SNP_check_resistance.sh`

Checks for resistance-conferring SNPs in ribosomal RNA genes (e.g., 23S rRNA mutations associated with macrolide resistance in *Campylobacter*).

#### Inputs

| Input | Description |
|-------|-------------|
| Assembled contig / genome FASTA files | Sample sequences to check |
| Reference rRNA sequences | Target rRNA sequences for alignment |

#### Outputs

Report of SNPs detected in rRNA loci associated with antimicrobial resistance.

#### Dependencies

```bash
# Alignment tool — e.g., minimap2 or BLASTN
conda install -c bioconda minimap2
```

#### Running

```bash
bash scripts/shell/rRNA_SNP_check_resistance.sh
```

---

### `resistance_colocation_pipeline_ASC.sh`

Identifies co-location of multiple resistance genes on the same contig or genomic region in post-ASC enrichment samples.

#### Inputs

| Input | Description |
|-------|-------------|
| AMR annotation outputs | ABRicate / hAMRonization TSV files |
| Assembled contigs | FASTA files for the same samples |

#### Outputs

Table of contigs carrying multiple resistance genes, with gene positions and inter-gene distances.

#### Dependencies

```bash
pip install pandas
conda install -c bioconda abricate
```

#### Running

```bash
bash scripts/shell/resistance_colocation_pipeline_ASC.sh
```

---

### `resistance_colocation_pipeline_ASC 2.sh`

Alternative version of the ASC resistance co-location pipeline. Review alongside the primary script to identify differences in filtering thresholds or input formats before choosing which to use.

#### Running

```bash
bash "scripts/shell/resistance_colocation_pipeline_ASC 2.sh"
```

*(Note the quotes — required because the filename contains a space.)*

---

### `run_abricate_MAGs.sh`

Runs ABRicate AMR gene screening on all MAG FASTA files in a specified directory.

#### Inputs

| Input | Description |
|-------|-------------|
| MAG FASTA files | One FASTA file per MAG in a single directory |
| ABRicate database | e.g., `card`, `resfinder`, `megares`, `ncbi` |

#### Outputs

Per-MAG ABRicate TSV result files, one per MAG per database run.

#### Dependencies

```bash
conda install -c bioconda abricate
abricate --setupdb   # download databases if not already present
```

#### Running

```bash
bash scripts/shell/run_abricate_MAGs.sh
```

---

### `run_abricate_chunks.sh`

Processes ABRicate analysis in batches to handle large numbers of input files without exhausting memory or runtime limits.

#### Inputs

Input FASTA files — the script splits them into chunks internally.

#### Outputs

Per-chunk ABRicate TSV result files.

#### Dependencies

```bash
conda install -c bioconda abricate
```

#### Running

```bash
bash scripts/shell/run_abricate_chunks.sh
```

---

### `run_checkm2_gtdbtk.sl`

Local wrapper (despite the `.sl` extension) that runs CheckM2 and GTDB-Tk sequentially on a set of MAG FASTA files. May be a draft or local equivalent of the NeSI version in `scripts/slurm/`.

#### Inputs

MAG FASTA files in a specified directory.

#### Outputs

CheckM2 quality report; GTDB-Tk taxonomy assignments.

#### Dependencies

```bash
conda install -c bioconda checkm2
conda install -c bioconda gtdbtk
# GTDB-Tk requires the reference database — download separately:
# https://ecogenomics.github.io/GTDBTk/installing/index.html
```

#### Running

```bash
bash scripts/shell/run_checkm2_gtdbtk.sl
```
