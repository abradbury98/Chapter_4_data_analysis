# scripts/rmd

R Markdown notebooks for chapter-level analyses, statistical testing, and figure generation.  
Each notebook produces publication-quality PNG figures and summary tables as HTML or PDF output.

---

## Contents

| Notebook | Purpose |
|----------|---------|
| [`Chapter_5_read_level_analysis_diversity.Rmd`](#chapter_5_read_level_analysis_diversityrmd) | Read-level diversity and composition from Kraken2 BPM data |
| [`Chapter_5_contigs_level_bacterial_diversity.Rmd`](#chapter_5_contigs_level_bacterial_diversityrmd) | Contig-level bacterial diversity |
| [`Chapter_5_contigs_AMR_resistance_genes.Rmd`](#chapter_5_contigs_amr_resistance_genesrmd) | Contig-level AMR gene content |
| [`Chapter_5_read_level_chord_diagram_resistance_taxonomy.Rmd`](#chapter_5_read_level_chord_diagram_resistance_taxonomyrmd) | Chord diagrams linking resistance classes to bacterial genera |
| [`Chapter_5_MAG_quality_summary_graphs.Rmd`](#chapter_5_mag_quality_summary_graphsrmd) | MAG completeness, contamination, and taxonomy |
| [`Chapter_5_base_and_read_tracking.Rmd`](#chapter_5_base_and_read_trackingrmd) | Data throughput across pipeline steps |
| [`qPCR_graphs_enrichment_campy+salm_method_check.Rmd`](#qpcr_graphs_enrichment_campysalm_method_checkrmd) | qPCR enrichment method validation |

---

## Running Any Notebook

```bash
# Render to HTML (default)
Rscript -e "rmarkdown::render('scripts/rmd/<notebook>.Rmd')"

# Render to PDF
Rscript -e "rmarkdown::render('scripts/rmd/<notebook>.Rmd', output_format = 'pdf_document')"
```

**Before rendering:** open the notebook and update the file paths in the `setup` chunk at the top.  
**Output location:** figures and tables are saved to the directory specified inside each notebook.

---

## Notebooks

### `Chapter_5_read_level_analysis_diversity.Rmd`

Analyses bacterial community composition and diversity from Kraken2 read-level BPM data, comparing post-EVI and post-ASC samples across three enrichment protocols.

#### Inputs

Normalised Kraken2 BPM tables (TSV) at genus, family, and species level — produced by `kraken2_bpm_normalised_reads.py`.

#### Analyses

- Stacked bar charts and heatmaps of community composition
- Alpha diversity: Shannon index, Simpson index, species richness
- Beta diversity: Bray-Curtis PCoA and PERMANOVA
- Per-taxon Wilcoxon signed-rank tests
- Sensitivity analysis excluding a potential outlier sample

#### Outputs

PNG figures and statistical summary tables saved to a designated output directory.

#### Dependencies

```r
install.packages(c("tidyverse", "ggplot2", "vegan", "patchwork", "reshape2", "rmarkdown", "knitr"))
```

---

### `Chapter_5_contigs_level_bacterial_diversity.Rmd`

Mirrors the read-level diversity analysis using Kraken2 contig-level BPM data from assembled metagenomes.

#### Inputs

Normalised Kraken2 BPM tables (TSV) at contig level — produced by `kraken2_bpm_normalise_contigs.py`.

#### Analyses

Equivalent to the read-level notebook: composition, alpha/beta diversity, and statistical comparisons.

#### Outputs

PNG figures and statistical tables.

#### Dependencies

```r
install.packages(c("tidyverse", "ggplot2", "vegan", "patchwork", "reshape2", "rmarkdown", "knitr"))
```

---

### `Chapter_5_contigs_AMR_resistance_genes.Rmd`

Examines AMR gene content at the contig level, integrating ABRicate/hAMRonization results with Kraken2 contig taxonomy.

#### Inputs

| File | Description |
|------|-------------|
| Deduplicated AMR TSV files | Output from `hamronize_and_dedup.py` |
| Kraken2 contig taxonomy | `.out` / `.report` files from `kraken2_nt_assembly_contigs.sl` |

#### Analyses

AMR gene prevalence, resistance class summaries, cross-sample comparisons, and integration of taxonomy with resistance.

#### Outputs

Figures and tables describing AMR gene distribution across samples and enrichment conditions.

#### Dependencies

```r
install.packages(c("tidyverse", "ggplot2", "patchwork", "reshape2", "rmarkdown", "knitr"))
```

---

### `Chapter_5_read_level_chord_diagram_resistance_taxonomy.Rmd`

Produces chord diagrams and heatmaps linking resistance gene families to bacterial taxonomy at the read level, using paired Kraken2 and ABRicate outputs.

#### Inputs

TSV files matching `*_all_databases_with_taxonomy.tsv` from the directory:  
`/Volumes/Alice/2_POULTRY_PROCESSING/Analysis/Read_level_kraken-to-abricate_associations/`

Expected columns: `READ_ID`, `DATABASE`, `GENE`, `PRODUCT`, `RESISTANCE`, `%COVERAGE`, `%IDENTITY`, `TAXONOMY`

File naming convention: `CH01_ASC_CAT_all_databases_with_taxonomy.tsv` (sample ID and treatment extracted from filename).

Update the directory path in the setup chunk before rendering.

#### Analyses

- Bipartite chord diagrams — co-occurrence of resistance classes and bacterial genera
- Clustered heatmaps of resistance–taxonomy associations
- Summary tables of read-level AMR hits

#### Outputs

Chord diagram and heatmap figures (PNG/SVG).

#### Dependencies

```r
install.packages(c(
  "tidyverse", "ggplot2",
  "circlize",       # chord diagrams
  "RColorBrewer",   # colour palettes
  "viridis",        # perceptually uniform scales
  "pheatmap",       # clustered heatmaps
  "reshape2",
  "ggrepel",
  "knitr", "kableExtra",
  "rmarkdown"
))
```

---

### `Chapter_5_MAG_quality_summary_graphs.Rmd`

Summarises and visualises MAG quality metrics from CheckM2 and GTDB-Tk outputs, producing completeness/contamination scatter plots, taxonomy summaries, and quality tier breakdowns.

#### Inputs

Four TSV files — one per enrichment condition — containing merged CheckM2 and GTDB-Tk results:

| File | Condition |
|------|-----------|
| `ASC_RVS_summary_checkm2_gtdbtk.tsv` | Post-ASC, RVS buffer |
| `EVI_RVS_summary_checkm2_gtdbtk.tsv` | Post-EVI, RVS buffer |
| `ASC_BFBB_summary_checkm2_gtdbtk.tsv` | Post-ASC, Bolton BFBB buffer |
| `EVI_BFBB_summary_checkm2_gtdbtk.tsv` | Post-EVI, Bolton BFBB buffer |

Expected columns: `completeness`, `contamination`, `gtdbtk_classification`, `bin`, `tool`

#### Outputs

All figures saved at 600 dpi:

| Figure | Description |
|--------|-------------|
| Two-panel quality scatter | Completeness vs. contamination for all MAGs |
| All-MAGs combined view | Encoded by sampling location and enrichment type |
| Faceted grid | Completeness vs. contamination by condition |
| High-quality MAGs (≥50% complete) | With species labels |
| Strict HQ MAGs (≥95% complete, ≤10% contamination) | With jitter |
| Summary tables | MAG counts by quality tier and condition |

#### Dependencies

```r
install.packages(c("tidyverse", "ggplot2", "scales", "patchwork", "rmarkdown", "knitr"))
```

---

### `Chapter_5_base_and_read_tracking.Rmd`

Tracks read and base counts through each processing step to quantify data loss across the pipeline (host depletion → trimming → classification).

#### Inputs

Per-step read/base count log files or FASTQ files from each processing stage.

#### Analyses

Waterfall/Sankey-style summaries of reads retained and removed at each step.

#### Outputs

Figures and tables showing pipeline throughput per sample.

#### Dependencies

```r
install.packages(c("tidyverse", "ggplot2", "patchwork", "rmarkdown", "knitr"))
```

---

### `qPCR_graphs_enrichment_campy+salm_method_check.Rmd`

Analyses qPCR data for *Campylobacter* and *Salmonella* enrichment method validation, comparing Cq values across protocols.

#### Inputs

qPCR Cq value tables (CSV or Excel) for each enrichment method comparison.

#### Analyses

- Cq comparisons across enrichment protocols
- Standard curve evaluation
- Method equivalence testing

#### Outputs

qPCR result figures and a statistical summary for method validation.

#### Dependencies

```r
install.packages(c("tidyverse", "ggplot2", "patchwork", "rmarkdown", "knitr"))
```
