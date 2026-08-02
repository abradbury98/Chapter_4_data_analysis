# scripts/R

Standalone R scripts for Nanopore sequencing quality visualisation and rarefaction analysis.  
Run directly with `Rscript` — no notebook rendering required.

---

## Contents

| Script | Purpose |
|--------|---------|
| [`pore_activity_plot.R`](#pore_activity_plotr) | Nanopore pore state activity across six sequencing runs |
| [`rarefaction_curves.R`](#rarefaction_curvesr) | Rarefaction curves from Kraken2 reports and AMRplusplus data |
| [`rarefaction_curves_all_conditions.Rmd`](#rarefaction_curves_all_conditionsrmd) | R Markdown version covering all enrichment conditions |

---

## Scripts

### `pore_activity_plot.R`

Generates publication-quality stacked bar charts of Nanopore pore state activity, comparing post-EVI and post-ASC enrichment conditions with two buffer treatments (Bolton base + CAT; BPW + RVS) across six sequencing runs.

#### Inputs

| File | Description |
|------|-------------|
| `pore_activity_data.csv` | Percentage-based pore state data per time point |
| `pore_scan_data.csv` | Absolute pore counts at 1.5-hour intervals |

Update the file paths near the top of the script before running.

#### Outputs

| File | Description |
|------|-------------|
| `pore_activity_all_runs.png` | Stacked bar chart — percentage of pores by state over run time |
| `pore_scan_all_runs.png` | Stacked bar chart — absolute pore counts over run time |

Both outputs: 300 dpi, 12 × 10 inches, faceted by run type, using a MinKNOW-compatible colour scheme for five pore states (Sequencing, Unavailable, etc.).

#### Dependencies

```r
install.packages(c("ggplot2", "dplyr"))
```

#### Running

```bash
Rscript scripts/R/pore_activity_plot.R
```

---

### `rarefaction_curves.R`

Plots rarefaction curves to assess whether sequencing depth was sufficient to capture community diversity. Reads Kraken2 species-level reports and AMRplusplus rarefaction data and generates combined figures.

#### Inputs

**Kraken2 reports** — tab-delimited `.report` files with columns: `pct`, `reads_covered`, `reads_direct`, `rank`, `taxid`, `name`. The script filters for species-level rows (`rank == "S"`).  
Expected directory structure:

```
~/Desktop/Kraken2_contigs/
├── ASC_CAT/    *.report
├── ASC_RVS/    *.report
├── EVI_CAT/    *.report
└── EVI_RVS/    *.report
```

**AMRplusplus rarefaction data** — tab-separated `.tsv` files named by barcode and classification level (e.g., `barcode01.gene.tsv`) with columns: read counts, unique ARG counts.  
Expected location: `~/Desktop/AMRplusplus_ASC_BFBB_results/ResistomeAnalysis/Rarefaction/Counts/`

Update all paths and sample labels at the top of the script before running.

#### Outputs

Rarefaction curve figures (PNG) saved to `~/Desktop/figures/`.

#### Dependencies

```r
install.packages(c("tidyverse", "vegan", "ggplot2", "patchwork"))
```

#### Running

```bash
Rscript scripts/R/rarefaction_curves.R
```

---

### `rarefaction_curves_all_conditions.Rmd`

R Markdown version of the rarefaction analysis covering all enrichment conditions. Combines code, narrative, and output figures in a single rendered document.

#### Inputs

Same input format as `rarefaction_curves.R`, extended across all sample conditions. Update file paths in the setup chunk.

#### Outputs

Rendered HTML (or PDF) report with embedded rarefaction plots.

#### Dependencies

```r
install.packages(c("tidyverse", "vegan", "ggplot2", "patchwork", "rmarkdown", "knitr"))
```

#### Running

```bash
# Render to HTML
Rscript -e "rmarkdown::render('scripts/R/rarefaction_curves_all_conditions.Rmd')"

# Render to PDF
Rscript -e "rmarkdown::render('scripts/R/rarefaction_curves_all_conditions.Rmd', output_format = 'pdf_document')"
```
