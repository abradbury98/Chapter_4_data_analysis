# scripts/rmd

R Markdown notebooks for chapter-level analyses and figure generation.

## Files

- `Chapter_5_MAG_quality_summary_graphs.Rmd`  
  MAG quality summaries and related figure generation.

- `Chapter_5_base_and_read_tracking.Rmd`  
  Base/read tracking analyses and summaries.

- `Chapter_5_contigs_AMR_resistance_genes.Rmd`  
  Contig-level AMR resistance gene analysis.

- `Chapter_5_contigs_level_bacterial_diversity.Rmd`  
  Contig-level bacterial diversity analyses.

- `Chapter_5_read_level_analysis_diversity.Rmd`  
  Read-level diversity analysis workflows.

- `Chapter_5_read_level_chord_diagram_resistance_taxonomy.Rmd`  
  Read-level resistance-taxonomy integration, including chord diagram visualization.

- `qPCR_graphs_enrichment_campy+salm_method_check.Rmd`  
  qPCR graphing and enrichment method checks for Campylobacter/Salmonella workflows.

## Typical usage

Render any notebook with:

```bash
Rscript -e "rmarkdown::render('scripts/rmd/<file>.Rmd')"
```

## Output

These notebooks typically produce HTML/PDF reports and publication-style figures/tables.
