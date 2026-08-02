# scripts/R

R and R Markdown files for plotting and rarefaction-focused analyses.

## Files

- `pore_activity_plot.R`  
  Generates pore activity plots (likely for sequencing run performance/activity over time).

- `rarefaction_curves.R`  
  Produces rarefaction curves, typically to assess sequencing depth versus observed diversity/features.

- `rarefaction_curves_all_conditions.Rmd`  
  R Markdown workflow for rarefaction analyses across all conditions, combining code, outputs, and narrative.

## Typical usage

Run `.R` scripts with:

```bash
Rscript scripts/R/<script_name>.R
```

Render the `.Rmd` report with:

```bash
Rscript -e "rmarkdown::render('scripts/R/rarefaction_curves_all_conditions.Rmd')"
```

## Dependencies

Check script headers for required packages (commonly tidyverse/ggplot2/rmarkdown-style dependencies).
