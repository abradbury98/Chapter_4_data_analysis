# =============================================================================
# Rarefaction Curves: Taxonomic Diversity (Kraken2) & ARGs (AMRplusplus)
# =============================================================================

library(tidyverse)
library(vegan)
library(ggplot2)
library(patchwork)

# --- Paths -------------------------------------------------------------------
kraken_base  <- "~/Desktop/Kraken2_contigs"
amr_base     <- "~/Desktop/AMRplusplus_ASC_BFBB_results/ResistomeAnalysis/Rarefaction/Counts"
out_dir      <- "~/Desktop/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Condition folders in Kraken2 directory (edit to add/remove)
kraken_conditions <- c("ASC_CAT", "ASC_RVS", "EVI_CAT", "EVI_RVS")

# Sample labels for AMRplusplus barcodes
amr_labels <- c(
  "all_barcode01_porechop_nanofilt_filtered" = "Barcode 01",
  "all_barcode03_porechop_nanofilt_filtered" = "Barcode 03",
  "all_barcode06_porechop_nanofilt_filtered" = "Barcode 06"
)

# AMR classification level to plot (gene / class / group / mech / type)
amr_level <- "gene"

# Colour palette (one colour per sample/condition)
palette_kraken <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
  "#FF7F00", "#A65628", "#F781BF", "#999999",
  "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
  "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",
  "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
  "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928"
)
palette_amr    <- c("#1B9E77", "#D95F02", "#7570B3")


# =============================================================================
# 1. TAXONOMIC RAREFACTION FROM KRAKEN2 REPORTS
# =============================================================================

#' Parse a single Kraken2 report and return species-level taxon counts
parse_kraken_report <- function(path) {
  df <- read.delim(path, header = FALSE, sep = "\t", quote = "",
                   col.names = c("pct", "reads_covered", "reads_direct",
                                 "rank", "taxid", "name"))
  # Keep only species-level rows with at least 1 direct read
  df %>%
    filter(rank == "S", reads_direct > 0) %>%
    mutate(name = trimws(name)) %>%
    select(name, reads_direct)
}

#' Build a community matrix (samples × species) from a condition folder
build_community_matrix <- function(condition, base = kraken_base) {
  folder  <- file.path(base, condition)
  files   <- list.files(folder, pattern = "\\.report$", full.names = TRUE)
  samples <- sub("_contigs_kraken2\\.report$", "", basename(files))
  
  all_taxa <- map(files, parse_kraken_report)
  all_names <- unique(unlist(map(all_taxa, "name")))
  
  mat <- map_dfc(all_taxa, function(d) {
    v <- setNames(integer(length(all_names)), all_names)
    v[d$name] <- d$reads_direct
    as.integer(v)
  }) %>%
    setNames(samples) %>%
    t()                   # rows = samples, cols = species
  
  colnames(mat) <- all_names
  mat
}

# Build rarefaction curves per condition and collect into a tidy data frame
kraken_rare_df <- map_dfr(kraken_conditions, function(cond) {
  message("Processing Kraken2: ", cond)
  mat <- build_community_matrix(cond)
  
  # Rarefy to the minimum library size in this condition
  min_reads <- min(rowSums(mat))
  message("  min reads = ", min_reads, " (used as rarefaction depth)")
  
  rare_out <- rarecurve(mat,
                        step     = max(1, floor(min_reads / 50)),
                        sample   = min_reads,
                        label    = FALSE,
                        tidy     = FALSE)
  
  # rarecurve returns a named list; convert to tidy
  map_dfr(seq_along(rare_out), function(i) {
    curve   <- rare_out[[i]]
    reads   <- as.integer(attr(curve, "Subsample"))
    tibble(
      condition = cond,
      sample    = rownames(mat)[i],
      reads     = reads,
      richness  = as.integer(curve)
    )
  })
})

# --- Plot --------------------------------------------------------------------
p_kraken <- ggplot(kraken_rare_df,
                   aes(x = reads, y = richness,
                       colour = sample, linetype = condition,
                       group = interaction(condition, sample))) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = palette_kraken) +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash", "twodash")) +
  labs(
    title    = "Taxonomic Rarefaction Curves (Kraken2 species-level)",
    x        = "Number of reads sampled",
    y        = "Species richness",
    colour   = "Sample",
    linetype = "Condition"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 13))

ggsave(file.path(out_dir, "rarefaction_taxonomic_kraken2.png"),
       p_kraken, width = 10, height = 6, dpi = 300)
message("Saved: rarefaction_taxonomic_kraken2.png")


# =============================================================================
# 2. ARG RAREFACTION FROM AMRPLUSPLUS PRE-COMPUTED COUNTS
# =============================================================================

#' Read a single AMRplusplus rarefaction .tsv (reads, unique_ARGs)
read_amr_rare <- function(path, sample_name, level) {
  df <- read.delim(path, header = FALSE, sep = "\t",
                   col.names = c("reads", "unique_args"))
  df$sample <- sample_name
  df$level  <- level
  df
}

# Find all rarefaction files for the chosen level
amr_files <- list.files(amr_base,
                        pattern = paste0("\\.", amr_level, "\\.tsv$"),
                        full.names = TRUE)

if (length(amr_files) == 0) {
  stop("No AMRplusplus rarefaction files found for level '", amr_level,
       "' in:\n  ", amr_base)
}

amr_rare_df <- map_dfr(amr_files, function(f) {
  raw_name    <- sub(paste0("\\.", amr_level, "\\.tsv$"), "", basename(f))
  sample_label <- if (raw_name %in% names(amr_labels)) amr_labels[[raw_name]] else raw_name
  read_amr_rare(f, sample_label, amr_level)
})

# --- Plot --------------------------------------------------------------------
p_amr <- ggplot(amr_rare_df,
                aes(x = reads, y = unique_args,
                    colour = sample, group = sample)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_colour_manual(values = palette_amr) +
  scale_y_continuous(breaks = scales::pretty_breaks(6)) +
  labs(
    title = paste0("ARG Rarefaction Curves — AMRplusplus (",
                   str_to_title(amr_level), " level)"),
    x      = "Number of reads sampled",
    y      = paste("Unique ARGs (", amr_level, "level)"),
    colour = "Sample"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 13))

ggsave(file.path(out_dir, "rarefaction_ARG_amrplusplus.png"),
       p_amr, width = 9, height = 5.5, dpi = 300)
message("Saved: rarefaction_ARG_amrplusplus.png")


# =============================================================================
# 3. COMBINED PANEL (taxonomy + ARGs side by side)
# =============================================================================

p_combined <- p_kraken + p_amr +
  plot_annotation(
    title   = "Rarefaction Curves: Taxonomic Diversity & Antimicrobial Resistance Genes",
    theme   = theme(plot.title = element_text(face = "bold", size = 14,
                                              hjust = 0.5))
  )

ggsave(file.path(out_dir, "rarefaction_combined.png"),
       p_combined, width = 18, height = 6.5, dpi = 300)
message("Saved: rarefaction_combined.png")


# =============================================================================
# 4. OPTIONAL: Per-condition Kraken2 plots (one panel per condition)
# =============================================================================

p_kraken_facet <- ggplot(kraken_rare_df,
                         aes(x = reads, y = richness,
                             colour = sample, group = sample)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ condition, scales = "free_x") +
  scale_colour_manual(values = palette_kraken) +
  labs(
    title  = "Taxonomic Rarefaction by Condition (Kraken2)",
    x      = "Reads sampled",
    y      = "Species richness",
    colour = "Sample"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 13))

ggsave(file.path(out_dir, "rarefaction_taxonomic_by_condition.png"),
       p_kraken_facet, width = 12, height = 9, dpi = 300)
message("Saved: rarefaction_taxonomic_by_condition.png")


# =============================================================================
# 5. OPTIONAL: All AMRplusplus levels on one plot
# =============================================================================

amr_levels_all <- c("gene", "class", "group", "mech", "type")

amr_all_df <- map_dfr(amr_levels_all, function(lvl) {
  files <- list.files(amr_base, pattern = paste0("\\.", lvl, "\\.tsv$"),
                      full.names = TRUE)
  map_dfr(files, function(f) {
    raw_name     <- sub(paste0("\\.", lvl, "\\.tsv$"), "", basename(f))
    sample_label <- if (raw_name %in% names(amr_labels)) amr_labels[[raw_name]] else raw_name
    read_amr_rare(f, sample_label, lvl)
  })
})

p_amr_all_levels <- ggplot(amr_all_df,
                           aes(x = reads, y = unique_args,
                               colour = sample, group = sample)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2, alpha = 0.6) +
  facet_wrap(~ level, scales = "free_y",
             labeller = labeller(level = str_to_title)) +
  scale_colour_manual(values = palette_amr) +
  labs(
    title  = "ARG Rarefaction — All AMRplusplus Classification Levels",
    x      = "Number of reads sampled",
    y      = "Unique ARGs",
    colour = "Sample"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 13))

ggsave(file.path(out_dir, "rarefaction_ARG_all_levels.png"),
       p_amr_all_levels, width = 14, height = 10, dpi = 300)
message("Saved: rarefaction_ARG_all_levels.png")

message("\nAll done. Figures saved to: ", out_dir)