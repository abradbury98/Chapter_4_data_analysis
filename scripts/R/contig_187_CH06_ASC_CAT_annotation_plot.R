## ============================================================
## Circular Contig Annotation Plot from Bakta GenBank Output
## Contig: contig_187 | Sample: ASC_CAT_CH06
## ============================================================

# install.packages(c("circlize", "dplyr", "stringr", "RColorBrewer"))
# BiocManager::install("Biostrings")

library(circlize)
library(dplyr)
library(stringr)
library(Biostrings)

# ── 1. CONFIGURATION ─────────────────────────────────────────

GBFF_FILE  <- "/Users/alicebradbury/Downloads/contig_187_CH06_ASC_CAT.fasta.gbff"
FASTA_FILE <- "/Users/alicebradbury/Downloads/contig_187_CH06_ASC_CAT.fasta.fna"

OUTPUT_PNG <- "/Users/alicebradbury/Desktop/contig_187_CH06_ASC_CAT_circular_annotation.png"
OUTPUT_PDF <- NULL

GC_WINDOW <- 2000
GC_STEP   <- 500

# ── 2. PARSE GENBANK FILE ────────────────────────────────────

parse_bakta_gbff <- function(gbff_path) {
  lines   <- readLines(gbff_path)
  records <- list()

  cur_contig  <- NULL
  cur_length  <- 0L
  in_features <- FALSE
  cur_feat    <- NULL
  cur_qual    <- NULL
  cur_val     <- ""

  flush_qualifier <- function() {
    if (!is.null(cur_qual) && !is.null(cur_feat)) {
      val <- trimws(cur_val)
      val <- gsub('^"|"$', "", val)
      cur_feat[[cur_qual]] <<- val
    }
    cur_qual <<- NULL
    cur_val  <<- ""
  }

  flush_feature <- function() {
    flush_qualifier()
    if (!is.null(cur_feat) && !is.null(cur_contig)) {
      records[[cur_contig]]$features <<- c(records[[cur_contig]]$features, list(cur_feat))
    }
    cur_feat <<- NULL
  }

  for (line in lines) {
    if (grepl("^LOCUS ", line)) {
      flush_feature()
      parts      <- strsplit(trimws(line), "\\s+")[[1]]
      cur_contig <- parts[2]
      cur_length <- as.integer(parts[3])
      records[[cur_contig]] <- list(name = cur_contig, length = cur_length, features = list())
      in_features <- FALSE
    }
    if (grepl("^FEATURES", line)) {
      flush_feature(); in_features <- TRUE; next
    }
    if (!in_features) next

    if (grepl("^     \\S", line)) {
      flush_feature()
      key_line   <- trimws(line)
      key        <- strsplit(key_line, "\\s+")[[1]][1]
      loc        <- strsplit(key_line, "\\s+")[[1]][2]
      complement <- grepl("complement", loc)
      loc_clean  <- gsub("complement\\(|join\\(|order\\(|\\)", "", loc)
      spans      <- strsplit(loc_clean, ",")[[1]]
      starts <- ends <- integer(length(spans))
      for (i in seq_along(spans)) {
        coords    <- as.integer(strsplit(gsub("[<>]", "", spans[i]), "\\.\\.")[[1]])
        starts[i] <- coords[1]
        ends[i]   <- if (length(coords) > 1) coords[2] else coords[1]
      }
      cur_feat <- list(
        type = key, start = min(starts), end = max(ends),
        strand = ifelse(complement, -1L, 1L), contig = cur_contig,
        gene = NA_character_, product = NA_character_, locus_tag = NA_character_
      )
    }

    if (grepl("^                     /", line)) {
      flush_qualifier()
      q_line <- trimws(sub("^\\s+", "", line))
      q_line <- sub("^/", "", q_line)
      eq_pos <- regexpr("=", q_line)
      if (eq_pos > 0) {
        cur_qual <- substr(q_line, 1, eq_pos - 1)
        cur_val  <- gsub('"', "", substr(q_line, eq_pos + 1, nchar(q_line)))
      } else {
        cur_qual <- q_line; cur_val <- "TRUE"
      }
    } else if (grepl("^                     [^\\/]", line) && !is.null(cur_qual)) {
      cur_val <- paste0(cur_val, " ", trimws(line))
    }

    if (grepl("^//", line)) flush_feature()
  }
  flush_feature()
  records
}

message("Parsing GenBank file …")
gbff_records <- parse_bakta_gbff(GBFF_FILE)

# ── 3. BUILD FEATURE TABLE ───────────────────────────────────

contig_sizes <- sapply(gbff_records, `[[`, "length") %>% sort(decreasing = TRUE)
top_contigs  <- names(contig_sizes)

extract_features <- function(records, contigs) {
  rows <- lapply(contigs, function(ctg) {
    feats <- records[[ctg]]$features
    if (length(feats) == 0) return(NULL)
    do.call(rbind, lapply(feats, function(f) {
      data.frame(
        contig  = ctg, type = f$type,
        start   = as.integer(f$start), end = as.integer(f$end),
        strand  = as.integer(f$strand),
        gene    = ifelse(is.null(f$gene),      NA, f$gene),
        product = ifelse(is.null(f$product),   NA, f$product),
        locus   = ifelse(is.null(f$locus_tag), NA, f$locus_tag),
        stringsAsFactors = FALSE
      )
    }))
  })
  bind_rows(rows)
}

feats_all <- extract_features(gbff_records, top_contigs)

genes <- feats_all %>%
  filter(type %in% c("CDS", "tRNA", "rRNA", "ncRNA", "tmRNA", "regulatory")) %>%
  mutate(label = coalesce(gene, locus))

# ── 4. FUNCTIONAL CATEGORIES ─────────────────────────────────

assign_category <- function(type, product, gene_name) {
  p <- tolower(coalesce(product, ""))
  t <- tolower(type)
  g <- tolower(coalesce(gene_name, ""))

  case_when(
    t %in% c("trna", "rrna", "tmrna", "ncrna")
      ~ "RNA gene",

    # Mobile elements — transposases, insertion sequences, integrases
    str_detect(p, "transposase|insertion element|is\\d|is[0-9]|integrase|resolvase|recombinase|relaxase|mobilization") |
      str_detect(g, "^tnp$|^tnpa$|^tnpb$")
      ~ "Mobile element",

    # Conjugal transfer — tra/pil conjugation machinery
    str_detect(p, "conjugal|conjugative|mating pair|type iv secretion|type ivb pilus|prepilin|pilus protein pil|pilp|piln|pilo|pilm|pilk|pilj|pilr|pils") |
      str_detect(g, "^tra[abcdefghijklmnopqrstuvwxyz]$|^pil[a-z]")
      ~ "Conjugal transfer",

    # Plasmid replication
    str_detect(p, "replication initiator|repa|repb|replication protein|plasmid replication") |
      str_detect(g, "^rep[abcz]$")
      ~ "Plasmid replication",

    # AMR
    str_detect(p, "resist|efflux|beta-lactam|carbapenem|aminoglycoside|tetracycline|quinolone|macrolide|colistin|sulfonamide") |
      str_detect(g, "^sul|^tet|^aph|^aac|^ant|^aad|^mcr|^bla|^cfr|^msr|^erm|^van")
      ~ "AMR",

    # Virulence / toxins
    str_detect(p, "virulence|toxin|colicin|hemolysin|fimbriae|fimbria|adhesin|invasion|flagell|type iii|type vi|siderophore|iron|addiction module")
      ~ "Virulence",

    str_detect(p, "ribosomal|translation|trna ligase|aminoacyl")
      ~ "Translation",

    str_detect(p, "dna polymerase|dna repair|helicase|gyrase|topoisomerase|dnaa|dnab|exonuclease")
      ~ "DNA replication/repair",

    str_detect(p, "rna polymerase|transcription|sigma factor|regulator|two-component|repressor|activator|nusg|antitermination")
      ~ "Transcription/regulation",

    str_detect(p, "transport|permease|abc transporter|efflux|outer membrane protein|porin|channel|atp-binding protein")
      ~ "Transport",

    str_detect(p, "metabol|synthase|synthetase|dehydrogenase|reductase|kinase|isomerase|transferase|oxidase|phosphatase|mutase|phosphohexomutase")
      ~ "Metabolism",

    str_detect(p, "hypothetical|uncharacterized|putative|unknown|predicted|duf")
      ~ "Hypothetical",

    TRUE ~ "Other"
  )
}

genes <- genes %>% mutate(category = assign_category(type, product, gene))

cat_colours <- c(
  "RNA gene"                 = "#E41A1C",  # vivid red
  "Mobile element"           = "#FF7F00",  # orange
  "Conjugal transfer"        = "#1F78B4",  # medium blue
  "Plasmid replication"      = "#33A02C",  # dark forest green
  "AMR"                      = "#C51B7D",  # deep magenta  (was near-identical red to RNA gene)
  "Virulence"                = "#984EA3",  # purple
  "Translation"              = "#FFD700",  # gold          (was green, too similar to Plasmid replication)
  "DNA replication/repair"   = "#8B4513",  # sienna brown  (was near-identical blue to Conjugal transfer)
  "Transcription/regulation" = "#17BECF",  # cyan
  "Transport"                = "#F781BF",  # light pink
  "Metabolism"               = "#66C2A5",  # sage/mint     (was periwinkle, drifting into blue family)
  "Hypothetical"             = "#B3B3B3",  # light grey
  "Other"                    = "#666666"   # dark grey
)

genes$colour <- cat_colours[genes$category]
genes$colour[is.na(genes$colour)] <- "#666666"

# ── 5. GC CONTENT & SKEW ────────────────────────────────────

compute_gc_tracks <- function(fasta_path, contigs, window = 2000, step = 500) {
  seqs        <- readDNAStringSet(fasta_path)
  names(seqs) <- sub(" .*", "", names(seqs))

  result <- lapply(contigs, function(ctg) {
    if (!ctg %in% names(seqs)) return(NULL)
    s   <- seqs[[ctg]]
    len <- length(s)
    if (len < window) return(NULL)
    starts <- seq(1, len - window + 1, by = step)
    gc_val <- skew <- midpts <- numeric(length(starts))
    for (i in seq_along(starts)) {
      sub_s   <- subseq(s, starts[i], min(starts[i] + window - 1, len))
      freq    <- letterFrequency(sub_s, c("G", "C", "A", "T"))
      g <- freq["G"]; c <- freq["C"]
      gc_val[i] <- (g + c) / sum(freq)
      skew[i]   <- if ((g + c) > 0) (g - c) / (g + c) else 0
      midpts[i] <- starts[i] + window / 2
    }
    data.frame(contig = ctg, pos = midpts, gc = gc_val, skew = skew)
  })
  bind_rows(result)
}

message("Computing GC tracks …")
gc_data <- tryCatch(
  compute_gc_tracks(FASTA_FILE, top_contigs, GC_WINDOW, GC_STEP),
  error = function(e) { message("GC track skipped: ", e$message); NULL }
)

# ── 6. GENOME LAYOUT ─────────────────────────────────────────

layout <- data.frame(
  contig = top_contigs,
  len    = as.integer(contig_sizes[top_contigs]),
  stringsAsFactors = FALSE
) %>% mutate(
  offset = cumsum(c(0, head(len, -1))),
  total  = offset + len
)

total_genome <- max(layout$total)

coords_for <- function(df, layout_df) {
  left_join(df, layout_df %>% select(contig, offset), by = "contig") %>%
    mutate(abs_start = start + offset, abs_end = end + offset)
}

genes_abs <- coords_for(genes, layout)
gc_abs <- if (!is.null(gc_data)) {
  left_join(gc_data, layout %>% select(contig, offset), by = "contig") %>%
    mutate(abs_pos = pos + offset)
} else NULL

# ── 7. LABELS — AMR, Mobile element, Conjugal transfer, Plasmid replication
label_cats <- c("AMR", "Mobile element", "Conjugal transfer", "Plasmid replication")

label_genes <- genes_abs %>%
  filter(category %in% label_cats, type == "CDS") %>%
  mutate(
    gene_label = case_when(
      !is.na(gene)                                                            ~ gene,
      str_detect(tolower(coalesce(product, "")), "transposase")               ~ "Tnp",
      str_detect(tolower(coalesce(product, "")), "relaxase")                  ~ "Relaxase",
      str_detect(tolower(coalesce(product, "")), "resolvase")                 ~ "Resolvase",
      str_detect(tolower(coalesce(product, "")), "recombinase")               ~ "Recombinase",
      str_detect(tolower(coalesce(product, "")), "prepilin|thin prepilin")    ~ "Prepilin",
      str_detect(tolower(coalesce(product, "")), "atp-binding")               ~ "ATPase",
      TRUE                                                                     ~ locus
    ),
    mid_abs = (abs_start + abs_end) / 2
  ) %>%
  # Keep one representative per closely-spaced cluster to avoid label clutter
  arrange(mid_abs) %>%
  group_by(gene_label) %>%
  dplyr::slice(1) %>%
  ungroup()

# ── 8. STAGGER HELPER ────────────────────────────────────────

stagger_levels <- function(positions, min_gap, n_levels = 6) {
  ord        <- order(positions)
  sorted_pos <- positions[ord]
  lv         <- integer(length(positions))
  last_pos   <- rep(-Inf, n_levels)
  for (i in seq_along(sorted_pos)) {
    chosen <- NA
    for (l in seq_len(n_levels)) {
      if (sorted_pos[i] - last_pos[l] >= min_gap) { chosen <- l; break }
    }
    if (is.na(chosen)) chosen <- ((i - 1) %% n_levels) + 1
    lv[ord[i]]       <- chosen
    last_pos[chosen] <- sorted_pos[i]
  }
  lv
}

# ── 9. DRAW CIRCULAR PLOT ────────────────────────────────────

if (!is.null(OUTPUT_PNG)) {
  png(OUTPUT_PNG, width = 14, height = 14, units = "in", res = 600, bg = "white")
} else if (!is.null(OUTPUT_PDF)) {
  pdf(OUTPUT_PDF, width = 14, height = 14, useDingbats = FALSE)
}

par(mar = c(1, 1, 1, 1))

circos.clear()
circos.par(
  "start.degree"            = 90,
  "clock.wise"              = TRUE,
  "track.height"            = 0.08,
  "gap.degree"              = 2,
  "cell.padding"            = c(0, 0, 0, 0),
  "points.overflow.warning" = FALSE,
  "canvas.xlim"             = c(-1.65, 1.65),
  "canvas.ylim"             = c(-1.65, 1.65)
)

circos.initialize(factors = "genome", xlim = c(0, total_genome))

sub_fwd <- genes_abs %>% filter(strand ==  1L, type == "CDS")
sub_rev <- genes_abs %>% filter(strand == -1L, type == "CDS")

LABEL_GAP <- total_genome * 0.05
y_out     <- c(1.20, 1.45, 1.70, 1.95, 2.20, 2.45)

fwd_labs <- label_genes %>% filter(strand ==  1L) %>% mutate(lty = 1L)
rev_labs <- label_genes %>% filter(strand == -1L) %>% mutate(lty = 2L)
all_labs <- bind_rows(fwd_labs, rev_labs) %>% arrange(mid_abs)

if (nrow(all_labs) > 0) {
  all_labs$level <- stagger_levels(all_labs$mid_abs, LABEL_GAP, n_levels = 6)
  rev_idx <- all_labs$lty == 2L
  all_labs$level[rev_idx] <- pmax(all_labs$level[rev_idx], 3L)
}

# ── Track 1: Forward CDS + labels ────────────────────────────
circos.track(
  ylim = c(0, 1), bg.col = "#F0F0F0", bg.border = NA, track.height = 0.09,
  panel.fun = function(x, y) {
    for (i in seq_len(nrow(sub_fwd)))
      circos.rect(sub_fwd$abs_start[i], 0.05, sub_fwd$abs_end[i], 0.95,
                  col = sub_fwd$colour[i], border = NA)
    hl <- all_labs[all_labs$lty == 1L, ]
    for (i in seq_len(nrow(hl)))
      circos.rect(hl$abs_start[i], 0.02, hl$abs_end[i], 0.98,
                  col = hl$colour[i], border = "white", lwd = 0.8)
    if (nrow(all_labs) > 0) {
      for (i in seq_len(nrow(all_labs))) {
        pos  <- all_labs$mid_abs[i]
        col  <- all_labs$colour[i]
        lv   <- all_labs$level[i]
        ylbl <- y_out[lv]
        lty  <- all_labs$lty[i]
        circos.segments(pos, 1.0, pos, ylbl, col = col, lwd = 1.2, lty = lty)
        circos.text(pos, ylbl, all_labs$gene_label[i],
                    facing = "clockwise", niceFacing = TRUE,
                    cex = 0.88, col = col,
                    font = ifelse(all_labs$category[i] == "AMR", 4L, 3L),
                    adj = c(0, 0.5))
      }
    }
  }
)

# ── Track 2: Reverse CDS ─────────────────────────────────────
circos.track(
  ylim = c(0, 1), bg.col = "#E8E8E8", bg.border = NA, track.height = 0.09,
  panel.fun = function(x, y) {
    for (i in seq_len(nrow(sub_rev)))
      circos.rect(sub_rev$abs_start[i], 0.05, sub_rev$abs_end[i], 0.95,
                  col = sub_rev$colour[i], border = NA)
    rl <- all_labs[all_labs$lty == 2L, ]
    for (i in seq_len(nrow(rl)))
      circos.rect(rl$abs_start[i], 0.02, rl$abs_end[i], 0.98,
                  col = rl$colour[i], border = "white", lwd = 0.8)
  }
)

# ── Track 3: GC skew ─────────────────────────────────────────
if (!is.null(gc_abs) && nrow(gc_abs) > 0) {
  skew_lim  <- max(abs(gc_abs$skew), na.rm = TRUE)
  if (skew_lim == 0) skew_lim <- 0.1
  skew_cols <- ifelse(gc_abs$skew >= 0, "#d7191c", "#1a9641")
  circos.track(
    ylim = c(-skew_lim, skew_lim), bg.col = "#EBEBEB", bg.border = NA,
    panel.fun = function(x, y) {
      for (i in seq_along(gc_abs$abs_pos))
        circos.lines(c(gc_abs$abs_pos[i], gc_abs$abs_pos[i]),
                     c(0, gc_abs$skew[i]), col = skew_cols[i], lwd = 0.5)
      circos.lines(CELL_META$xlim, c(0, 0), col = "grey50", lwd = 0.5, lty = 2)
    }
  )
}

# ── Track 4: RNA genes (skipped if none) ─────────────────────
rna_genes <- genes_abs %>% filter(type %in% c("tRNA", "rRNA", "tmRNA", "ncRNA"))
if (nrow(rna_genes) > 0) {
  circos.track(
    ylim = c(0, 1), bg.col = "#F8F8F8", bg.border = NA, track.height = 0.04,
    panel.fun = function(x, y) {
      for (i in seq_len(nrow(rna_genes)))
        circos.rect(rna_genes$abs_start[i], 0.05,
                    rna_genes$abs_end[i],   0.95, col = "#E41A1C", border = NA)
    }
  )
}

# ── Spacer ───────────────────────────────────────────────────
circos.track(ylim = c(0, 1), bg.col = NA, bg.border = NA, track.height = 0.03)

# ── Track 5: kb scale ring ───────────────────────────────────
circos.track(
  ylim = c(0, 1), track.height = 0.05, bg.border = "grey70", bg.col = "white",
  panel.fun = function(x, y) {
    tick_at <- seq(0, total_genome, by = 10000)
    circos.axis(
      h                 = "bottom",
      major.at          = tick_at,
      labels            = paste0(round(tick_at / 1e3, 0), " kb"),
      labels.cex        = 0.75,
      labels.facing     = "clockwise",
      labels.niceFacing = TRUE,
      col               = "grey30",
      labels.col        = "grey30",
      lwd               = 0.8,
      minor.ticks       = 4
    )
  }
)

# ── Centre text ──────────────────────────────────────────────
contig_size_kb <- round(sum(layout$len) / 1e3, 1)
n_cds   <- sum(genes$type == "CDS")
n_trna  <- sum(genes$type == "tRNA")
n_amr   <- sum(genes$category == "AMR")
n_mob   <- sum(genes$category == "Mobile element")
n_conj  <- sum(genes$category == "Conjugal transfer")
gc_pct  <- if (!is.null(gc_abs) && nrow(gc_abs) > 0) round(mean(gc_abs$gc) * 100, 1) else "n/a"

symbols(0, 0, circles = 0.44, inches = FALSE, bg = "grey15", fg = NA, add = TRUE)
text(0,  0.22, "contig_187",                       cex = 1.40, font = 3, col = "white")
text(0,  0.09, "ASC_CAT_CH06",                     cex = 1.10, font = 1, col = "grey80")
text(0, -0.05, paste0(n_cds, " CDSs | ", n_trna, " tRNAs"),
               cex = 0.95, font = 1, col = "grey70")
text(0, -0.17, paste0(n_amr, " AMR | ", n_mob, " mobile | ", n_conj, " conjugal"),
               cex = 0.88, font = 1, col = "grey70")
text(0, -0.30, paste0(contig_size_kb, " kb  |  GC ", gc_pct, "%"),
               cex = 1.00, font = 1, col = "grey70")

# ── Legends ──────────────────────────────────────────────────
LEG_CEX   <- 1.05
LEG_INTER <- 0.70
USR       <- par("usr")

cats_present <- unique(genes$category)
legend_cols  <- cat_colours[intersect(names(cat_colours), cats_present)]

# Right: functional categories
legend(
  x = USR[2], y = USR[3] + 0.05, xjust = 1, yjust = 0,
  legend = names(legend_cols), fill = legend_cols, border = NA,
  bty = "n", cex = LEG_CEX, text.font = 1,
  title = "Functional category", title.font = 2, title.col = "grey15",
  y.intersp = LEG_INTER
)

# Left: strand key
legend(
  x = USR[1], y = USR[3] + 0.05 + (2 + 1.8) * LEG_INTER * LEG_CEX * 0.182,
  xjust = 0, yjust = 0,
  legend = c("(+) Forward strand", "(-) Reverse strand"),
  col = c("grey35", "grey35"), lwd = 1.5, lty = c(1, 2),
  bty = "n", cex = LEG_CEX * 0.90,
  title = "Strand (line style)", title.font = 2, title.col = "grey15",
  y.intersp = LEG_INTER
)

circos.clear()

if (!is.null(OUTPUT_PNG) || !is.null(OUTPUT_PDF)) {
  dev.off()
  message("Plot saved: ", coalesce(OUTPUT_PNG, OUTPUT_PDF))
}
