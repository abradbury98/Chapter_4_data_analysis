## ============================================================
## Circular Genome Annotation Plot from Bakta GenBank Output
## Genome: ASC_CAT_CH06 Escherichia fergusonii
## ============================================================

# Install packages if needed (run once):
# install.packages(c("circlize", "dplyr", "stringr", "RColorBrewer", "BiocManager"))
# BiocManager::install("Biostrings")

library(circlize)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(Biostrings)
library(readxl)

# ── 1. CONFIGURATION ────────────────────────────────────────

GBFF_FILE   <- "/Users/alicebradbury/Downloads/ASC_CAT_CH06_Escherichia_fergusonii.gbff"
FASTA_FILE  <- "/Users/alicebradbury/Downloads/ASC_CAT_CH06_Escherichia_fergusonii.fna"
ARG_FILE    <- "/Users/alicebradbury/Downloads/ASC_CAT_CH06_Escherichia_fergusonii_ARG_summary.csv"
VIR_FILE    <- "/Users/alicebradbury/Downloads/ASC_CAT_CH06_Escherichia_fergusonii_virulence_summary.xlsx"

# Show all contigs so no AMR/virulence genes are missed (e.g. blaCTX-M-1 is on a small contig)
TOP_N_CONTIGS <- Inf

# GC window size (bp)
GC_WINDOW   <- 5000
GC_STEP     <- 1000

# Output PDF (set to NULL to plot interactively)
OUTPUT_PDF  <- "/Users/alicebradbury/Desktop/ASC_CAT_CH06_circular_annotation.pdf"
OUTPUT_PNG  <- "/Users/alicebradbury/Desktop/ASC_CAT_CH06_circular_annotation.png"

# ── 2. PARSE GENBANK FILE ───────────────────────────────────

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
    # New locus
    if (grepl("^LOCUS ", line)) {
      flush_feature()
      parts      <- strsplit(trimws(line), "\\s+")[[1]]
      cur_contig <- parts[2]
      cur_length <- as.integer(parts[3])
      records[[cur_contig]] <- list(name = cur_contig, length = cur_length, features = list())
      in_features <- FALSE
    }

    # Start of FEATURES block
    if (grepl("^FEATURES", line)) {
      flush_feature()
      in_features <- TRUE
      next
    }

    if (!in_features) next

    # Feature key line: 5-space indent + feature type + location
    if (grepl("^     \\S", line)) {
      flush_feature()
      key_line <- trimws(line)
      key   <- strsplit(key_line, "\\s+")[[1]][1]
      loc   <- strsplit(key_line, "\\s+")[[1]][2]

      complement <- grepl("complement", loc)
      loc_clean  <- gsub("complement\\(|join\\(|order\\(|\\)", "", loc)
      spans      <- strsplit(loc_clean, ",")[[1]]

      starts <- ends <- integer(length(spans))
      for (i in seq_along(spans)) {
        coords     <- as.integer(strsplit(gsub("[<>]", "", spans[i]), "\\.\\.")[[1]])
        starts[i]  <- coords[1]
        ends[i]    <- if (length(coords) > 1) coords[2] else coords[1]
      }

      cur_feat <- list(
        type       = key,
        start      = min(starts),
        end        = max(ends),
        strand     = ifelse(complement, -1L, 1L),
        contig     = cur_contig,
        gene       = NA_character_,
        product    = NA_character_,
        locus_tag  = NA_character_
      )
    }

    # Qualifier lines: 21-space indent
    if (grepl("^                     /", line)) {
      flush_qualifier()
      q_line   <- trimws(sub("^\\s+", "", line))
      q_line   <- sub("^/", "", q_line)
      eq_pos   <- regexpr("=", q_line)
      if (eq_pos > 0) {
        cur_qual <- substr(q_line, 1, eq_pos - 1)
        cur_val  <- gsub('"', "", substr(q_line, eq_pos + 1, nchar(q_line)))
      } else {
        cur_qual <- q_line
        cur_val  <- "TRUE"
      }
    } else if (grepl("^                     [^\\/]", line) && !is.null(cur_qual)) {
      # Continuation of qualifier value
      cur_val <- paste0(cur_val, " ", trimws(line))
    }

    if (grepl("^//", line)) flush_feature()
  }

  flush_feature()
  records
}

message("Parsing GenBank file …")
gbff_records <- parse_bakta_gbff(GBFF_FILE)

# ── 3. BUILD FEATURE TABLE ──────────────────────────────────

contig_sizes <- sapply(gbff_records, `[[`, "length") %>%
  sort(decreasing = TRUE)

top_contigs <- names(contig_sizes)[seq_len(min(TOP_N_CONTIGS, length(contig_sizes)))]

extract_features <- function(records, contigs) {
  rows <- lapply(contigs, function(ctg) {
    feats <- records[[ctg]]$features
    if (length(feats) == 0) return(NULL)
    do.call(rbind, lapply(feats, function(f) {
      data.frame(
        contig   = ctg,
        type     = f$type,
        start    = as.integer(f$start),
        end      = as.integer(f$end),
        strand   = as.integer(f$strand),
        gene     = ifelse(is.null(f$gene),      NA, f$gene),
        product  = ifelse(is.null(f$product),   NA, f$product),
        locus    = ifelse(is.null(f$locus_tag), NA, f$locus_tag),
        stringsAsFactors = FALSE
      )
    }))
  })
  bind_rows(rows)
}

feats_all <- extract_features(gbff_records, top_contigs)

# Keep only CDS, tRNA, rRNA, ncRNA features
genes <- feats_all %>%
  filter(type %in% c("CDS", "tRNA", "rRNA", "ncRNA", "tmRNA", "regulatory")) %>%
  mutate(label = coalesce(gene, locus))

# ── 4. ASSIGN FUNCTIONAL CATEGORIES ─────────────────────────

assign_category <- function(type, product) {
  p <- tolower(coalesce(product, ""))
  t <- tolower(type)

  case_when(
    t %in% c("trna", "rrna", "tmrna", "ncrna") ~ "RNA gene",
    str_detect(p, "transposase|insertion element|is\\d|mobile")         ~ "Mobile element",
    str_detect(p, "resist|efflux|beta-lactam|carbapenem|aminoglycoside|tetracycline|quinolone|macrolide") ~ "AMR",
    str_detect(p, "virulence|toxin|hemolysin|fimbriae|fimbria|adhesin|invasion|flagell|type iii|type vi|siderophore|iron") ~ "Virulence",
    str_detect(p, "ribosomal|translation|trna ligase|aminoacyl")        ~ "Translation",
    str_detect(p, "dna polymerase|dna repair|recombinase|helicase|gyrase|topoisomerase|dnaa|dnab|recomb") ~ "DNA replication/repair",
    str_detect(p, "rna polymerase|transcription|sigma factor|regulator|two-component|repressor|activator") ~ "Transcription/regulation",
    str_detect(p, "transport|permease|abc transporter|efflux|outer membrane protein|porin|channel") ~ "Transport",
    str_detect(p, "metabol|synthase|synthetase|dehydrogenase|reductase|kinase|isomerase|transferase|oxidase|phosphatase") ~ "Metabolism",
    str_detect(p, "hypothetical|uncharacterized|putative|unknown|predicted")  ~ "Hypothetical",
    TRUE                                                                  ~ "Other"
  )
}

genes <- genes %>%
  mutate(category = assign_category(type, product))

# ── 5. FUNCTIONAL CATEGORY COLOURS ──────────────────────────

cat_colours <- c(
  "RNA gene"              = "#E41A1C",
  "Mobile element"        = "#FF7F00",
  "AMR"                   = "#A65628",
  "Virulence"             = "#984EA3",
  "Translation"           = "#4DAF4A",
  "DNA replication/repair"= "#377EB8",
  "Transcription/regulation" = "#17BECF",
  "Transport"             = "#F781BF",
  "Metabolism"            = "#8DA0CB",
  "Hypothetical"          = "#B3B3B3",
  "Other"                 = "#666666"
)

genes$colour <- cat_colours[genes$category]
genes$colour[is.na(genes$colour)] <- "#666666"

# ── 6. GC CONTENT & GC SKEW ─────────────────────────────────

compute_gc_tracks <- function(fasta_path, contigs, window = 5000, step = 1000) {
  seqs <- readDNAStringSet(fasta_path)
  names(seqs) <- sub(" .*", "", names(seqs))

  result <- lapply(contigs, function(ctg) {
    if (!ctg %in% names(seqs)) return(NULL)
    s   <- seqs[[ctg]]
    len <- length(s)
    if (len < window) return(NULL)   # contig shorter than window
    starts <- seq(1, len - window + 1, by = step)

    gc_val  <- numeric(length(starts))
    skew    <- numeric(length(starts))
    midpts  <- numeric(length(starts))

    for (i in seq_along(starts)) {
      sub_s   <- subseq(s, starts[i], min(starts[i] + window - 1, len))
      freq    <- letterFrequency(sub_s, c("G", "C", "A", "T"))
      g <- freq["G"]; c <- freq["C"]
      a <- freq["A"]; t <- freq["T"]
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
  error = function(e) {
    message("GC track skipped: ", e$message)
    NULL
  }
)

# ── 7. BUILD GENOME LAYOUT ──────────────────────────────────

# Create a cumulative offset for each contig so they tile around the circle.
# Small spacer gaps separate contigs visually.
SPACER_BP <- 0   # no spacers — keeps Mb ring consistent with genome total

layout <- data.frame(
  contig   = top_contigs,
  len      = as.integer(contig_sizes[top_contigs]),
  stringsAsFactors = FALSE
) %>%
  mutate(
    offset = cumsum(c(0, head(len + SPACER_BP, -1))),
    total  = offset + len
  )

total_genome <- max(layout$total) + SPACER_BP

coords_for <- function(df, layout_df) {
  left_join(df, layout_df %>% select(contig, offset), by = "contig") %>%
    mutate(
      abs_start = start + offset,
      abs_end   = end   + offset
    )
}

genes_abs <- coords_for(genes, layout)

gc_abs <- if (!is.null(gc_data)) {
  left_join(gc_data, layout %>% select(contig, offset), by = "contig") %>%
    mutate(abs_pos = pos + offset)
} else NULL

# ── 8. LOAD & MAP AMR / VIRULENCE GENES ─────────────────────

# --- AMR ---
arg_raw <- read.csv(ARG_FILE, stringsAsFactors = FALSE)

amr_class_col <- c(
  "Beta-Lactamase"            = "#D62728",
  "Efflux pump"               = "#1F77B4",
  "Resistance regulator"      = "#FF7F0E",
  "Aminoglycoside resistance" = "#9467BD",
  "Peptide resistance"        = "#8C564B",
  "Other AMR"                 = "#7F7F7F"
)

assign_amr_class <- function(gene) {
  g <- tolower(gene)
  dplyr::case_when(
    str_detect(g, "bla|amp|ctx")  ~ "Beta-Lactamase",
    str_detect(g, "acrb|acra|acrd|acrf|mdtc|mdtb|mdta|mdtf|mdth|mdtm|emra|emrb|msba|tolc|baer|baes") ~ "Efflux pump",
    str_detect(g, "mara|hns|crp|emrr|kdpe|cpxa|baca|yoji|pmrf") ~ "Resistance regulator",
    str_detect(g, "aph|aac|ant|aad|rmtb|rmtc|arm")              ~ "Aminoglycoside resistance",
    str_detect(g, "mcr|pmr|arn|bac")                             ~ "Peptide resistance",
    TRUE ~ "Other AMR"
  )
}

amr_abs <- arg_raw %>%
  transmute(
    contig = SEQUENCE,
    start  = as.integer(START),
    end    = as.integer(END),
    strand = ifelse(STRAND == "+", 1L, -1L),
    gene   = str_replace(GENE, "Escherichia_coli_", "") %>% str_replace_all("_", "")
  ) %>%
  group_by(contig, gene) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  mutate(
    amr_class = assign_amr_class(gene),
    colour    = amr_class_col[amr_class]
  ) %>%
  left_join(layout %>% select(contig, offset), by = "contig") %>%
  filter(!is.na(offset)) %>%
  mutate(mid_abs = (start + end) / 2 + offset)

# --- Virulence ---
vir_raw <- read_excel(VIR_FILE, skip = 2,
  col_names = c("DATABASE","MAG","FILE","contig","START","END","STRAND",
                "gene","COVERAGE","COVERAGE_MAP","GAPS","pct_cov","pct_id",
                "DATABASE2","ACCESSION","product","VIR_CAT"))

vir_class_col <- c(
  "Iron acquisition"   = "#2CA02C",   # green
  "Curli / Biofilm"    = "#BCBD22",   # yellow-green
  "Flagella"           = "#17BECF",   # sky-teal
  "T2SS"               = "#AEC7E8",   # pale blue
  "T6SS"               = "#00838F",   # dark teal  (was #9467BD — same as Aminoglycoside resistance)
  "O-antigen"          = "#E377C2",   # light pink
  "Invasin"            = "#880E4F",   # deep magenta (was #D62728 — same as Beta-Lactamase)
  "Other virulence"    = "#AAAAAA"    # light grey  (was #7F7F7F — same as Other AMR)
)

assign_vir_class <- function(gene) {
  g <- tolower(gene)
  dplyr::case_when(
    str_detect(g, "^ent|^fep|^fes")           ~ "Iron acquisition",
    str_detect(g, "^csg")                      ~ "Curli / Biofilm",
    str_detect(g, "^fli|^flg|^flh|^che|^mot") ~ "Flagella",
    str_detect(g, "^gsp|^pppd")               ~ "T2SS",
    str_detect(g, "hcp|clpv|^aec|vgr")        ~ "T6SS",
    str_detect(g, "^gtr")                      ~ "O-antigen",
    str_detect(g, "ibe|nada|nadb")             ~ "Invasin",
    TRUE                                        ~ "Other virulence"
  )
}

# Curated label list — one representative per cluster + key standalone genes
vir_label_genes <- c(
  "fepA",   # enterobactin uptake — iron acquisition representative
  "csgD",   # curli master regulator — biofilm representative
  # flhD removed — co-localises with msbA in layout; flagella still shown in legend/ring colours
  "gspD",   # T2SS secretin — secretion representative
  "hcp",    # T6SS tube protein — T6SS representative
  "gtrA",   # O-antigen modification
  "ibeB"    # invasion of brain endothelium
)

vir_abs <- vir_raw %>%
  filter(!is.na(gene), !is.na(contig)) %>%
  mutate(start = as.numeric(START), end = as.numeric(END),
         strand = ifelse(STRAND == "+", 1L, -1L)) %>%
  filter(!is.na(start)) %>%
  dplyr::group_by(contig, gene) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  mutate(
    vir_class = assign_vir_class(gene),
    colour    = vir_class_col[vir_class],
    show_label = gene %in% vir_label_genes
  ) %>%
  left_join(layout %>% select(contig, offset), by = "contig") %>%
  filter(!is.na(offset)) %>%
  mutate(mid_abs = (start + end) / 2 + offset)

# ── 9. DRAW CIRCULAR PLOT ──────────────────────────────────

if (!is.null(OUTPUT_PNG)) {
  png(OUTPUT_PNG, width = 14, height = 14, units = "in", res = 600, bg = "white")
} else if (!is.null(OUTPUT_PDF)) {
  pdf(OUTPUT_PDF, width = 14, height = 14, useDingbats = FALSE)
}

# Give extra canvas space so labels outside the outermost ring aren't clipped
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

# Single sector spanning the whole genome
circos.initialize(factors = "genome", xlim = c(0, total_genome))

# ── Label stagger helper ──────────────────────────────────────
stagger_levels <- function(positions, min_gap, n_levels = 5) {
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
    lv[ord[i]]      <- chosen
    last_pos[chosen] <- sorted_pos[i]
  }
  lv
}

# Strictly one label per genomic region — removes mdtC/mdtH (same cluster as baeR)
amr_label_genes <- c(
  "blaCTX-M-1",   # ESBL
  "blaEC",         # AmpC
  "acrB",          # RND efflux (main pump)
  "acrD",          # RND efflux (aminoglycosides)
  "emrB",          # EmrAB efflux
  "tolC",          # outer membrane channel
  "marA",          # global MDR regulator
  "cpxA",          # CpxAR two-component system
  "baeR",          # BaeSR two-component system
  "msbA"           # lipid A / MDR transporter
)
amr_abs <- amr_abs %>% dplyr::filter(gene %in% amr_label_genes)

LABEL_GAP <- total_genome * 0.075
vir_lab   <- vir_abs %>% dplyr::filter(show_label)

# Build per-strand label tables with abs coords; stagger independently
make_labs <- function(amr, vir, sv) {
  dplyr::bind_rows(
    amr %>% dplyr::filter(strand == sv) %>%
      dplyr::mutate(lgroup = "AMR", lfont = 4,
                    abs_start = start + offset, abs_end = end + offset),
    vir %>% dplyr::filter(strand == sv) %>%
      dplyr::mutate(lgroup = "VIR", lfont = 3,
                    colour = vir_class_col[vir_class],
                    abs_start = start + offset, abs_end = end + offset)
  ) %>% dplyr::arrange(mid_abs)
}

# ── All labels go OUTWARD — solid = forward strand, dashed = reverse ──
# Merge both strands into one stagger so they compete for levels together
fwd_labs <- make_labs(amr_abs, vir_lab,  1L) %>% dplyr::mutate(lty = 1L)
rev_labs <- make_labs(amr_abs, vir_lab, -1L) %>% dplyr::mutate(lty = 2L)

all_labs <- dplyr::bind_rows(fwd_labs, rev_labs) %>% dplyr::arrange(mid_abs)
if (nrow(all_labs) > 0) {
  all_labs$level <- stagger_levels(all_labs$mid_abs, LABEL_GAP, n_levels = 6)
  # Reverse-strand labels get bumped at least one level deeper so ticks are always longer
  rev_idx <- all_labs$lty == 2L
  all_labs$level[rev_idx] <- pmax(all_labs$level[rev_idx], 3L)  # reverse ticks always long enough to show dashes
}

# y levels outside the forward CDS ring (ylim c(0,1); y > 1 is outside)
y_out <- c(1.20, 1.45, 1.70, 1.95, 2.20, 2.45)   # wider spacing; reverse ≥ level 3 = min 3 dashes

sub_fwd <- genes_abs %>% dplyr::filter(strand ==  1L, type == "CDS")
sub_rev <- genes_abs %>% dplyr::filter(strand == -1L, type == "CDS")

# ── Track 1: Forward CDS — outermost, all labels drawn here ───
circos.track(
  ylim = c(0, 1), bg.col = "#F0F0F0", bg.border = NA, track.height = 0.09,
  panel.fun = function(x, y) {
    # CDS bars
    for (i in seq_len(nrow(sub_fwd)))
      circos.rect(sub_fwd$abs_start[i], 0.05, sub_fwd$abs_end[i], 0.95,
                  col = sub_fwd$colour[i], border = NA)
    # Highlighted AMR/VIR genes — bright bar + white border
    for (i in seq_len(nrow(all_labs[all_labs$lty == 1L, ]))) {
      hl <- all_labs[all_labs$lty == 1L, ][i, ]
      if (!is.na(hl$abs_start))
        circos.rect(hl$abs_start, 0.02, hl$abs_end, 0.98,
                    col = hl$colour, border = "white", lwd = 0.8)
    }
    # Labels: only angle when a close neighbour exists, otherwise tick is straight
    ANGLE_STEP <- total_genome * 0.013   # larger step = more separation when fanning
    NEARBY     <- total_genome * 0.04    # catch labels within ~4% of genome (~180 kb)

    if (nrow(all_labs) > 0) {
      for (i in seq_len(nrow(all_labs))) {
        pos  <- all_labs$mid_abs[i]
        col  <- all_labs$colour[i]
        lv   <- all_labs$level[i]
        ylbl <- y_out[lv]
        lty  <- all_labs$lty[i]
        circos.segments(pos, 1.0, pos, ylbl, col = col, lwd = 1.2, lty = lty)
        circos.text(pos, ylbl, all_labs$gene[i],
                    facing = "clockwise", niceFacing = TRUE,
                    cex = 1.0, col = col, font = all_labs$lfont[i],
                    adj = c(0, 0.5))
      }
    }
  }
)

# ── Track 2: Reverse CDS ──────────────────────────────────────
circos.track(
  ylim = c(0, 1), bg.col = "#E8E8E8", bg.border = NA, track.height = 0.09,
  panel.fun = function(x, y) {
    for (i in seq_len(nrow(sub_rev)))
      circos.rect(sub_rev$abs_start[i], 0.05, sub_rev$abs_end[i], 0.95,
                  col = sub_rev$colour[i], border = NA)
    # Highlight reverse-strand AMR/VIR genes on this ring too
    rl <- all_labs[all_labs$lty == 2L, ]
    for (i in seq_len(nrow(rl))) {
      if (!is.na(rl$abs_start[i]))
        circos.rect(rl$abs_start[i], 0.02, rl$abs_end[i], 0.98,
                    col = rl$colour[i], border = "white", lwd = 0.8)
    }
  }
)

# ── Track 3: GC skew ─────────────────────────────────────────
if (!is.null(gc_abs)) {
  skew_lim  <- max(abs(gc_abs$skew), na.rm = TRUE)
  skew_cols <- ifelse(gc_abs$skew >= 0, "#d7191c", "#1a9641")
  circos.track(
    ylim = c(-skew_lim, skew_lim), bg.col = "#EBEBEB", bg.border = NA,
    panel.fun = function(x, y) {
      pos  <- gc_abs$abs_pos
      skew <- gc_abs$skew
      for (i in seq_along(pos))
        circos.lines(c(pos[i], pos[i]), c(0, skew[i]), col = skew_cols[i], lwd = 0.5)
      circos.lines(CELL_META$xlim, c(0, 0), col = "grey50", lwd = 0.5, lty = 2)
    }
  )
}

# ── Track 4: RNA genes ────────────────────────────────────────
rna_genes <- genes_abs %>% dplyr::filter(type %in% c("tRNA", "rRNA", "tmRNA", "ncRNA"))
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

# ── Spacer between RNA ring and Mb ring ──────────────────────
circos.track(ylim = c(0, 1), bg.col = NA, bg.border = NA, track.height = 0.03)

# ── Track 5: Mb scale ring ────────────────────────────────────
circos.track(
  ylim = c(0, 1), track.height = 0.05, bg.border = "grey70", bg.col = "white",
  panel.fun = function(x, y) {
    tick_at <- seq(0, total_genome, by = 500000)
    circos.axis(
      h                 = "bottom",        # labels face inward — away from RNA ring
      major.at          = tick_at,
      labels            = paste0(round(tick_at / 1e6, 1), " Mb"),
      labels.cex        = 0.80,
      labels.facing     = "clockwise",
      labels.niceFacing = TRUE,
      col               = "grey30",
      labels.col        = "grey30",
      lwd               = 0.8,
      minor.ticks       = 4
    )
  }
)

# ── Centre text ──
genome_size_mb <- round(sum(layout$len) / 1e6, 2)
n_cds  <- sum(genes$type == "CDS")
n_trna <- sum(genes$type == "tRNA")
n_rrna <- sum(genes$type == "rRNA")
gc_pct <- if (!is.null(gc_abs)) round(mean(gc_abs$gc) * 100, 1) else "n/a"

# Filled circle background so centre text is always readable
symbols(0, 0, circles = 0.46, inches = FALSE,
        bg = "grey15", fg = NA, add = TRUE)

text(0,  0.16, "Escherichia fergusonii", cex = 1.45, font = 3, col = "white")
text(0,  0.04, "ASC_CAT_CH06",           cex = 1.15, font = 1, col = "grey80")
text(0, -0.10, paste0(n_cds, " CDSs | ", n_trna, " tRNAs | ", n_rrna, " rRNAs"),
               cex = 1.00, font = 1, col = "grey70")
text(0, -0.24, paste0(genome_size_mb, " Mb  |  GC ", gc_pct, "%"),
               cex = 1.05, font = 1, col = "grey70")

# ── Legends: right column = Functional; left column = AMR + Virulence + Strand ──
LEG_CEX   <- 1.05
LEG_INTER <- 0.70   # tighter to prevent overflow with large cex
USR       <- par("usr")   # c(xmin, xmax, ymin, ymax)

cats_present <- unique(genes$category)
legend_cols  <- cat_colours[cats_present]

# Height estimate (canvas units) for a legend block
leg_h <- function(n) (n + 1.8) * LEG_INTER * LEG_CEX * 0.182

# ── RIGHT: Functional category (bottom-right, growing upward) ──
legend(
  x = USR[2], y = USR[3] + 0.05,
  xjust = 1, yjust = 0,
  legend    = names(legend_cols),
  fill      = legend_cols, border = NA,
  bty = "n", cex = LEG_CEX, text.font = 1,
  title = "Functional category", title.font = 2, title.col = "grey15",
  y.intersp = LEG_INTER
)

# ── LEFT: stack from bottom — Strand key, then AMR, then Virulence ──
lx <- USR[1]   # left edge
ly <- USR[3] + 0.05

# Strand key (bottommost on left)
legend(x = lx, y = ly + leg_h(2), xjust = 0, yjust = 0,
  legend = c("(+) Forward strand", "(-) Reverse strand"),
  col = c("grey35","grey35"), lwd = 1.5, lty = c(1,2),
  bty = "n", cex = LEG_CEX * 0.88, text.font = 1,
  title = "Strand (line style)", title.font = 2, title.col = "grey15",
  y.intersp = LEG_INTER)

# AMR class
amr_bot <- ly + leg_h(2) + 0.04
legend(x = lx, y = amr_bot + leg_h(length(amr_class_col)), xjust = 0, yjust = 0,
  legend = names(amr_class_col), col = amr_class_col, lwd = 3,
  bty = "n", cex = LEG_CEX, text.font = 1,
  title = "AMR class", title.font = 2, title.col = "grey15",
  y.intersp = LEG_INTER)

# Virulence category
vir_bot <- amr_bot + leg_h(length(amr_class_col)) + 0.04
legend(x = lx, y = vir_bot + leg_h(length(vir_class_col)), xjust = 0, yjust = 0,
  legend = names(vir_class_col), col = vir_class_col, lwd = 3,
  bty = "n", cex = LEG_CEX, text.font = 1,
  title = "Virulence category", title.font = 2, title.col = "grey15",
  y.intersp = LEG_INTER)

circos.clear()

if (!is.null(OUTPUT_PNG) || !is.null(OUTPUT_PDF)) {
  dev.off()
  message("Plot saved to: ", coalesce(OUTPUT_PNG, OUTPUT_PDF))
}
