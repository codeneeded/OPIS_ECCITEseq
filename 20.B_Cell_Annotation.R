# =============================================================================
# OPIS ECCITEseq — B-cell marker heatmaps (request from Ben)
#
# Purpose: render ADT + RNA marker heatmaps across the annotated B-cell
#          subclusters so Ben can judge whether the naive clusters (and the
#          transitional clusters) can be collapsed for downstream analysis.
#
# Reads:  OPIS_BCELL_Annotated.qs2   (output of the B-cell annotation pipeline)
# Writes: <Bcells>/heatmaps_Ben/  -> PNG heatmaps + CSV tables (nothing deleted)
#
# Conventions carried over from the subclustering pipeline:
#   - Grouping variable = Bcell_Annotation (the manual annotation Ben reviewed)
#   - ADT averaged as plain rowMeans on the DSB "data" slot  (NO expm1)
#   - RNA averaged via AverageExpression on the log1p "data" slot (expm1 internal
#     back-transform is correct there)
#   - Each row min-max scaled 0–1 across clusters for the colour scale; raw
#     (unscaled) values are written to CSV and printed in ADT cells.
#
# ---- RUNTIME ASSUMPTIONS (check these for your session) ---------------------
#   [A1] OPIS_BCELL_Annotated.qs2 lives in `load.path` below.
#   [A2] The grouping column is "Bcell_Annotation". If you'd rather group by the
#        numeric "Subcluster_ID", set group_var <- "Subcluster_ID".
#   [A3] ADT features for CD20/IgD/IgM are physically stored under the gene-style
#        names MS4A1/IGHD/IGHM in the ADT assay (per Ben's rename notes). The
#        resolver below maps Ben's display labels -> whatever is actually present,
#        case/punctuation-insensitively, and reports anything it can't find.
#   [A4] RNA "CD10" and "CD21" are resolved to MME and CR2 respectively.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(viridis)
  library(grid)
})

# ---- Paths -------------------------------------------------------------------
load.path     <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/saved_R_data"
subclust.root <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/subclustering/pre_annotation"
bcell.save.path <- file.path(subclust.root, "Bcells")
out.dir       <- file.path(bcell.save.path, "heatmaps_Ben")
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load annotated B-cell object -------------------------------------------
message("Loading OPIS_BCELL_Annotated.qs2 ...")
OPIS_BCELL <- qs_read(file.path(load.path, "OPIS_BCELL_Annotated.qs2"))

group_var <- "Bcell_Annotation"        # [A2]
stopifnot(group_var %in% colnames(OPIS_BCELL@meta.data))

grp <- droplevels(factor(OPIS_BCELL@meta.data[[group_var]]))
lvls <- levels(grp)
cell_counts <- table(grp)
message("Grouping by '", group_var, "' across: ", paste(lvls, collapse = ", "))

# =============================================================================
# 1.  Marker panels (exactly as Ben listed; label = what shows on the heatmap)
# =============================================================================
adt_markers <- tibble::tribble(
  ~category,       ~label,    ~aliases,
  "Lineage",       "CD19",    "",
  "Lineage",       "CD20",    "MS4A1",
  "Lineage",       "CD22",    "",
  "Lineage",       "CD79B",   "",
  "Lineage",       "HLA-DR",  "HLA.DR,HLADR,HLA-DRA",
  "Naive/Memory",  "IgD",     "IGHD",
  "Naive/Memory",  "IgM",     "IGHM",
  "Naive/Memory",  "CD27",    "",
  "Naive/Memory",  "CD38",    "",
  "Naive/Memory",  "CD24",    "",
  "Naive/Memory",  "CD86",    "",
  "Atypical",      "CXCR5",   "",
  "Atypical",      "LILRB1",  "CD85J",
  "Activation",    "CD40",    "",
  "Activation",    "CD69",    ""
)
adt_cat_levels <- c("Lineage", "Naive/Memory", "Atypical", "Activation")

rna_markers <- tibble::tribble(
  ~category,                ~label,      ~aliases,
  "Lineage",                "MS4A1",     "",
  "Lineage",                "CD79A",     "",
  "Lineage",                "CD79B",     "",
  "Lineage",                "CD74",      "",
  "Lineage",                "BANK1",     "",
  "Lineage",                "CD19",      "",
  "Lineage",                "BLK",       "",
  "Naive B (resting)",      "IGHD",      "",
  "Naive B (resting)",      "IGHM",      "",
  "Naive B (resting)",      "TCL1A",     "",
  "Naive B (resting)",      "IL4R",      "",
  "Naive B (resting)",      "FCER2",     "",
  "Naive B (resting)",      "CCR7",      "",
  "Naive B (resting)",      "CXCR4",     "",
  "Naive B (resting)",      "BACH2",     "",
  "Transitional B",         "CD38",      "",
  "Transitional B",         "CD10",      "MME",
  "Transitional B",         "SOX4",      "",
  "Transitional B",         "VPREB3",    "",
  "Transitional B",         "IGLL5",     "",
  "Transitional B",         "TCL1A",     "",
  "Transitional B",         "CD24",      "",
  "Transitional B",         "IL7R",      "",
  "Transitional B",         "CD21",      "CR2",
  "IEG-High",               "FOS",       "",
  "IEG-High",               "JUN",       "",
  "IEG-High",               "JUNB",      "",
  "IEG-High",               "DUSP1",     "",
  "IEG-High",               "EGR1",      "",
  "IEG-High",               "IER2",      "",
  "IEG-High",               "NR4A1",     "",
  "Intermediate HLA-G+",    "HLA-G",     "HLA.G",
  "Intermediate HLA-G+",    "HLA-DRA",   "HLA.DRA",
  "Intermediate HLA-G+",    "HLA-DPA1",  "HLA.DPA1",
  "Intermediate HLA-G+",    "HLA-DPB1",  "HLA.DPB1",
  "Intermediate HLA-G+",    "HLA-DRB1",  "HLA.DRB1",
  "Intermediate HLA-G+",    "AIM2",      "",
  "Atypical B",             "TBX21",     "",
  "Atypical B",             "ITGAX",     "",
  "Atypical B",             "FCRL5",     "",
  "Atypical B",             "FCRL3",     "",
  "Atypical B",             "ZEB2",      "",
  "Atypical B",             "LILRB1",    "",
  "Atypical B",             "LILRB2",    "",
  "Atypical B",             "CXCR3",     "",
  "Memory B",               "CD27",      "",
  "Memory B",               "COCH",      "",
  "Memory B",               "LTB",       "",
  "Memory B",               "CD40",      "",
  "Memory B",               "CD83",      "",
  "Memory B",               "TNFRSF13B", "",
  "Plasmablast exclusion",  "JCHAIN",    "",
  "Plasmablast exclusion",  "XBP1",      "",
  "Plasmablast exclusion",  "MZB1",      "",
  "Plasmablast exclusion",  "PRDM1",     ""
)
rna_cat_levels <- c("Lineage", "Naive B (resting)", "Transitional B", "IEG-High",
                    "Intermediate HLA-G+", "Atypical B", "Memory B",
                    "Plasmablast exclusion")

# =============================================================================
# 2.  Helpers
# =============================================================================
.norm <- function(x) toupper(gsub("[-_. ]", "", x))

# Resolve a display label (+ aliases) to whatever feature is actually present.
resolve_feature <- function(label, aliases, available) {
  cand <- c(label, if (nzchar(aliases)) strsplit(aliases, ",")[[1]] else character(0))
  cand <- trimws(cand)
  for (c0 in cand) if (c0 %in% available) return(c0)        # exact
  av_norm <- .norm(available)
  for (c0 in cand) {                                        # punctuation/case-insensitive
    hit <- available[av_norm == .norm(c0)]
    if (length(hit)) return(hit[1])
  }
  NA_character_
}

# min-max scale each row to 0–1 across columns (constant rows -> 0)
scale_01 <- function(m) {
  m <- as.matrix(m)
  out <- t(apply(m, 1, function(x) {
    rng <- range(x, na.rm = TRUE)
    if (!is.finite(diff(rng)) || diff(rng) == 0) return(rep(0, length(x)))
    (x - rng[1]) / (rng[2] - rng[1])
  }))
  dimnames(out) <- dimnames(m)
  out
}

# Resolve a marker table against an assay's available features; drop missing,
# de-duplicate repeated features (keep first category), preserve category order.
prepare_panel <- function(markers, cat_levels, available, modality) {
  res <- markers %>%
    mutate(feature = mapply(resolve_feature, label, aliases,
                            MoreArgs = list(available = available)))
  missing <- res %>% filter(is.na(feature))
  if (nrow(missing))
    message("[", modality, "] not found, skipped: ",
            paste(missing$label, collapse = ", "))
  res <- res %>%
    filter(!is.na(feature)) %>%
    mutate(category = factor(category, levels = cat_levels)) %>%
    arrange(category)
  dup <- duplicated(res$feature)
  if (any(dup))
    message("[", modality, "] duplicate feature(s) kept once (first category): ",
            paste(unique(res$feature[dup]), collapse = ", "))
  res[!dup, , drop = FALSE]
}

# Map annotation level -> coarse developmental stage (for the column annotation)
stage_of <- function(x) dplyr::recode(as.character(x),
                                      "Transitional B"               = "Transitional",
                                      "Transitional/Naive (T2-like)" = "Transitional",
                                      "Naive B (resting)"            = "Naive",
                                      "Naive B (IEG-high)"           = "Naive",
                                      "Intermediate B (HLA-G+)"      = "Intermediate",
                                      "Memory B"                     = "Memory",
                                      "Atypical B (ABC/DN2)"         = "Atypical",
                                      .default = "Other")

stage_cols <- c(Transitional = "#4C72B0", Naive = "#55A868",
                Intermediate = "#C44E52", Memory = "#8172B3",
                Atypical = "#CCB974", Other = "#999999")

col_lab   <- paste0(lvls, "\n(n=", as.integer(cell_counts[lvls]), ")")
col_stage <- stage_of(lvls)

# =============================================================================
# 3.  Average-expression matrices
# =============================================================================
# ---- ADT: manual rowMeans on DSB "data" slot (NO expm1) ---------------------
DefaultAssay(OPIS_BCELL) <- "ADT"
adt_avail <- rownames(OPIS_BCELL[["ADT"]])
adt_panel <- prepare_panel(adt_markers, adt_cat_levels, adt_avail, "ADT")

adt_data <- GetAssayData(OPIS_BCELL, assay = "ADT", slot = "data")[adt_panel$feature, , drop = FALSE]
adt_avg  <- sapply(lvls, function(cl) {
  idx <- which(grp == cl)
  rowMeans(adt_data[, idx, drop = FALSE])
})
rownames(adt_avg) <- adt_panel$label          # display labels (CD20, IgD, IgM ...)
colnames(adt_avg) <- lvls
adt_scaled <- scale_01(adt_avg)

# ---- RNA: AverageExpression on log1p "data" slot ----------------------------
DefaultAssay(OPIS_BCELL) <- "RNA"
rna_avail <- rownames(OPIS_BCELL[["RNA"]])
rna_panel <- prepare_panel(rna_markers, rna_cat_levels, rna_avail, "RNA")

rna_avg_raw <- AverageExpression(OPIS_BCELL, assays = "RNA",
                                 features = rna_panel$feature,
                                 group.by = group_var, slot = "data")$RNA
# robustly align AverageExpression columns to our level order
ci <- match(.norm(lvls), .norm(colnames(rna_avg_raw)))
if (anyNA(ci)) ci <- seq_along(lvls)          # fall back to native order if names mangle
rna_avg <- rna_avg_raw[rna_panel$feature, ci, drop = FALSE]
rownames(rna_avg) <- rna_panel$label
colnames(rna_avg) <- lvls
rna_scaled <- scale_01(rna_avg)

# =============================================================================
# 4.  Heatmap renderer
# =============================================================================
draw_panel_heatmap <- function(scaled, raw, row_cat, title, file, show_values) {
  col_fun <- colorRamp2(seq(0, 1, length.out = 9), viridis(9, option = "A"))
  top_anno <- HeatmapAnnotation(
    Stage = col_stage,
    col = list(Stage = stage_cols[unique(col_stage)]),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 9),
    show_legend = TRUE
  )
  cell_fun <- NULL
  if (show_values) cell_fun <- function(j, i, x, y, w, h, fill)
    grid.text(sprintf("%.2f", raw[i, j]), x, y, gp = gpar(fontsize = 6, col = "grey25"))
  
  ht <- Heatmap(
    scaled,
    name              = "Scaled\n(0–1)",
    col               = col_fun,
    cluster_rows      = FALSE,
    row_split         = row_cat,
    row_title_rot     = 0,
    row_title_gp      = gpar(fontsize = 10, fontface = "bold"),
    row_gap           = unit(2, "mm"),
    row_names_gp      = gpar(fontsize = 9),
    cluster_columns   = TRUE,                 # column dendrogram = the merge question
    column_dend_height = unit(12, "mm"),
    column_labels     = col_lab,
    column_names_rot  = 45,
    column_names_gp   = gpar(fontsize = 9),
    top_annotation    = top_anno,
    rect_gp           = gpar(col = "white", lwd = 0.5),
    cell_fun          = cell_fun,
    column_title      = title,
    column_title_gp   = gpar(fontsize = 13, fontface = "bold"),
    heatmap_legend_param = list(at = c(0, 0.5, 1))
  )
  
  h_in <- 2.2 + 0.20 * nrow(scaled)
  w_in <- 4.5 + 0.95 * ncol(scaled)
  png(file, width = w_in, height = h_in, units = "in", res = 400)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
       merge_legend = TRUE)
  dev.off()
  message("Saved: ", file)
}

draw_panel_heatmap(adt_scaled, adt_avg, adt_panel$category,
                   "B cells | ADT panel (row-scaled 0–1; cell = mean DSB)",
                   file.path(out.dir, "Bcell_ADT_marker_heatmap.png"),
                   show_values = TRUE)

draw_panel_heatmap(rna_scaled, rna_avg, rna_panel$category,
                   "B cells | RNA panel (row-scaled 0–1 mean expression)",
                   file.path(out.dir, "Bcell_RNA_marker_heatmap.png"),
                   show_values = FALSE)

# =============================================================================
# 5.  Cluster-similarity heatmap  (directly addresses the merge question)
#     Pearson correlation of the combined scaled marker profile between clusters.
#     High correlation between the two naive clusters / the two transitional
#     clusters -> good candidates to collapse.
# =============================================================================
combined_scaled <- rbind(adt_scaled, rna_scaled)
cl_cor <- cor(combined_scaled, method = "pearson")

cor_fun <- colorRamp2(c(min(cl_cor), (min(cl_cor) + 1) / 2, 1),
                      c("#f7fbff", "#6baed6", "#08306b"))
cor_anno <- HeatmapAnnotation(Stage = col_stage,
                              col = list(Stage = stage_cols[unique(col_stage)]),
                              show_annotation_name = FALSE)
ht_cor <- Heatmap(
  cl_cor, name = "Pearson r", col = cor_fun,
  column_labels = col_lab, row_labels = col_lab,
  column_names_rot = 45, column_names_gp = gpar(fontsize = 9),
  row_names_gp = gpar(fontsize = 9),
  top_annotation = cor_anno,
  rect_gp = gpar(col = "white", lwd = 1),
  cell_fun = function(j, i, x, y, w, h, fill)
    grid.text(sprintf("%.2f", cl_cor[i, j]), x, y, gp = gpar(fontsize = 8)),
  column_title = "B-cell subcluster similarity (combined ADT+RNA marker profile)",
  column_title_gp = gpar(fontsize = 13, fontface = "bold")
)
png(file.path(out.dir, "Bcell_cluster_similarity_heatmap.png"),
    width = 4.5 + 0.95 * ncol(cl_cor), height = 4.0 + 0.55 * ncol(cl_cor),
    units = "in", res = 400)
draw(ht_cor, merge_legend = TRUE)
dev.off()
message("Saved: ", file.path(out.dir, "Bcell_cluster_similarity_heatmap.png"))

# =============================================================================
# 6.  CSV exports (numbers behind the figures) + feature-coverage report
# =============================================================================
write.csv(round(adt_avg, 4),    file.path(out.dir, "AvgExpr_ADT_panel_raw.csv"))
write.csv(round(adt_scaled, 4), file.path(out.dir, "AvgExpr_ADT_panel_scaled.csv"))
write.csv(round(rna_avg, 4),    file.path(out.dir, "AvgExpr_RNA_panel_raw.csv"))
write.csv(round(rna_scaled, 4), file.path(out.dir, "AvgExpr_RNA_panel_scaled.csv"))
write.csv(round(cl_cor, 4),     file.path(out.dir, "Cluster_similarity_correlation.csv"))

coverage <- bind_rows(
  adt_markers %>% mutate(modality = "ADT",
                         feature = mapply(resolve_feature, label, aliases,
                                          MoreArgs = list(available = adt_avail))),
  rna_markers %>% mutate(modality = "RNA",
                         feature = mapply(resolve_feature, label, aliases,
                                          MoreArgs = list(available = rna_avail)))
) %>%
  transmute(modality, category, requested = label, resolved_to = feature,
            found = !is.na(feature))
write.csv(coverage, file.path(out.dir, "Marker_coverage_report.csv"), row.names = FALSE)

message("\nFeature coverage:")
print(coverage %>% count(modality, found))
message("\nAll outputs written to: ", out.dir)
message("Done.")

