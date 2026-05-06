###############################################################################
# 18 — OPIS Manuscript-Figure TCR / TCRex analyses
# ============================================================================
# Builds on the objects produced by:
#   16_vdj.R                 -> OPIS_CD8_TCR_Combined.qs2
#                                (OPIS_CD8 with cloneSize, OPIS_CD8_Trex ED0/1/2)
#   17_Module_Scoring.R       -> module-score helpers (replicated below)
#
# All outputs go to:
#   <project.root>/Manuscript_Figures/TCR/
#
# Subfolders created automatically (one per analysis block):
#   - Cluster_UMAPs/              individual UMAPs highlighting TEMRA, Innate, PRF1+ Tem
#   - Clonal_Distribution/        clones distributed across (current) clusters
#   - DGE_PRF1Tem/                clones vs non-clones in PRF1+ Tem
#   - Modules_HIV_CMV_PRF1Tem/    Exhaustion / Cytotoxic / IFN modules, HIV vs CMV
#   - DGE_HIV_vs_CMV_PRF1Tem/     HIV vs CMV DGE within PRF1+ Tem
#   - Specificity_vs_CloneSize/   xy plot + clone-size-bin bar chart
#   - HIV_OUD/                    DGE & modules in CD8 HIV-specific by OUD
#   - Shared_Clones/              clones shared across participants
#   - Stats/                      3 stats CSVs requested
#
# IMPORTANT — verify these labels against your Seurat objects before running:
#   prf1_tem_label <- "CD8+ Tem (PRF1+)"   # search-replace if different
#   temra_label    <- "CD8+ TEMRA"
#   innate_label   <- "Innate-like T"
###############################################################################

# ---- Libraries -------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(scRepertoire)
  library(Trex)
  library(tidyverse)
  library(cowplot)
  library(scales)
  library(scCustomize)
  library(RColorBrewer)
  library(Polychrome)
  library(viridis)
  library(patchwork)
  library(qs2)
})

# ---- Paths -----------------------------------------------------------------
# load.path = where script 16 saved the qs2 packages (input)
# project.root = top-level OPIS folder
# manu.root = Manuscript_Figures/TCR sits at the project root, NOT nested
#             inside VDJ/TCR/, to avoid the awkward .../TCR/Manuscript_Figures/TCR/ path.
load.path    <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/saved_R_data"
project.root <- "~/Documents/OPIS_ECCITEseq"
manu.root    <- file.path(project.root, "Manuscript_Figures", "TCR")

dir.create(manu.root, recursive = TRUE, showWarnings = FALSE)

sub_dirs <- c(
  "Cluster_UMAPs",
  "Clonal_Distribution",
  "DGE_PRF1Tem",
  "Modules_HIV_CMV_PRF1Tem",
  "DGE_HIV_vs_CMV_PRF1Tem",
  "Specificity_vs_CloneSize",
  "HIV_OUD",
  "Shared_Clones",
  "Stats"
)
for (sd in sub_dirs) dir.create(file.path(manu.root, sd),
                                recursive = TRUE, showWarnings = FALSE)

cluster_umap_dir   <- file.path(manu.root, "Cluster_UMAPs")
clonal_dist_dir    <- file.path(manu.root, "Clonal_Distribution")
dge_prf1_dir       <- file.path(manu.root, "DGE_PRF1Tem")
mod_hivcmv_dir     <- file.path(manu.root, "Modules_HIV_CMV_PRF1Tem")
dge_hivcmv_dir     <- file.path(manu.root, "DGE_HIV_vs_CMV_PRF1Tem")
spec_size_dir      <- file.path(manu.root, "Specificity_vs_CloneSize")
hiv_oud_dir        <- file.path(manu.root, "HIV_OUD")
shared_clones_dir  <- file.path(manu.root, "Shared_Clones")
stats_dir          <- file.path(manu.root, "Stats")

###############################################################################
# 1. CONFIG --------------------------------------------------------------
###############################################################################

# Metadata columns
oud_col      <- "OUD_status"
oud_pos      <- "OUD+"
oud_neg      <- "OUD-"
oud_palette  <- c(`OUD-` = "#3B7FB8", `OUD+` = "#D95F02")

celltype_col <- "celltype_annotation"
red_name     <- "wnn.umap"

# Cluster labels — VERIFY against unique(OPIS_CD8@meta.data[[celltype_col]])
prf1_tem_label <- "CD8+ Tem (PRF1+)"
temra_label    <- "CD8+ TEMRA"
innate_label   <- "Innate-like T"

# Highlight colors for the individual cluster UMAPs
cluster_highlight_palette <- c(
  `CD8+ TEMRA`        = "#C4463A",
  `Innate-like T`     = "#7570B3",
  `CD8+ Tem (PRF1+)`  = "#D95F02"
)

# Specificity colors (used in several plots)
specificity_palette <- c(
  HIV   = "#C4463A",
  CMV   = "#3B7FB8",
  Other = "grey75"
)

# Clone-size bin definitions (must match what was passed to combineExpression
# in script 16 so the cloneSize meta column lines up).
clone_size_bins <- c(Single = 1, Small = 5, Medium = 20,
                     Large = 100, Hyperexpanded = 500)

# Edit-distance choice for specificity calls. ED1 used in script 16.
trex_ed_to_use <- "ed1"

###############################################################################
# 2. LOAD --------------------------------------------------------------
###############################################################################

cat("\nLoading saved CD8 TCR package...\n")
cd8_pkg <- qs2::qs_read(file.path(load.path, "OPIS_CD8_TCR_Combined.qs2"))
OPIS_CD8 <- cd8_pkg$OPIS_CD8

trex_obj_cd8 <- switch(trex_ed_to_use,
                       ed0 = cd8_pkg$OPIS_CD8_ED0,
                       ed1 = cd8_pkg$OPIS_CD8_ED1,
                       ed2 = cd8_pkg$OPIS_CD8_ED2)
stopifnot(!is.null(trex_obj_cd8))

# CD4 is loaded for the clonal-distribution block only (Block 5).
# Other blocks operate exclusively on the CD8 object / Trex CD8.
cat("\nLoading saved CD4 TCR package (clonal distribution only)...\n")
cd4_pkg  <- qs2::qs_read(file.path(load.path, "OPIS_CD4_TCR_Combined.qs2"))
OPIS_CD4 <- cd4_pkg$OPIS_CD4

# Sanity print
cat("CD8 cells:", ncol(OPIS_CD8), "\n")
cat("CD8 clusters present:\n")
print(table(OPIS_CD8@meta.data[[celltype_col]], useNA = "ifany"))
cat("\nCD4 cells:", ncol(OPIS_CD4), "\n")
cat("CD4 clusters present:\n")
print(table(OPIS_CD4@meta.data[[celltype_col]], useNA = "ifany"))
cat("\nReductions on OPIS_CD8: ", paste(Reductions(OPIS_CD8), collapse = ", "), "\n")
if (!red_name %in% Reductions(OPIS_CD8))
  warning("Reduction '", red_name, "' not found on OPIS_CD8 — UMAPs will fail.")

# Cluster sanity
miss <- setdiff(c(prf1_tem_label, temra_label, innate_label),
                unique(OPIS_CD8@meta.data[[celltype_col]]))
if (length(miss))
  warning("Cluster label(s) not present in OPIS_CD8: ",
          paste(miss, collapse = " | "),
          "\n  -> Edit the labels at the top of the script to match `levels()`.")

###############################################################################
# 3. HELPERS ------------------------------------------------------------
###############################################################################

# --- UMAP arrows in bottom-left corner --------------------------------------
# Strips axis line/text/ticks/title from a DimPlot and overlays two small
# arrows (UMAP 1 →, UMAP 2 ↑) plus labels in the bottom-left of the canvas.
# Coordinates are in npc (0–1 across the whole ggdraw canvas) so the arrows
# sit consistently regardless of UMAP scale or axis units.
# Returns a ggdraw object — still saveable with ggsave().
add_umap_arrows <- function(p,
                            arrow_len   = 0.13,   # length of each arrow in npc
                            x_origin    = 0.04,   # arrow base, x (npc)
                            y_origin    = 0.05,   # arrow base, y (npc)
                            label_size  = 9,
                            line_width  = 0.6,
                            line_color  = "black") {
  p <- p + theme(
    axis.line   = element_blank(),
    axis.text   = element_blank(),
    axis.ticks  = element_blank(),
    axis.title  = element_blank()
  )
  cowplot::ggdraw(p) +
    # Horizontal arrow → UMAP 1
    cowplot::draw_line(
      x = c(x_origin, x_origin + arrow_len),
      y = c(y_origin, y_origin),
      arrow = grid::arrow(length = grid::unit(2, "mm"),
                          type = "closed", angle = 25),
      color = line_color, size = line_width) +
    # Vertical arrow ↑ UMAP 2
    cowplot::draw_line(
      x = c(x_origin, x_origin),
      y = c(y_origin, y_origin + arrow_len),
      arrow = grid::arrow(length = grid::unit(2, "mm"),
                          type = "closed", angle = 25),
      color = line_color, size = line_width) +
    cowplot::draw_text("UMAP 1",
                       x = x_origin + arrow_len / 2,
                       y = y_origin - 0.025,
                       hjust = 0.5, vjust = 1, size = label_size) +
    cowplot::draw_text("UMAP 2",
                       x = x_origin - 0.025,
                       y = y_origin + arrow_len / 2,
                       hjust = 0.5, vjust = 0,
                       size = label_size, angle = 90)
}

# Mark cells as having a TCR / being in an expanded clone (script 16 helper)
mark_expanded <- function(obj) {
  cs <- as.character(obj$cloneSize)
  obj$has_TCR     <- !is.na(cs)
  obj$is_expanded <- obj$has_TCR & !grepl("^Single", cs)
  obj
}

# Classify each cell as HIV / CMV / Other (script 16 helper, slight tweak).
# "Other" here = has TCR but no HIV/CMV epitope hit, OR no TCR.
classify_specificity <- function(obj) {
  meta <- obj@meta.data
  spec <- rep("Other", nrow(meta))
  if ("TRB_Epitope.species" %in% colnames(meta)) {
    sp <- as.character(meta[["TRB_Epitope.species"]])
    is_hiv <- grepl("HIV|Human immunodeficiency virus", sp, ignore.case = TRUE)
    is_cmv <- grepl("CMV|Cytomegalovirus",              sp, ignore.case = TRUE)
    spec[is_hiv]            <- "HIV"
    spec[is_cmv & !is_hiv]  <- "CMV"
  } else {
    warning("TRB_Epitope.species not found; all cells default to 'Other'.")
  }
  obj$specificity <- factor(spec, levels = c("HIV","CMV","Other"))
  obj
}

OPIS_CD8     <- mark_expanded(OPIS_CD8)
OPIS_CD4     <- mark_expanded(OPIS_CD4)
trex_obj_cd8 <- mark_expanded(trex_obj_cd8)
trex_obj_cd8 <- classify_specificity(trex_obj_cd8)

# Bring specificity back onto OPIS_CD8 by cell barcode (so we can do all
# downstream subsetting on a single object)
spec_lookup <- setNames(as.character(trex_obj_cd8$specificity),
                        rownames(trex_obj_cd8@meta.data))
OPIS_CD8$specificity <- factor(
  unname(spec_lookup[Cells(OPIS_CD8)]),
  levels = c("HIV","CMV","Other"))
# Cells in OPIS_CD8 not in trex object (shouldn't happen but be safe) get NA -> "Other"
OPIS_CD8$specificity[is.na(OPIS_CD8$specificity)] <- "Other"
OPIS_CD8$specificity <- factor(OPIS_CD8$specificity,
                               levels = c("HIV","CMV","Other"))

cat("\nSpecificity table on OPIS_CD8:\n")
print(table(OPIS_CD8$specificity, useNA = "ifany"))

# Cells in a given cluster
cells_in_cluster <- function(obj, cluster_label) {
  rownames(obj@meta.data)[
    !is.na(obj@meta.data[[celltype_col]]) &
      obj@meta.data[[celltype_col]] == cluster_label
  ]
}

# Gini coefficient (standard formula, no external pkg)
gini <- function(x) {
  x <- as.numeric(x)
  x <- x[!is.na(x) & x >= 0]
  n <- length(x)
  if (n < 2 || sum(x) == 0) return(NA_real_)
  x <- sort(x)
  G <- (2 * sum(seq_len(n) * x) / (n * sum(x))) - (n + 1) / n
  G
}

###############################################################################
# 4. INDIVIDUAL CLUSTER UMAPs ---------------------------------------
#    (TEMRA, Innate-like, PRF1+ Tem within the CD8+ UMAP)
###############################################################################

cat("\n=== Block 4: cluster-highlight UMAPs ===\n")

highlight_one_cluster <- function(obj, cluster_label, color, file_tag) {
  highlight_cells <- cells_in_cluster(obj, cluster_label)
  if (length(highlight_cells) == 0) {
    message("No cells in '", cluster_label, "' — skipping.")
    return(invisible())
  }
  p <- DimPlot(obj, reduction = red_name,
               cells.highlight = list(`Cluster` = highlight_cells),
               cols.highlight  = color, cols = "grey85",
               sizes.highlight = 0.9, pt.size = 0.4) +
    ggtitle(paste0(cluster_label, "  (n = ", length(highlight_cells), ")")) +
    theme_cowplot(font_size = 14) +
    theme(plot.title      = element_text(face = "bold", hjust = 0.5),
          legend.position = "none")
  p <- add_umap_arrows(p)
  ggsave(file.path(cluster_umap_dir, paste0("CD8_UMAP_", file_tag, ".png")),
         p, width = 7, height = 6, dpi = 300, bg = "white")
}

highlight_one_cluster(OPIS_CD8, temra_label,
                      cluster_highlight_palette[temra_label],    "TEMRA")
highlight_one_cluster(OPIS_CD8, innate_label,
                      cluster_highlight_palette[innate_label],   "InnateLike")
highlight_one_cluster(OPIS_CD8, prf1_tem_label,
                      cluster_highlight_palette[prf1_tem_label], "PRF1pos_Tem")

# Combined panel for convenience
combo_cells <- list(
  `CD8+ TEMRA`       = cells_in_cluster(OPIS_CD8, temra_label),
  `Innate-like T`    = cells_in_cluster(OPIS_CD8, innate_label),
  `CD8+ Tem (PRF1+)` = cells_in_cluster(OPIS_CD8, prf1_tem_label)
)
combo_cells <- combo_cells[lengths(combo_cells) > 0]
if (length(combo_cells) > 0) {
  p <- DimPlot(OPIS_CD8, reduction = red_name,
               cells.highlight = combo_cells,
               cols.highlight  = cluster_highlight_palette[names(combo_cells)],
               cols = "grey85", sizes.highlight = 0.7, pt.size = 0.3) +
    ggtitle("CD8: TEMRA / Innate-like / PRF1+ Tem") +
    theme_cowplot(font_size = 14) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  p <- add_umap_arrows(p)
  ggsave(file.path(cluster_umap_dir, "CD8_UMAP_HighlightAllThree.png"),
         p, width = 8, height = 6, dpi = 300, bg = "white")
}

###############################################################################
# 5. CLONAL DISTRIBUTION ACROSS CLUSTERS ------------------------------
#    Updated bar charts using current celltype_annotation labels.
#    Run for both CD8 and CD4. All outputs share the same Clonal_Distribution
#    folder, distinguished by CD8_ / CD4_ prefix.
###############################################################################

cat("\n=== Block 5: clonal distribution across clusters ===\n")

run_clonal_distribution <- function(obj, tag) {
  cat("  -- ", tag, "\n", sep = "")
  
  # (a) scRepertoire's clonalOccupy with celltype_annotation
  #     (counts per cluster, stacked by cloneSize bin)
  p_count <- clonalOccupy(obj, x.axis = celltype_col, label = FALSE) +
    theme_cowplot(font_size = 14) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"),
          plot.title  = element_text(face = "bold", hjust = 0.5)) +
    ggtitle(paste0(tag,
                   ": clonal distribution across clusters (cell counts)"))
  ggsave(file.path(clonal_dist_dir,
                   paste0(tag, "_ClonalOccupancy_byCluster_count.png")),
         p_count, width = 13, height = 7, dpi = 300, bg = "white")
  
  p_prop <- clonalOccupy(obj, x.axis = celltype_col,
                         proportion = TRUE, label = FALSE) +
    theme_cowplot(font_size = 14) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"),
          plot.title  = element_text(face = "bold", hjust = 0.5)) +
    ggtitle(paste0(tag,
                   ": clonal distribution across clusters (proportions)"))
  ggsave(file.path(clonal_dist_dir,
                   paste0(tag, "_ClonalOccupancy_byCluster_proportion.png")),
         p_prop, width = 13, height = 7, dpi = 300, bg = "white")
  
  # Export underlying table
  occ_tab <- clonalOccupy(obj, x.axis = celltype_col, exportTable = TRUE)
  write.csv(occ_tab,
            file.path(clonal_dist_dir,
                      paste0(tag, "_ClonalOccupancy_byCluster_table.csv")),
            row.names = FALSE)
  
  # (b) Custom: % of clonal vs non-clonal cells per cluster
  clone_status_df <- obj@meta.data %>%
    dplyr::filter(!is.na(.data[[celltype_col]]), has_TCR) %>%
    dplyr::mutate(clone_status = ifelse(is_expanded, "Clonal", "Non-clonal")) %>%
    dplyr::count(.data[[celltype_col]], clone_status) %>%
    dplyr::group_by(.data[[celltype_col]]) %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    dplyr::ungroup()
  
  write.csv(clone_status_df,
            file.path(clonal_dist_dir,
                      paste0(tag, "_ClonalStatus_byCluster.csv")),
            row.names = FALSE)
  
  p_status <- ggplot(clone_status_df,
                     aes(x = .data[[celltype_col]], y = prop, fill = clone_status)) +
    geom_col(width = 0.75, color = "grey15", linewidth = 0.3) +
    scale_fill_manual(values = c(Clonal = "#C4463A", `Non-clonal` = "#9DB7C9"),
                      name = NULL) +
    scale_y_continuous(labels = scales::percent_format(),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = "Proportion of cells (with TCR)",
         title = paste0(tag, ": clonal vs non-clonal cells per cluster")) +
    theme_cowplot(font_size = 14) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"),
          plot.title  = element_text(face = "bold", hjust = 0.5))
  ggsave(file.path(clonal_dist_dir,
                   paste0(tag, "_ClonalStatus_byCluster_proportion.png")),
         p_status, width = 11, height = 6, dpi = 300, bg = "white")
}

run_clonal_distribution(OPIS_CD8, "CD8")
run_clonal_distribution(OPIS_CD4, "CD4")

###############################################################################
# 6. DGE: CLONES vs NON-CLONES IN CD8+ PRF1+ Tem ----------------------
###############################################################################

cat("\n=== Block 6: DGE clones vs non-clones in PRF1+ Tem ===\n")

dge_clones_in_cluster <- function(obj, cluster_label, file_tag, out_dir) {
  cells <- cells_in_cluster(obj, cluster_label)
  if (length(cells) < 20) {
    message("Cluster '", cluster_label, "' has only ", length(cells),
            " cells; skipping."); return(invisible())
  }
  obj_sub <- subset(obj, cells = cells)
  obj_sub <- subset(obj_sub,
                    cells = rownames(obj_sub@meta.data)[obj_sub$has_TCR])
  
  n_clo <- sum(obj_sub$is_expanded, na.rm = TRUE)
  n_non <- sum(!obj_sub$is_expanded, na.rm = TRUE)
  cat("  Cluster ", cluster_label, ": clonal=", n_clo,
      " non-clonal=", n_non, "\n", sep = "")
  if (n_clo < 5 || n_non < 5) {
    message("Too few cells; skipping DGE."); return(invisible())
  }
  
  obj_sub$clone_status <- ifelse(obj_sub$is_expanded, "Clonal", "NonClonal")
  Idents(obj_sub) <- obj_sub$clone_status
  DefaultAssay(obj_sub) <- "RNA"
  m <- FindMarkers(obj_sub, ident.1 = "Clonal", ident.2 = "NonClonal",
                   logfc.threshold = 0.1, min.pct = 0.1)
  m$gene <- rownames(m)
  m <- m %>% dplyr::arrange(p_val_adj, dplyr::desc(abs(avg_log2FC))) %>%
    dplyr::select(gene, dplyr::everything())
  write.csv(m, file.path(out_dir,
                         paste0("CD8_", file_tag, "_Clonal_vs_NonClonal_DGE.csv")),
            row.names = FALSE)
  
  # Quick volcano for the figure folder (purely illustrative; user can re-style)
  m_vol <- m %>%
    dplyr::mutate(neg_log10_padj = -log10(p_val_adj + 1e-300),
                  signif = p_val_adj < 0.05 & abs(avg_log2FC) > 0.25,
                  direction = dplyr::case_when(
                    signif & avg_log2FC > 0 ~ "Up in Clonal",
                    signif & avg_log2FC < 0 ~ "Down in Clonal",
                    TRUE                    ~ "ns"))
  p <- ggplot(m_vol, aes(x = avg_log2FC, y = neg_log10_padj, color = direction)) +
    geom_point(alpha = 0.7, size = 1.4) +
    scale_color_manual(values = c(`Up in Clonal`   = "#C4463A",
                                  `Down in Clonal` = "#3B7FB8",
                                  ns               = "grey75")) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed",
               color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "grey50") +
    labs(x = "avg log2FC (Clonal vs NonClonal)",
         y = expression(-log[10](p[adj])),
         title = paste0("CD8 ", cluster_label,
                        ": Clonal vs Non-clonal DGE")) +
    theme_cowplot(font_size = 14) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  ggsave(file.path(out_dir,
                   paste0("CD8_", file_tag, "_Clonal_vs_NonClonal_volcano.png")),
         p, width = 8, height = 6, dpi = 300, bg = "white")
  invisible(m)
}

dge_clones_in_cluster(OPIS_CD8, prf1_tem_label, "PRF1posTem", dge_prf1_dir)

###############################################################################
# 7. MODULE SCORING — HIV vs CMV in PRF1+ Tem ------------------------
#    Modules: Exhaustion, Cytotoxic_Effector, Interferon_Response (script 17)
###############################################################################

cat("\n=== Block 7: Module scoring (HIV vs CMV) in PRF1+ Tem ===\n")

# Module gene lists — copied from script 17 so this script is self-contained.
modules_hivcmv <- list(
  Cytotoxic_Effector  = c("NKG7","PRF1","GZMB","GZMH","GZMA","GNLY","CCL5",
                          "CST7","FGFBP2"),
  Exhaustion          = c("PDCD1","LAG3","TIGIT","HAVCR2","TOX","CTLA4",
                          "ENTPD1","CXCL13"),
  Interferon_Response = c("ISG15","IFITM1","IFITM2","IFI6","MX1","OAS1",
                          "OASL","STAT1","IRF7")
)
module_labels_hivcmv <- c(
  Cytotoxic_Effector  = "Cytotoxic effector",
  Exhaustion          = "Exhaustion",
  Interferon_Response = "Interferon response"
)

# --- Score modules + return per-cell long table -----------------------------
score_modules_simple <- function(seu, modules, label_col, group_col,
                                 keep_groups = NULL, min_genes = 3) {
  set.seed(42)
  DefaultAssay(seu) <- "RNA"
  modules_filt <- lapply(modules, function(g) intersect(g, rownames(seu)))
  too_small    <- vapply(modules_filt, length, integer(1)) < min_genes
  if (any(too_small)) {
    message("  Dropping modules with <", min_genes, " genes detected: ",
            paste(names(modules_filt)[too_small], collapse = ", "))
    modules_filt <- modules_filt[!too_small]
  }
  if (length(modules_filt) == 0) stop("No usable modules.")
  message("  Scoring ", length(modules_filt), " modules: ",
          paste(names(modules_filt), collapse = ", "))
  seu <- AddModuleScore(seu, features = modules_filt,
                        name = "MODULE_", nbin = 24, ctrl = 100, seed = 42)
  score_cols <- paste0("MODULE_", seq_along(modules_filt))
  df <- data.frame(
    cell    = colnames(seu),
    PID     = seu@meta.data$orig.ident,
    Cluster = as.character(seu@meta.data[[label_col]]),
    Group   = as.character(seu@meta.data[[group_col]]),
    seu@meta.data[, score_cols, drop = FALSE],
    check.names = FALSE, stringsAsFactors = FALSE)
  colnames(df)[match(score_cols, colnames(df))] <- names(modules_filt)
  if (!is.null(keep_groups)) df <- df %>% dplyr::filter(Group %in% keep_groups)
  list(seurat = seu, scores = df, modules_used = names(modules_filt))
}

# --- Cohen's d helper (script 17) -------------------------------------------
compute_cohens_d_two <- function(scores_df, modules_used,
                                 g_pos, g_neg, min_n = 3,
                                 group_col = "Group") {
  results <- list()
  for (mod in modules_used) {
    x <- scores_df %>% dplyr::filter(.data[[group_col]] == g_pos) %>% dplyr::pull(!!sym(mod))
    y <- scores_df %>% dplyr::filter(.data[[group_col]] == g_neg) %>% dplyr::pull(!!sym(mod))
    if (length(x) >= min_n && length(y) >= min_n) {
      nx <- length(x); ny <- length(y)
      sx <- sd(x);     sy <- sd(y)
      pooled <- sqrt(((nx-1)*sx^2 + (ny-1)*sy^2) / (nx+ny-2))
      cd <- if (!is.na(pooled) && pooled > 0)
        (mean(x) - mean(y)) / pooled else NA_real_
      pv <- tryCatch(wilcox.test(x, y)$p.value, error = function(e) NA_real_)
      results[[mod]] <- data.frame(
        Module   = mod,
        n_pos    = nx,        n_neg    = ny,
        mean_pos = mean(x),   mean_neg = mean(y),
        cohens_d = cd,        p_value  = pv,
        stringsAsFactors = FALSE)
    }
  }
  if (!length(results)) return(NULL)
  out <- do.call(rbind, results)
  out$p_adj <- p.adjust(out$p_value, method = "BH")
  out$star  <- ifelse(is.na(out$p_adj), "",
                      ifelse(out$p_adj < 0.001, "***",
                             ifelse(out$p_adj < 0.01, "**",
                                    ifelse(out$p_adj < 0.05, "*", ""))))
  rownames(out) <- NULL
  out
}

# Subset to PRF1+ Tem, then to HIV / CMV cells
prf1_cells <- cells_in_cluster(trex_obj_cd8, prf1_tem_label)
if (length(prf1_cells) == 0) {
  warning("No cells in PRF1+ Tem cluster — skipping Block 7.")
} else {
  prf1_obj <- subset(trex_obj_cd8, cells = prf1_cells)
  cat("  PRF1+ Tem: ", ncol(prf1_obj), " cells\n", sep = "")
  cat("  Specificity in PRF1+ Tem:\n"); print(table(prf1_obj$specificity))
  
  # Score on the full PRF1+ Tem subset, then filter to HIV/CMV
  hivcmv_out <- score_modules_simple(
    seu         = prf1_obj,
    modules     = modules_hivcmv,
    label_col   = celltype_col,
    group_col   = "specificity",
    keep_groups = c("HIV","CMV"))
  
  # CSV: per-cell module scores
  write.csv(hivcmv_out$scores,
            file.path(mod_hivcmv_dir,
                      "CD8_PRF1posTem_HIVvsCMV_ModuleScores_per_cell.csv"),
            row.names = FALSE)
  
  # CSV: per-PID summary (mean scores by PID x specificity)
  per_pid <- hivcmv_out$scores %>%
    tidyr::pivot_longer(cols = dplyr::all_of(hivcmv_out$modules_used),
                        names_to  = "Module",
                        values_to = "Score") %>%
    dplyr::group_by(PID, Group, Module) %>%
    dplyr::summarise(n_cells = dplyr::n(),
                     mean_score   = mean(Score, na.rm = TRUE),
                     median_score = median(Score, na.rm = TRUE),
                     .groups = "drop")
  write.csv(per_pid,
            file.path(mod_hivcmv_dir,
                      "CD8_PRF1posTem_HIVvsCMV_ModuleScores_per_PID.csv"),
            row.names = FALSE)
  
  # CSV: Cohen's d HIV vs CMV
  if (all(c("HIV","CMV") %in% hivcmv_out$scores$Group)) {
    cd_hivcmv <- compute_cohens_d_two(
      scores_df    = hivcmv_out$scores,
      modules_used = hivcmv_out$modules_used,
      g_pos        = "HIV", g_neg = "CMV")
    if (!is.null(cd_hivcmv))
      write.csv(cd_hivcmv,
                file.path(mod_hivcmv_dir,
                          "CD8_PRF1posTem_HIVvsCMV_CohenD.csv"),
                row.names = FALSE)
  } else {
    cd_hivcmv <- NULL
    message("Need both HIV and CMV cells for Cohen's d; got: ",
            paste(unique(hivcmv_out$scores$Group), collapse = ", "))
  }
  
  # Figure 1: violin plot (one per module), HIV vs CMV
  long_scores <- hivcmv_out$scores %>%
    tidyr::pivot_longer(cols = dplyr::all_of(hivcmv_out$modules_used),
                        names_to  = "Module",
                        values_to = "Score") %>%
    dplyr::mutate(Module_Label = module_labels_hivcmv[Module])
  
  p_violin <- ggplot(long_scores,
                     aes(x = Group, y = Score, fill = Group)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.7,
                color = "grey20", linewidth = 0.3) +
    geom_boxplot(width = 0.18, fill = "white", alpha = 0.9,
                 outlier.shape = NA) +
    facet_wrap(~ Module_Label, scales = "free_y") +
    scale_fill_manual(values = specificity_palette[c("HIV","CMV")]) +
    labs(x = NULL, y = "Module score",
         title = "PRF1+ Tem: module scores — HIV vs CMV") +
    theme_cowplot(font_size = 14) +
    theme(legend.position = "none",
          strip.background = element_rect(fill = "grey92", color = NA),
          strip.text       = element_text(face = "bold"),
          plot.title       = element_text(face = "bold", hjust = 0.5))
  ggsave(file.path(mod_hivcmv_dir,
                   "CD8_PRF1posTem_HIVvsCMV_ModuleScores_violin.png"),
         p_violin, width = 12, height = 5, dpi = 300, bg = "white")
  
  # Figure 2: Cohen's d bar
  if (!is.null(cd_hivcmv)) {
    cd_hivcmv$Module_Label <- factor(module_labels_hivcmv[cd_hivcmv$Module],
                                     levels = module_labels_hivcmv)
    p_cd <- ggplot(cd_hivcmv,
                   aes(x = Module_Label, y = cohens_d,
                       fill = ifelse(cohens_d > 0, "HIV", "CMV"))) +
      geom_col(width = 0.65, color = "grey15", linewidth = 0.35) +
      geom_text(aes(label = star),
                vjust = ifelse(cd_hivcmv$cohens_d > 0, -0.4, 1.2),
                size = 6, fontface = "bold") +
      scale_fill_manual(values = specificity_palette[c("HIV","CMV")],
                        name = "Higher in") +
      geom_hline(yintercept = 0, color = "grey30") +
      labs(x = NULL, y = "Cohen's d (HIV vs CMV)",
           title = "PRF1+ Tem: HIV vs CMV — Cohen's d") +
      theme_cowplot(font_size = 14) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))
    ggsave(file.path(mod_hivcmv_dir,
                     "CD8_PRF1posTem_HIVvsCMV_CohenD_bar.png"),
           p_cd, width = 7, height = 5, dpi = 300, bg = "white")
  }
}

###############################################################################
# 8. DGE: HIV vs CMV within PRF1+ Tem (CSV only) ----------------------
###############################################################################

cat("\n=== Block 8: DGE HIV vs CMV in PRF1+ Tem ===\n")

if (length(prf1_cells) > 0) {
  prf1_obj_full <- subset(trex_obj_cd8, cells = prf1_cells)
  hiv_n <- sum(prf1_obj_full$specificity == "HIV", na.rm = TRUE)
  cmv_n <- sum(prf1_obj_full$specificity == "CMV", na.rm = TRUE)
  cat("  PRF1+ Tem HIV cells = ", hiv_n,
      "  CMV cells = ", cmv_n, "\n", sep = "")
  
  if (hiv_n >= 5 && cmv_n >= 5) {
    Idents(prf1_obj_full) <- prf1_obj_full$specificity
    DefaultAssay(prf1_obj_full) <- "RNA"
    dge_hc <- FindMarkers(prf1_obj_full,
                          ident.1 = "HIV", ident.2 = "CMV",
                          logfc.threshold = 0.1, min.pct = 0.1)
    dge_hc$gene <- rownames(dge_hc)
    dge_hc <- dge_hc %>%
      dplyr::arrange(p_val_adj, dplyr::desc(abs(avg_log2FC))) %>%
      dplyr::select(gene, dplyr::everything())
    write.csv(dge_hc,
              file.path(dge_hivcmv_dir,
                        "CD8_PRF1posTem_HIV_vs_CMV_DGE.csv"),
              row.names = FALSE)
  } else {
    message("Skipping DGE (need >=5 cells per group).")
  }
}

###############################################################################
# 9. CLONE SIZE vs SPECIFICITY ----------------------------------
#    (a) xy plot: x = clone size, y = % HIV / % CMV / % Other
#    (b) bar chart: % HIV per clone-size bin
###############################################################################

cat("\n=== Block 9: Clone size vs specificity ===\n")

# We use the full CD8 object (every TCR-bearing cell with its specificity call)
size_spec <- OPIS_CD8@meta.data %>%
  dplyr::filter(has_TCR, !is.na(clonalFrequency)) %>%
  dplyr::select(orig.ident, CTstrict, clonalFrequency, cloneSize, specificity)

write.csv(size_spec,
          file.path(spec_size_dir, "CD8_CellLevel_CloneSize_Specificity.csv"),
          row.names = FALSE)

# (a) Per-clone-size aggregate (each clone size bin = unique clonalFrequency value)
agg_by_size <- size_spec %>%
  dplyr::count(clonalFrequency, specificity) %>%
  dplyr::group_by(clonalFrequency) %>%
  dplyr::mutate(total = sum(n),
                pct   = 100 * n / total) %>%
  dplyr::ungroup()
write.csv(agg_by_size,
          file.path(spec_size_dir,
                    "CD8_PercentSpecificity_per_CloneSize.csv"),
          row.names = FALSE)

p_xy <- ggplot(agg_by_size,
               aes(x = clonalFrequency, y = pct, color = specificity)) +
  geom_point(aes(size = total), alpha = 0.75) +
  geom_smooth(method = "loess", se = FALSE, span = 0.7,
              linewidth = 0.8, alpha = 0.6, na.rm = TRUE) +
  scale_x_log10(breaks = c(1,2,5,10,20,50,100,200,500),
                labels = c(1,2,5,10,20,50,100,200,500)) +
  scale_color_manual(values = specificity_palette, name = "Specificity") +
  scale_size_continuous(name = "Cells at\nthis clone size",
                        range = c(1, 6)) +
  labs(x = "Clone size (cells per clone, log10)",
       y = "% of cells at clone size",
       title = "CD8: clone size vs antigen specificity") +
  theme_cowplot(font_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
ggsave(file.path(spec_size_dir, "CD8_CloneSize_vs_Specificity_xy.png"),
       p_xy, width = 9, height = 6, dpi = 300, bg = "white")

# Stacked area version (clearer trend)
p_area <- ggplot(agg_by_size,
                 aes(x = clonalFrequency, y = pct, fill = specificity)) +
  geom_area(position = "fill", alpha = 0.85, color = "white",
            linewidth = 0.2) +
  scale_x_log10() +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_fill_manual(values = specificity_palette, name = "Specificity") +
  labs(x = "Clone size (log10)",
       y = "% of cells",
       title = "CD8: specificity composition by clone size") +
  theme_cowplot(font_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
ggsave(file.path(spec_size_dir, "CD8_CloneSize_vs_Specificity_area.png"),
       p_area, width = 9, height = 6, dpi = 300, bg = "white")

# (b) Bar chart: % HIV (and CMV) per clonal-expansion size bin
#     Bins follow the cloneSize labels assigned in script 16:
#     Single / Small / Medium / Large / Hyperexpanded
bin_spec <- size_spec %>%
  dplyr::filter(!is.na(cloneSize)) %>%
  dplyr::mutate(bin = sapply(strsplit(as.character(cloneSize), " "),
                             `[`, 1)) %>%
  dplyr::count(bin, specificity) %>%
  dplyr::group_by(bin) %>%
  dplyr::mutate(total = sum(n),
                pct   = 100 * n / total) %>%
  dplyr::ungroup()

bin_order <- c("Single","Small","Medium","Large","Hyperexpanded")
bin_spec$bin <- factor(bin_spec$bin,
                       levels = bin_order[bin_order %in% bin_spec$bin])
write.csv(bin_spec,
          file.path(spec_size_dir,
                    "CD8_PercentSpecificity_per_CloneSizeBin.csv"),
          row.names = FALSE)

# %HIV per bin
hiv_bin <- bin_spec %>% dplyr::filter(specificity == "HIV")
p_hiv_bin <- ggplot(hiv_bin, aes(x = bin, y = pct)) +
  geom_col(fill = specificity_palette["HIV"],
           color = "grey15", linewidth = 0.35, width = 0.7) +
  geom_text(aes(label = sprintf("%.2f%%\n(n=%d / %d)", pct, n, total)),
            vjust = -0.2, size = 3.6, lineheight = 0.95) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = "Clonal-expansion size bin",
       y = "% HIV-specific cells",
       title = "CD8: % HIV-specific cells per clone-size bin") +
  theme_cowplot(font_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
ggsave(file.path(spec_size_dir,
                 "CD8_PercentHIV_per_CloneSizeBin.png"),
       p_hiv_bin, width = 8, height = 6, dpi = 300, bg = "white")

# Optional companion: HIV + CMV side-by-side per bin
both_bin <- bin_spec %>% dplyr::filter(specificity %in% c("HIV","CMV"))
p_hc_bin <- ggplot(both_bin,
                   aes(x = bin, y = pct, fill = specificity)) +
  geom_col(position = position_dodge(width = 0.8),
           width = 0.7, color = "grey15", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.2f%%", pct)),
            position = position_dodge(width = 0.8),
            vjust = -0.4, size = 3.4) +
  scale_fill_manual(values = specificity_palette[c("HIV","CMV")],
                    name = "Specificity") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = "Clonal-expansion size bin",
       y = "% specific cells",
       title = "CD8: % HIV / CMV cells per clone-size bin") +
  theme_cowplot(font_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
ggsave(file.path(spec_size_dir,
                 "CD8_PercentHIV_CMV_per_CloneSizeBin.png"),
       p_hc_bin, width = 9, height = 6, dpi = 300, bg = "white")

###############################################################################
# 10. WITHIN CD8 HIV-SPECIFIC: DGE & MODULES BY OUD STATUS -------------
###############################################################################

cat("\n=== Block 10: HIV-specific CD8 cells by OUD status ===\n")

hiv_cells_all <- rownames(trex_obj_cd8@meta.data)[
  trex_obj_cd8$specificity == "HIV"]
cat("  Total CD8 HIV-specific cells: ", length(hiv_cells_all), "\n", sep = "")

if (length(hiv_cells_all) >= 10) {
  
  hiv_obj <- subset(trex_obj_cd8, cells = hiv_cells_all)
  
  # Per-OUD counts
  oud_tab <- table(hiv_obj@meta.data[[oud_col]], useNA = "ifany")
  cat("  HIV cells per OUD:\n"); print(oud_tab)
  write.csv(as.data.frame(oud_tab),
            file.path(hiv_oud_dir, "CD8_HIVspecific_perOUD_counts.csv"),
            row.names = FALSE)
  
  # ---- 10a. DGE: OUD+ vs OUD- among HIV-specific CD8 cells ---------------
  n_pos <- sum(hiv_obj@meta.data[[oud_col]] == oud_pos, na.rm = TRUE)
  n_neg <- sum(hiv_obj@meta.data[[oud_col]] == oud_neg, na.rm = TRUE)
  if (n_pos >= 5 && n_neg >= 5) {
    Idents(hiv_obj) <- hiv_obj@meta.data[[oud_col]]
    DefaultAssay(hiv_obj) <- "RNA"
    dge_hivoud <- FindMarkers(hiv_obj,
                              ident.1 = oud_pos, ident.2 = oud_neg,
                              logfc.threshold = 0.1, min.pct = 0.1)
    dge_hivoud$gene <- rownames(dge_hivoud)
    dge_hivoud <- dge_hivoud %>%
      dplyr::arrange(p_val_adj, dplyr::desc(abs(avg_log2FC))) %>%
      dplyr::select(gene, dplyr::everything())
    write.csv(dge_hivoud,
              file.path(hiv_oud_dir,
                        "CD8_HIVspecific_OUDpos_vs_OUDneg_DGE.csv"),
              row.names = FALSE)
    
    # Volcano
    dge_vol <- dge_hivoud %>%
      dplyr::mutate(neg_log10_padj = -log10(p_val_adj + 1e-300),
                    direction = dplyr::case_when(
                      p_val_adj < 0.05 & avg_log2FC >  0.25 ~ "Up in OUD+",
                      p_val_adj < 0.05 & avg_log2FC < -0.25 ~ "Down in OUD+",
                      TRUE ~ "ns"))
    p <- ggplot(dge_vol,
                aes(x = avg_log2FC, y = neg_log10_padj, color = direction)) +
      geom_point(alpha = 0.7, size = 1.4) +
      scale_color_manual(values = c(`Up in OUD+`   = "#D95F02",
                                    `Down in OUD+` = "#3B7FB8",
                                    ns             = "grey75")) +
      geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed",
                 color = "grey50") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed",
                 color = "grey50") +
      labs(x = "avg log2FC (OUD+ vs OUD-)",
           y = expression(-log[10](p[adj])),
           title = "CD8 HIV-specific: OUD+ vs OUD- DGE") +
      theme_cowplot(font_size = 14) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))
    ggsave(file.path(hiv_oud_dir,
                     "CD8_HIVspecific_OUDpos_vs_OUDneg_volcano.png"),
           p, width = 8, height = 6, dpi = 300, bg = "white")
  } else {
    message("  Too few HIV cells per OUD group for DGE (OUD+=", n_pos,
            ", OUD-=", n_neg, ").")
  }
  
  # ---- 10b. Module scoring on HIV-specific CD8, OUD+ vs OUD- ------------
  hiv_mods_out <- score_modules_simple(
    seu         = hiv_obj,
    modules     = modules_hivcmv,
    label_col   = celltype_col,
    group_col   = oud_col,
    keep_groups = c(oud_pos, oud_neg))
  
  write.csv(hiv_mods_out$scores,
            file.path(hiv_oud_dir,
                      "CD8_HIVspecific_ModuleScores_per_cell.csv"),
            row.names = FALSE)
  
  # Per-PID summary
  per_pid_hiv <- hiv_mods_out$scores %>%
    tidyr::pivot_longer(cols = dplyr::all_of(hiv_mods_out$modules_used),
                        names_to  = "Module",
                        values_to = "Score") %>%
    dplyr::group_by(PID, Group, Module) %>%
    dplyr::summarise(n_cells = dplyr::n(),
                     mean_score   = mean(Score, na.rm = TRUE),
                     median_score = median(Score, na.rm = TRUE),
                     .groups = "drop")
  write.csv(per_pid_hiv,
            file.path(hiv_oud_dir,
                      "CD8_HIVspecific_ModuleScores_per_PID.csv"),
            row.names = FALSE)
  
  if (all(c(oud_pos, oud_neg) %in% hiv_mods_out$scores$Group)) {
    cd_hiv <- compute_cohens_d_two(
      scores_df    = hiv_mods_out$scores,
      modules_used = hiv_mods_out$modules_used,
      g_pos        = oud_pos, g_neg = oud_neg)
    if (!is.null(cd_hiv))
      write.csv(cd_hiv,
                file.path(hiv_oud_dir,
                          "CD8_HIVspecific_ModuleScores_CohenD.csv"),
                row.names = FALSE)
    
    long_hiv <- hiv_mods_out$scores %>%
      tidyr::pivot_longer(cols = dplyr::all_of(hiv_mods_out$modules_used),
                          names_to  = "Module",
                          values_to = "Score") %>%
      dplyr::mutate(Module_Label = module_labels_hivcmv[Module])
    
    p_v <- ggplot(long_hiv,
                  aes(x = Group, y = Score, fill = Group)) +
      geom_violin(scale = "width", trim = TRUE, alpha = 0.7,
                  color = "grey20", linewidth = 0.3) +
      geom_boxplot(width = 0.18, fill = "white", alpha = 0.9,
                   outlier.shape = NA) +
      facet_wrap(~ Module_Label, scales = "free_y") +
      scale_fill_manual(values = oud_palette) +
      labs(x = NULL, y = "Module score",
           title = "CD8 HIV-specific: module scores by OUD") +
      theme_cowplot(font_size = 14) +
      theme(legend.position  = "none",
            strip.background = element_rect(fill = "grey92", color = NA),
            strip.text       = element_text(face = "bold"),
            plot.title       = element_text(face = "bold", hjust = 0.5))
    ggsave(file.path(hiv_oud_dir,
                     "CD8_HIVspecific_ModuleScores_byOUD_violin.png"),
           p_v, width = 12, height = 5, dpi = 300, bg = "white")
  }
}

###############################################################################
# 11. SHARED HIV-SPECIFIC CLONES ACROSS PARTICIPANTS -----------------
###############################################################################

cat("\n=== Block 11: shared (HIV-specific) clones across participants ===\n")

shared_summary <- function(meta, file_tag, n_min_pids = 2) {
  meta <- meta %>%
    dplyr::filter(!is.na(CTstrict), !is.na(orig.ident))
  if (nrow(meta) == 0) {
    message("No clones for ", file_tag); return(invisible())
  }
  
  per_clone <- meta %>%
    dplyr::group_by(CTstrict) %>%
    dplyr::summarise(n_cells = dplyr::n(),
                     n_PIDs  = dplyr::n_distinct(orig.ident),
                     PIDs    = paste(sort(unique(orig.ident)), collapse = ";"),
                     n_OUDpos = sum(.data[[oud_col]] == oud_pos, na.rm = TRUE),
                     n_OUDneg = sum(.data[[oud_col]] == oud_neg, na.rm = TRUE),
                     .groups  = "drop") %>%
    dplyr::arrange(dplyr::desc(n_PIDs), dplyr::desc(n_cells))
  
  write.csv(per_clone,
            file.path(shared_clones_dir,
                      paste0(file_tag, "_PerClone_Sharing.csv")),
            row.names = FALSE)
  
  shared <- per_clone %>% dplyr::filter(n_PIDs >= n_min_pids)
  write.csv(shared,
            file.path(shared_clones_dir,
                      paste0(file_tag, "_SharedClones_min", n_min_pids,
                             "PIDs.csv")),
            row.names = FALSE)
  
  cat("  ", file_tag, ": ", nrow(per_clone), " unique clones, ",
      nrow(shared), " shared across >= ", n_min_pids, " PIDs\n", sep = "")
  
  # Sharing distribution plot (how many clones are seen in 1, 2, ... PIDs)
  pid_dist <- per_clone %>%
    dplyr::count(n_PIDs, name = "n_clones") %>%
    dplyr::arrange(n_PIDs)
  if (nrow(pid_dist) > 0) {
    p <- ggplot(pid_dist, aes(x = factor(n_PIDs), y = n_clones)) +
      geom_col(fill = "#3B7FB8", color = "grey15", linewidth = 0.35,
               width = 0.7) +
      geom_text(aes(label = n_clones), vjust = -0.3, size = 3.8) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
      labs(x = "# of participants sharing the clone",
           y = "# of unique clones",
           title = paste0(file_tag, ": clone sharing distribution")) +
      theme_cowplot(font_size = 14) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))
    ggsave(file.path(shared_clones_dir,
                     paste0(file_tag, "_SharingDistribution.png")),
           p, width = 7, height = 5, dpi = 300, bg = "white")
  }
  
  # If any shared, write a sample-by-clone wide table
  if (nrow(shared) > 0) {
    cell_lvl <- meta %>%
      dplyr::filter(CTstrict %in% shared$CTstrict) %>%
      dplyr::count(CTstrict, orig.ident) %>%
      tidyr::pivot_wider(names_from = orig.ident, values_from = n,
                         values_fill = 0)
    write.csv(cell_lvl,
              file.path(shared_clones_dir,
                        paste0(file_tag, "_SharedClones_byPID_matrix.csv")),
              row.names = FALSE)
  }
  invisible(per_clone)
}

# (a) HIV-specific clones
hiv_meta <- trex_obj_cd8@meta.data[
  trex_obj_cd8$specificity == "HIV", ,drop = FALSE]
shared_summary(hiv_meta, "CD8_HIVspecific")

# (b) CMV-specific (for context)
cmv_meta <- trex_obj_cd8@meta.data[
  trex_obj_cd8$specificity == "CMV", ,drop = FALSE]
shared_summary(cmv_meta, "CD8_CMVspecific")

# (c) All CD8 clones
all_meta <- trex_obj_cd8@meta.data
shared_summary(all_meta, "CD8_All")

###############################################################################
# 12. STATS CSVs --------------------------------------------------
#     (a) % HIV-specific cells in CD8+ PRF1+ Tem per PID
#     (b) Cytotoxic module scores in CD8+ PRF1+ Tem HIV-specific cells per PID
#     (c) CD8+ Gini coefficient of clonality per PID
###############################################################################

cat("\n=== Block 12: Stats CSVs ===\n")

# (a) % HIV-specific in PRF1+ Tem per PID -----------------------------------
prf1_meta <- OPIS_CD8@meta.data %>%
  dplyr::filter(.data[[celltype_col]] == prf1_tem_label)

stat_a <- prf1_meta %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::summarise(
    n_PRF1posTem = dplyr::n(),
    n_HIV        = sum(specificity == "HIV", na.rm = TRUE),
    n_CMV        = sum(specificity == "CMV", na.rm = TRUE),
    pct_HIV      = 100 * n_HIV / dplyr::n(),
    pct_CMV      = 100 * n_CMV / dplyr::n(),
    .groups      = "drop") %>%
  dplyr::left_join(
    prf1_meta %>% dplyr::distinct(orig.ident, .data[[oud_col]]),
    by = "orig.ident") %>%
  dplyr::arrange(orig.ident)
write.csv(stat_a,
          file.path(stats_dir,
                    "Stat_A_PercentHIV_in_PRF1posTem_perParticipant.csv"),
          row.names = FALSE)
cat("  (a) wrote Stat_A_PercentHIV_in_PRF1posTem_perParticipant.csv (",
    nrow(stat_a), " participants)\n", sep = "")

# (b) Cytotoxic module score for HIV-specific cells in PRF1+ Tem per PID ----
# Reuse the per-cell scores produced in Block 7 if they exist; else recompute.
if (exists("hivcmv_out") && !is.null(hivcmv_out)) {
  
  hiv_scores <- hivcmv_out$scores %>% dplyr::filter(Group == "HIV")
  
  # Always-present: Cytotoxic_Effector (script aborts in Block 7 if missing).
  # Exhaustion / Interferon are also requested but may have been dropped if too
  # few genes were detected — handle each independently.
  available_mods <- intersect(
    c("Cytotoxic_Effector", "Exhaustion", "Interferon_Response"),
    names(hiv_scores))
  
  cyto_per_pid <- hiv_scores %>%
    dplyr::group_by(PID) %>%
    dplyr::summarise(
      n_cells          = dplyr::n(),
      mean_Cytotoxic   = if ("Cytotoxic_Effector" %in% available_mods)
        mean(Cytotoxic_Effector, na.rm = TRUE) else NA_real_,
      median_Cytotoxic = if ("Cytotoxic_Effector" %in% available_mods)
        median(Cytotoxic_Effector, na.rm = TRUE) else NA_real_,
      sd_Cytotoxic     = if ("Cytotoxic_Effector" %in% available_mods)
        sd(Cytotoxic_Effector, na.rm = TRUE) else NA_real_,
      mean_Exhaustion  = if ("Exhaustion" %in% available_mods)
        mean(Exhaustion, na.rm = TRUE) else NA_real_,
      mean_Interferon  = if ("Interferon_Response" %in% available_mods)
        mean(Interferon_Response, na.rm = TRUE) else NA_real_,
      .groups = "drop")
  
  pid_oud <- OPIS_CD8@meta.data %>%
    dplyr::distinct(orig.ident, .data[[oud_col]]) %>%
    dplyr::rename(PID = orig.ident)
  cyto_per_pid <- cyto_per_pid %>% dplyr::left_join(pid_oud, by = "PID")
  
  write.csv(cyto_per_pid,
            file.path(stats_dir,
                      "Stat_B_PRF1posTem_HIV_CytotoxicModuleScore_perParticipant.csv"),
            row.names = FALSE)
  cat("  (b) wrote Stat_B_PRF1posTem_HIV_CytotoxicModuleScore_perParticipant.csv (",
      nrow(cyto_per_pid), " participants)\n", sep = "")
} else {
  message("Skipping stat (b) — Block 7 scoring did not run.")
}

# (c) CD8+ Gini coefficient of clonality per PID ---------------------------
gini_df <- OPIS_CD8@meta.data %>%
  dplyr::filter(has_TCR, !is.na(CTstrict)) %>%
  dplyr::group_by(orig.ident, CTstrict) %>%
  dplyr::summarise(clone_n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::summarise(
    n_cells_with_TCR = sum(clone_n),
    n_unique_clones  = dplyr::n(),
    gini_clonality   = gini(clone_n),
    .groups = "drop")

pid_oud <- OPIS_CD8@meta.data %>%
  dplyr::distinct(orig.ident, .data[[oud_col]])
gini_df <- gini_df %>% dplyr::left_join(pid_oud, by = "orig.ident")

write.csv(gini_df,
          file.path(stats_dir, "Stat_C_CD8_Gini_perParticipant.csv"),
          row.names = FALSE)
cat("  (c) wrote Stat_C_CD8_Gini_perParticipant.csv (",
    nrow(gini_df), " participants)\n", sep = "")

# Optional: quick visual of Gini by OUD for QC
if (!all(is.na(gini_df[[oud_col]]))) {
  p_g <- ggplot(gini_df %>% dplyr::filter(!is.na(.data[[oud_col]])),
                aes(x = .data[[oud_col]], y = gini_clonality,
                    fill = .data[[oud_col]])) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.55,
                 color = "grey20") +
    geom_jitter(width = 0.15, size = 2.5, alpha = 0.85, color = "grey25") +
    scale_fill_manual(values = oud_palette) +
    labs(x = NULL, y = "Gini coefficient (CD8 clonality)",
         title = "CD8 clonality Gini per participant, by OUD") +
    theme_cowplot(font_size = 14) +
    theme(legend.position = "none",
          plot.title      = element_text(face = "bold", hjust = 0.5))
  ggsave(file.path(stats_dir, "Stat_C_CD8_Gini_byOUD_boxplot.png"),
         p_g, width = 6, height = 5, dpi = 300, bg = "white")
}

###############################################################################
# DONE
###############################################################################

cat("\n", strrep("=", 60), "\n", sep = "")
cat("Script 18 complete.\n")
cat("All outputs under:\n  ", manu.root, "\n", sep = "")
cat(strrep("=", 60), "\n", sep = "")