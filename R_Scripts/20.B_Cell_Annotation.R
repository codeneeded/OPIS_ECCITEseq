# =============================================================================
# OPIS ECCITEseq — B-cell FINAL MERGE + publication heatmaps + DGE + pathways
#
# Request from Ben:
#   1. Naive B (resting) + Naive B (IEG-high)            -> "Naive B"
#   2. Transitional B + Transitional/Naive (T2-like)     -> "Transitional B"
#   3. keep Intermediate B (HLA-G+)
#   4. keep Memory B
#   5. keep Atypical B (ABC/DN2)
#   -> re-do the same-feature heatmaps (publication ready) to justify the calls
#   -> DGE on the merged clusters using the MAST test
#   -> pathway analysis on the merged clusters
#
# Reads:  OPIS_BCELL_Annotated.qs2
# Writes: <Bcells>/final_merged/  (heatmaps, DGE_MAST, pathway) + a new .qs2
#         Nothing is overwritten or deleted on disk.
#
# ---- RUNTIME ASSUMPTIONS / NOTES --------------------------------------------
#   [A1] OPIS_BCELL_Annotated.qs2 lives in `load.path`; it carries the
#        Bcell_Annotation column with the 7 pre-merge labels and a wnn.umap.
#   [A2] DGE/pathway are RNA-based. MAST must be installed (BiocManager::install
#        ("MAST")). Pathway needs clusterProfiler + org.Hs.eg.db (+ msigdbr,
#        enrichplot for GSEA/plots). Missing optional pkgs are skipped with a
#        message rather than crashing the run.
#   [A3] MAST DGE is run ONCE with logfc.threshold = 0 on expressed genes so the
#        same table feeds both the DGE report (filtered) and GSEA (full ranking).
#        This is the slow step (~minutes to tens of minutes on ~6–7k cells).
#   [A4] Heatmaps reuse the EXACT ADT/RNA panels Ben reviewed, same 0–1 row
#        scaling, columns in fixed developmental order (no dendrogram) for a
#        clean publication figure.
# =============================================================================

library(Seurat)
library(SeuratExtend)
library(qs2)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(grid)

# ---- Paths -------------------------------------------------------------------
load.path       <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/saved_R_data"
subclust.root   <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/subclustering/pre_annotation"
bcell.save.path <- file.path(subclust.root, "Bcells")

final.dir   <- file.path(bcell.save.path, "final_merged")
hm.dir      <- file.path(final.dir, "heatmaps")
dge.dir     <- file.path(final.dir, "DGE_MAST")
ora.dir     <- file.path(final.dir, "pathway", "ORA_GO_BP")
gsea.dir    <- file.path(final.dir, "pathway", "GSEA_Hallmark")
for (d in c(hm.dir, dge.dir, ora.dir, gsea.dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- Load object -------------------------------------------------------------
message("Loading OPIS_BCELL_Annotated.qs2 ...")
OPIS_BCELL <- qs_read(file.path(load.path, "OPIS_BCELL_Annotated.qs2"))

# Seurat v5 split-layer safety: join before differential testing
if (inherits(OPIS_BCELL[["RNA"]], "Assay5")) {
  message("Joining RNA layers (Seurat v5) ...")
  OPIS_BCELL[["RNA"]] <- tryCatch(JoinLayers(OPIS_BCELL[["RNA"]]),
                                  error = function(e) OPIS_BCELL[["RNA"]])
}

# =============================================================================
# 1.  MERGE clusters into the 5 final lineages
# =============================================================================
remap <- c(
  "Naive B (resting)"            = "Naive B",
  "Naive B (IEG-high)"           = "Naive B",
  "Transitional B"               = "Transitional B",
  "Transitional/Naive (T2-like)" = "Transitional B",
  "Intermediate B (HLA-G+)"      = "Intermediate B (HLA-G+)",
  "Memory B"                     = "Memory B",
  "Atypical B (ABC/DN2)"         = "Atypical B (ABC/DN2)"
)
final.levels <- c("Transitional B", "Naive B",
                  "Intermediate B (HLA-G+)", "Memory B", "Atypical B (ABC/DN2)")

src <- as.character(OPIS_BCELL$Bcell_Annotation)
if (any(is.na(remap[src])))
  message("WARNING: unmapped annotation(s): ",
          paste(unique(src[is.na(remap[src])]), collapse = ", "))
OPIS_BCELL$Bcell_Final <- factor(unname(remap[src]), levels = final.levels)
Idents(OPIS_BCELL) <- "Bcell_Final"

group_var   <- "Bcell_Final"
grp         <- droplevels(factor(OPIS_BCELL@meta.data[[group_var]]))
lvls        <- levels(grp)
cell_counts <- table(grp)
col_lab     <- paste0(lvls, "\n(n=", as.integer(cell_counts[lvls]), ")")

message("Final clusters:")
print(cell_counts)
print(table(OPIS_BCELL$Bcell_Annotation, OPIS_BCELL$Bcell_Final))

# Cluster colours (consistent across UMAP / heatmaps)
clust_cols <- setNames(
  c("#4C72B0", "#55A868", "#C44E52", "#8172B3", "#CCB974")[seq_along(lvls)], lvls)

# Merged UMAP (existing wnn.umap embedding; no re-clustering)
# SeuratExtend boxed style, matching the subclustering pipeline's DimPlot2 calls
umap.final <- DimPlot2(OPIS_BCELL, reduction = "wnn.umap",
                       group.by = "Bcell_Final", cols = clust_cols,
                       label = TRUE, repel = TRUE, box = TRUE, label.size = 5) +
  ggtitle("B cells | Final merged annotation")
ggsave(file.path(final.dir, "Bcell_Final_UMAP.png"), umap.final,
       dpi = 500, width = 10, height = 8, bg = "white")

# Save merged object (new file; PreAnnotation / Annotated untouched)
qs_save(OPIS_BCELL, file = file.path(load.path, "OPIS_BCELL_FinalMerged.qs2"))
message("Saved: ", file.path(load.path, "OPIS_BCELL_FinalMerged.qs2"))

# =============================================================================
# 2.  Marker panels (identical to the panels Ben reviewed)
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
# 3.  Helpers (feature resolver, scaling, panel prep, averaging)
# =============================================================================
.norm <- function(x) toupper(gsub("[-_. ]", "", x))

resolve_feature <- function(label, aliases, available) {
  cand <- c(label, if (nzchar(aliases)) strsplit(aliases, ",")[[1]] else character(0))
  cand <- trimws(cand)
  for (c0 in cand) if (c0 %in% available) return(c0)
  av_norm <- .norm(available)
  for (c0 in cand) {
    hit <- available[av_norm == .norm(c0)]
    if (length(hit)) return(hit[1])
  }
  NA_character_
}

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

prepare_panel <- function(markers, cat_levels, available, modality) {
  res <- markers %>%
    mutate(feature = mapply(resolve_feature, label, aliases,
                            MoreArgs = list(available = available)))
  miss <- res %>% filter(is.na(feature))
  if (nrow(miss)) message("[", modality, "] not found, skipped: ",
                          paste(miss$label, collapse = ", "))
  res <- res %>% filter(!is.na(feature)) %>%
    mutate(category = factor(category, levels = cat_levels)) %>% arrange(category)
  dup <- duplicated(res$feature)
  if (any(dup)) message("[", modality, "] duplicate feature kept once: ",
                        paste(unique(res$feature[dup]), collapse = ", "))
  res[!dup, , drop = FALSE]
}

# =============================================================================
# 4.  Average-expression matrices for the curated panels
# =============================================================================
# ADT: plain rowMeans on DSB "data" (NO expm1)
DefaultAssay(OPIS_BCELL) <- "ADT"
adt_avail <- rownames(OPIS_BCELL[["ADT"]])
adt_panel <- prepare_panel(adt_markers, adt_cat_levels, adt_avail, "ADT")
adt_data  <- GetAssayData(OPIS_BCELL, assay = "ADT", slot = "data")[adt_panel$feature, , drop = FALSE]
adt_avg   <- sapply(lvls, function(cl) rowMeans(adt_data[, which(grp == cl), drop = FALSE]))
rownames(adt_avg) <- adt_panel$label; colnames(adt_avg) <- lvls
adt_scaled <- scale_01(adt_avg)

# RNA: AverageExpression on log1p "data"
DefaultAssay(OPIS_BCELL) <- "RNA"
rna_avail <- rownames(OPIS_BCELL[["RNA"]])
rna_panel <- prepare_panel(rna_markers, rna_cat_levels, rna_avail, "RNA")
rna_avg_raw <- AverageExpression(OPIS_BCELL, assays = "RNA",
                                 features = rna_panel$feature,
                                 group.by = group_var, slot = "data")$RNA
ci <- match(.norm(lvls), .norm(colnames(rna_avg_raw)))
if (anyNA(ci)) ci <- seq_along(lvls)
rna_avg <- rna_avg_raw[rna_panel$feature, ci, drop = FALSE]
rownames(rna_avg) <- rna_panel$label; colnames(rna_avg) <- lvls
rna_scaled <- scale_01(rna_avg)

# CSVs behind the panels
write.csv(round(adt_avg, 4),    file.path(hm.dir, "AvgExpr_ADT_panel_raw.csv"))
write.csv(round(adt_scaled, 4), file.path(hm.dir, "AvgExpr_ADT_panel_scaled.csv"))
write.csv(round(rna_avg, 4),    file.path(hm.dir, "AvgExpr_RNA_panel_raw.csv"))
write.csv(round(rna_scaled, 4), file.path(hm.dir, "AvgExpr_RNA_panel_scaled.csv"))

# =============================================================================
# 5.  Publication heatmaps (fixed developmental column order, no cell text)
# =============================================================================
pub_heatmap <- function(scaled, row_cat, title, file_base, italic_rows) {
  col_fun <- colorRamp2(seq(0, 1, length.out = 9), viridis(9, option = "A"))
  rn_gp   <- gpar(fontsize = 10, fontface = if (italic_rows) "italic" else "plain")
  ht <- Heatmap(
    scaled, name = "Scaled\n(0–1)", col = col_fun,
    cluster_rows = FALSE, row_split = row_cat,
    row_title_rot = 0, row_title_gp = gpar(fontsize = 11, fontface = "bold"),
    row_gap = unit(2.2, "mm"), row_names_gp = rn_gp,
    cluster_columns = FALSE, column_order = lvls,
    column_labels = col_lab, column_names_rot = 30,
    column_names_gp = gpar(fontsize = 11), column_names_centered = TRUE,
    rect_gp = gpar(col = "white", lwd = 0.6),
    border = TRUE,
    column_title = title, column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    heatmap_legend_param = list(at = c(0, 0.5, 1),
                                legend_height = unit(2.4, "cm"))
  )
  h_in <- 2.4 + 0.205 * nrow(scaled)
  w_in <- 4.0 + 1.05 * ncol(scaled)
  png(paste0(file_base, ".png"), width = w_in, height = h_in, units = "in", res = 500)
  draw(ht, merge_legend = TRUE); dev.off()
  message("Saved: ", file_base, ".png")
}

pub_heatmap(adt_scaled, adt_panel$category,
            "B cells | ADT marker panel",
            file.path(hm.dir, "Bcell_ADT_marker_heatmap_FINAL"), italic_rows = FALSE)
pub_heatmap(rna_scaled, rna_panel$category,
            "B cells | RNA marker panel",
            file.path(hm.dir, "Bcell_RNA_marker_heatmap_FINAL"), italic_rows = TRUE)

# =============================================================================
# 6.  Expressed-gene universe (shared by MAST + pathway ORA)
# =============================================================================
DefaultAssay(OPIS_BCELL) <- "RNA"
Idents(OPIS_BCELL) <- "Bcell_Final"

counts <- GetAssayData(OPIS_BCELL, assay = "RNA", slot = "counts")
det    <- Matrix::rowSums(counts > 0)
expressed_genes <- names(det)[det >= 0.01 * ncol(OPIS_BCELL)]
writeLines(expressed_genes, file.path(dge.dir, "expressed_genes_universe.txt"))
message("Genes tested (>=1% detection): ", length(expressed_genes))

# =============================================================================
# 7.  DGE — MAST (one-vs-rest), full pass
#     logfc.threshold = 0 so this one table serves both the DGE report and GSEA.
#     This is the slow step (minutes to tens of minutes).
# =============================================================================
if (!requireNamespace("MAST", quietly = TRUE))
  stop("MAST not installed. Install with: BiocManager::install('MAST')")

mast.latent <- NULL   # set to "nFeature_RNA" to add CDR as a covariate (MAST best practice)

message("Running FindAllMarkers (MAST) ...")
t0 <- Sys.time()
dge_all <- FindAllMarkers(
  OPIS_BCELL, assay = "RNA", slot = "data",
  test.use        = "MAST",
  latent.vars     = mast.latent,
  features        = expressed_genes,
  only.pos        = FALSE,
  logfc.threshold = 0,
  min.pct         = 0.01,
  verbose         = TRUE
)
message("MAST done in ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min.")

dge_all <- dge_all %>% mutate(cluster = factor(cluster, levels = lvls))
write.csv(dge_all, file.path(dge.dir, "DGE_MAST_full_allClusters.csv"), row.names = FALSE)

sig <- dge_all %>% filter(p_val_adj < 0.05, avg_log2FC > 0.25) %>%
  arrange(cluster, desc(avg_log2FC))
write.csv(sig, file.path(dge.dir, "DGE_MAST_significant_up.csv"), row.names = FALSE)
for (cl in lvls) {
  write.csv(sig %>% filter(cluster == cl),
            file.path(dge.dir, paste0("DGE_MAST_up_",
                                      gsub("[^A-Za-z0-9]+", "_", cl), ".csv")),
            row.names = FALSE)
}

# =============================================================================
# 8.  Top-DGE marker heatmap (row z-score) — supports the annotations
#     Independent: reloads the significant DGE table from CSV if not in memory.
# =============================================================================
if (!exists("sig"))
  sig <- read.csv(file.path(dge.dir, "DGE_MAST_significant_up.csv")) %>%
  mutate(cluster = factor(cluster, levels = lvls))

top_n     <- 10
top_genes <- sig %>% group_by(cluster) %>% slice_head(n = top_n) %>% pull(gene) %>% unique()

if (length(top_genes) >= 2) {
  tg_avg <- AverageExpression(OPIS_BCELL, assays = "RNA", features = top_genes,
                              group.by = group_var, slot = "data")$RNA
  ci2 <- match(.norm(lvls), .norm(colnames(tg_avg)))
  if (anyNA(ci2)) ci2 <- seq_along(lvls)
  tg_avg <- tg_avg[, ci2, drop = FALSE]; colnames(tg_avg) <- lvls
  tg_z <- t(scale(t(as.matrix(tg_avg)))); tg_z[is.na(tg_z)] <- 0
  row_assign <- factor(lvls[max.col(tg_z, ties.method = "first")], levels = lvls)
  ord <- order(row_assign)
  ht_top <- Heatmap(
    tg_z[ord, , drop = FALSE], name = "Row z",
    col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
    cluster_rows = FALSE, row_split = row_assign[ord], row_title_rot = 0,
    cluster_columns = FALSE, column_order = lvls,
    column_labels = col_lab, column_names_rot = 30, column_names_centered = TRUE,
    row_names_gp = gpar(fontsize = 8, fontface = "italic"),
    rect_gp = gpar(col = "white", lwd = 0.4), border = TRUE,
    column_title = paste0("B cells | Top ", top_n, " MAST markers / cluster"))
  h_in <- 2.5 + 0.16 * nrow(tg_z); w_in <- 4 + 1.05 * ncol(tg_z)
  png(file.path(hm.dir, "Bcell_topDGE_MAST_heatmap.png"),
      width = w_in, height = h_in, units = "in", res = 500)
  draw(ht_top); dev.off()
  message("Saved top-DGE marker heatmap.")
}
