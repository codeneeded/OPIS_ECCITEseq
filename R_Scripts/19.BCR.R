# =============================================================================
# OPIS ECCITEseq - Subclustering: B CELLS
# Pre-Annotation Pipeline
# (mirrors the CD4 / NK_CD8 subclustering workflow)
# =============================================================================

# ---- Libraries ---------------------------------------------------------------
library(Seurat)
library(scater)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(hdf5r)
library(data.table)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(patchwork)
library(reticulate)
library(circlize)
library(ComplexHeatmap)
library(readxl)
library(scCustomize)
library(Polychrome)
library(viridis)
library(SeuratExtend)
library(qs2)

# ---- Paths -------------------------------------------------------------------
load.path  <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/saved_R_data"

# Output root: subclustering/pre_annotation/<subset>/
subclust.root <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/subclustering/pre_annotation"

bcell.save.path <- file.path(subclust.root, "Bcells")

# Create all output subdirectories
for (sub in c("plots/ADT/violin", "plots/ADT/violin_OUD",
              "plots/ADT/feature",
              "plots/RNA/violin", "plots/RNA/violin_OUD",
              "plots/RNA/feature",
              "plots/UMAP",
              "qs2_objects")) {
  dir.create(file.path(bcell.save.path, sub), recursive = TRUE, showWarnings = FALSE)
}

# ---- Load full Seurat object -------------------------------------------------
message("Loading OPIS_ALL ...")
OPIS_ALL <- qs_read(file.path(load.path, "OPIS_ALL_PostAnnotation.qs2"))

# ---- Feature lists -----------------------------------------------------------
# Main RNA panel (same broad panel used in the CD4 / NK_CD8 script -- it already
# carries the key B / plasma markers: CD19, MS4A1, CD79A, CR2, FCRL5, IGH*, JCHAIN,
# XBP1, PRDM1, TCF4, FCGR2B, BACH2, CD200, etc.)
rna.features <- c(
  'ASCL2','BATF','BATF3','BCL6','BACH2','C1QBP','CCL2','CCL3','CCL4L2','CCL5','CCND3','CD14','CD19','CD1C',
  'CD200','CD27','CD3D','CD3E','CD36','CD4','CD40','CD40LG','CD70','CD7','CD79A','CD8A','CD8B', "CCR5",
  'CLEC9A','CR2','CTLA4','CTSW','CXCL8','CXCR3','CXCR5','EBI3','ENTPD1','FABP5','FCGR2B','FCGR3A',
  'FCRL5','FOXP3','GNLY','GP1BA','GP9','GATA3','GZMK','HAVCR2','HIF1A','HIST1H4C','HLA-DPA1',
  'HLA-DRA','HLA-DRB1','ICOS','IFI30','IFNG','IGFBP2','IGFBP4','IGHA1','IGHA2','IGHG1','IGHG2','IGHG3','IGHG4','IGHM','IKZF2','IL10',
  'IL17A','IL18BP','IL18RAP','IL1B','IL21','IL2RA','IL2RB','IRF4','IRF8','ITGAX','JCHAIN','KLRB1',
  'KLRD1','KLRC1','LAG3','LDHA','LGALS1','LTA','LTB','MAF','MAL','MALAT1','MIR155HG','MKI67',
  'MT-ND1','MT-ND5','MS4A1','NELL2','NCAM1','NKG7','NR4A1','PDCD1','PF4','PPBP','PRDM1','PRF1',
  'RORC','SELL','SERPINA1','SERPING1','SH2D1A','TCF4','TIGIT','TNF','TNFAIP2','TNFRSF18',
  'TNFRSF4','TNFRSF9','TOX','TBX21','TRBC1','TRDC','TRDV1','TRDV2','TRGC1','TRGC2','TRGV9','XBP1',
  'XCL1','XCL2','ZBTB16','ZEB2'
)

prots <- rownames(OPIS_ALL@assays$ADT)

# =============================================================================
# HELPER FUNCTION: Generate violin + feature plots for a Seurat subset
# =============================================================================
generate_annotation_plots <- function(obj, save.base, subset.label) {
  
  prots.present <- rownames(obj@assays$ADT)
  
  adt_vln_dir     <- file.path(save.base, "plots/ADT/violin")
  adt_vln_OUD_dir <- file.path(save.base, "plots/ADT/violin_OUD")
  adt_feature_dir <- file.path(save.base, "plots/ADT/feature")
  rna_vln_dir     <- file.path(save.base, "plots/RNA/violin")
  rna_vln_OUD_dir <- file.path(save.base, "plots/RNA/violin_OUD")
  rna_feature_dir <- file.path(save.base, "plots/RNA/feature")
  
  # ---- ADT plots -------------------------------------------------------------
  message(paste0("[", subset.label, "] Generating ADT plots ..."))
  DefaultAssay(obj) <- "ADT"
  
  for (i in prots.present) {
    
    # 1) Split by OUD_status
    tryCatch({
      vln.split <- VlnPlot2(obj, features = i, cols = "default",
                            split.by = "OUD_status",
                            stat.method = "wilcox.test",
                            show.mean = TRUE,
                            mean_colors = c("red","blue")) +
        ggtitle(paste(subset.label, "| ADT |", i, "| Split by OUD_status"))
      ggsave(file.path(adt_vln_OUD_dir, paste0(i, "_ADT_VLNplot_splitBy_OUD_status.png")),
             plot = vln.split, dpi = 500, width = 14, height = 8, bg = "white")
    }, error = function(e) message("  Skipping ADT split vln for ", i, ": ", e$message))
    
    # 2) Non-split violin
    tryCatch({
      vln.pl <- VlnPlot2(obj, features = i, cols = "default",
                         show.mean = TRUE, mean_colors = c("red","blue")) +
        ggtitle(paste(subset.label, "| ADT |", i))
      ggsave(file.path(adt_vln_dir, paste0(i, "_ADT_VLNplot.png")),
             plot = vln.pl, dpi = 500, width = 14, height = 8, bg = "white")
    }, error = function(e) message("  Skipping ADT vln for ", i, ": ", e$message))
    
    # 3) Feature plot
    tryCatch({
      pal <- viridis(n = 10, option = "A")
      fea.pl <- FeaturePlot_scCustom(obj, reduction = "wnn.umap",
                                     features = i, colors_use = pal, order = TRUE)
      ggsave(file.path(adt_feature_dir, paste0(i, "_ADT_Featureplot_Magma.png")),
             plot = fea.pl, dpi = 500, width = 8, bg = "white")
    }, error = function(e) message("  Skipping ADT feature for ", i, ": ", e$message))
  }
  
  # ---- RNA plots -------------------------------------------------------------
  message(paste0("[", subset.label, "] Generating RNA plots ..."))
  DefaultAssay(obj) <- "RNA"
  
  for (i in rna.features) {
    if (!i %in% rownames(obj[["RNA"]])) next
    
    # 1) Split by OUD_status
    tryCatch({
      vln.split <- VlnPlot2(obj, features = i, cols = "default",
                            split.by = "OUD_status",
                            stat.method = "wilcox.test",
                            show.mean = TRUE,
                            mean_colors = c("red","blue")) +
        ggtitle(paste(subset.label, "| RNA |", i, "| Split by OUD_status"))
      ggsave(file.path(rna_vln_OUD_dir, paste0(i, "_RNA_VLNplot_splitBy_OUD_status.png")),
             plot = vln.split, dpi = 500, width = 14, height = 8, bg = "white")
    }, error = function(e) message("  Skipping RNA split vln for ", i, ": ", e$message))
    
    # 2) Non-split violin
    tryCatch({
      vln.pl <- VlnPlot2(obj, features = i, cols = "default",
                         show.mean = TRUE, mean_colors = c("red","blue")) +
        ggtitle(paste(subset.label, "| RNA |", i))
      ggsave(file.path(rna_vln_dir, paste0(i, "_RNA_VLNplot.png")),
             plot = vln.pl, dpi = 500, width = 14, height = 8, bg = "white")
    }, error = function(e) message("  Skipping RNA vln for ", i, ": ", e$message))
    
    # 3) Feature plot
    tryCatch({
      pal <- viridis(n = 10, option = "A")
      fea.pl <- FeaturePlot_scCustom(obj, reduction = "wnn.umap",
                                     features = i, colors_use = pal, order = TRUE)
      ggsave(file.path(rna_feature_dir, paste0(i, "_RNA_Featureplot_Magma.png")),
             plot = fea.pl, dpi = 500, width = 8, bg = "white")
    }, error = function(e) message("  Skipping RNA feature for ", i, ": ", e$message))
  }
  
  message(paste0("[", subset.label, "] Plot generation complete."))
}

# =============================================================================
# HELPER FUNCTION: Run subclustering (WNN-aware)
# =============================================================================
run_subclustering <- function(obj,
                              resolution = 0.3,   # start low, walk up
                              n.pcs.rna  = 15,    # fewer dims for a homogeneous subset
                              n.pcs.adt  = 10,
                              n.hvg      = 2000) {
  
  # ---- RNA: RECOMPUTE HVGs ON THE SUBSET, then scale + PCA -------------------
  # Critical: after subsetting, VariableFeatures() still holds the genes that
  # separated the MAJOR lineages (T vs B vs NK vs mono). A PCA built on those
  # captures between-lineage / contaminant structure rather than variation
  # *within* B cells -> spurious small clusters. Recompute on the subset.
  DefaultAssay(obj) <- "RNA"
  obj <- FindVariableFeatures(obj, nfeatures = n.hvg, verbose = FALSE)
  obj <- ScaleData(obj, verbose = TRUE)
  obj <- RunPCA(obj, npcs = 30, verbose = TRUE)
  
  # ---- ADT: scale + PCA on all proteins -------------------------------------
  DefaultAssay(obj) <- "ADT"
  VariableFeatures(obj) <- rownames(obj)
  obj <- ScaleData(obj, verbose = TRUE)
  obj <- RunPCA(obj, reduction.name = "apca",
                npcs = min(nrow(obj) - 1, 30), verbose = TRUE)
  
  # ---- WNN graph (fewer dims than the whole-dataset run) --------------------
  obj <- FindMultiModalNeighbors(
    obj,
    reduction.list = list("pca", "apca"),
    dims.list      = list(1:n.pcs.rna, 1:n.pcs.adt),
    modality.weight.name = "RNA.weight",
    verbose        = FALSE
  )
  
  # Re-run UMAP on WNN graph
  obj <- RunUMAP(obj, nn.name = "weighted.nn",
                 reduction.name = "wnn.umap",
                 reduction.key  = "wnnUMAP_",
                 verbose = TRUE)
  
  # Clustering
  obj <- FindClusters(obj, graph.name = "wsnn",
                      algorithm = 3,
                      resolution = resolution,
                      verbose = TRUE)
  
  # Store subcluster IDs cleanly
  obj$Subcluster_ID <- as.character(Idents(obj))
  return(obj)
}

# =============================================================================
# 1.  B-CELL SUBCLUSTER
# =============================================================================
message("===== B-cell Subclustering =====")

# -----------------------------------------------------------------------------
# B-cell clusters confirmed from table(Cluster_ID, Manual_Annotation):
#   3  = B Naïve (4902)
#   12 = Transitional B cells (2168)
#   21 = B intermediate (279)
# (No separate plasma/plasmablast cluster in this annotation.)
bcell.clusters <- c("3","12","21")
# -----------------------------------------------------------------------------

OPIS_BCELL <- subset(OPIS_ALL,
                     subset = Cluster_ID %in% bcell.clusters)

message(paste0("B-cell subset: ", ncol(OPIS_BCELL), " cells"))

# Run subclustering (tune these: start low and walk up)
OPIS_BCELL <- run_subclustering(OPIS_BCELL,
                                resolution = 0.3,
                                n.pcs.rna  = 15,
                                n.pcs.adt  = 10)

# UMAP coloured by subcluster (boxed SeuratExtend style)
umap.bcell <- DimPlot2(OPIS_BCELL, reduction = "wnn.umap",
                       group.by = "Subcluster_ID", label = TRUE,
                       repel = TRUE, box = TRUE, label.size = 5) +
  ggtitle("B-cell | Subclusters (res 0.3)")

ggsave(file.path(bcell.save.path, "plots/UMAP", "Bcell_Subclusters_UMAP.png"),
       plot = umap.bcell, dpi = 500, width = 10, height = 8, bg = "white")

# Also show original Cluster_ID for reference
umap.bcell.orig <- DimPlot2(OPIS_BCELL, reduction = "wnn.umap",
                            group.by = "Cluster_ID", label = TRUE,
                            repel = TRUE, box = TRUE, label.size = 5) +
  ggtitle("B-cell | Original Cluster IDs")

ggsave(file.path(bcell.save.path, "plots/UMAP", "Bcell_OriginalClusterIDs_UMAP.png"),
       plot = umap.bcell.orig, dpi = 500, width = 10, height = 8, bg = "white")

# Generate annotation plots
generate_annotation_plots(OPIS_BCELL, bcell.save.path, "Bcell")

# Save object to saved_R_data
qs_save(OPIS_BCELL,
        file = file.path(load.path, "OPIS_BCELL_PreAnnotation.qs2"))
message("B-cell object saved.")

# =============================================================================
# RNA_Extra plots — B-cell-focused panel
# =============================================================================

# ---------------------------- #
# 1) (Re)load object
# ---------------------------- #
opis_save_dir <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/saved_R_data"
OPIS_BCELL <- qs_read(file.path(opis_save_dir, "OPIS_BCELL_PreAnnotation.qs2"))

# ---------------------------- #
# 2) B-cell extra RNA panel
#    (subset-relevant markers not central to the broad panel above;
#     missing genes are skipped automatically)
# ---------------------------- #
rna.features <- c(
  # Pan-B
  'CD19','MS4A1','CD79A','CD79B','CD74','PAX5','EBF1',
  # Naive / transitional
  'IGHD','TCL1A','FCER2','IL4R','CD72','CD24','CD38','MME',
  # Memory / activated
  'CD27','CD80','CD86','TNFRSF13B','CD83','NR4A1','EGR1','MYC',
  # Germinal-center
  'BCL6','AICDA','RGS13','MEF2B','LRMP','MKI67','TOP2A',
  # Atypical / age-associated (ABC / DN2)
  'TBX21','ITGAX','FCRL4','FCRL5','ZEB2','FGR',
  # Marginal-zone-like
  'CD1C','CR2',
  # Transcription factors
  'IRF4','IRF8','SPIB','POU2AF1','PRDM1',
  # Plasmablast / plasma cell
  'XBP1','MZB1','JCHAIN','SDC1','SLAMF7','DERL3','TNFRSF17',
  # Isotypes
  'IGHM','IGHG1','IGHA1'
)

rna.features.2 <- character(0)

# ---------------------------- #
# 3) Output folders
# ---------------------------- #
bcell_plot_base <- file.path(bcell.save.path, "plots")

bcell_rna_extra_dir <- file.path(bcell_plot_base, "RNA_Extra")
dir.create(bcell_rna_extra_dir, recursive = TRUE, showWarnings = FALSE)

bcell_rna_vln_dir      <- file.path(bcell_rna_extra_dir, "Violin")
bcell_rna_vln_OUD_dir  <- file.path(bcell_rna_extra_dir, "Violin_split_by_OUD_status")
bcell_rna_feature_dir  <- file.path(bcell_rna_extra_dir, "FeaturePlot")

dir.create(bcell_rna_vln_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(bcell_rna_vln_OUD_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(bcell_rna_feature_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------- #
# 4) Plotting function
# ---------------------------- #
make_rna_extra_plots <- function(
    seu,
    feature_vec,
    out_vln_dir,
    out_vln_oud_dir,
    out_feature_dir,
    obj_name = "SeuratObject"
) {
  
  DefaultAssay(seu) <- "RNA"
  pal <- viridis(n = 10, option = "A")
  
  for (i in feature_vec) {
    if (i %in% rownames(seu[["RNA"]])) {
      
      message("[", obj_name, "] Plotting: ", i)
      
      # 1) SPLIT-BY OUD STATUS VIOLIN
      vln.pl.split <- VlnPlot2(
        seu,
        features    = i,
        cols        = "default",
        split.by    = "OUD_status",
        stat.method = "wilcox.test",
        show.mean   = TRUE,
        mean_colors = c("red", "blue")
      ) + ggtitle(paste("RNA |", i, "| Split by OUD_status"))
      
      ggsave(
        filename = file.path(out_vln_oud_dir, paste0(i, "_RNA_VLNplot_splitBy_OUD_status.png")),
        plot     = vln.pl.split,
        dpi      = 500,
        width    = 14,
        height   = 8,
        bg       = "white"
      )
      
      # 2) NON-SPLIT VIOLIN PLOT
      vln.pl <- VlnPlot2(
        seu,
        features    = i,
        cols        = "default",
        show.mean   = TRUE,
        mean_colors = c("red", "blue")
      ) + ggtitle(paste("RNA |", i))
      
      ggsave(
        filename = file.path(out_vln_dir, paste0(i, "_RNA_VLNplot.png")),
        plot     = vln.pl,
        dpi      = 500,
        width    = 14,
        height   = 8,
        bg       = "white"
      )
      
      # 3) FEATURE PLOT
      fea.pl <- FeaturePlot_scCustom(
        seu,
        reduction  = "wnn.umap",
        features   = i,
        colors_use = pal,
        order      = TRUE
      )
      
      ggsave(
        filename = file.path(out_feature_dir, paste0(i, "_RNA_Featureplot_Magma.png")),
        plot     = fea.pl,
        dpi      = 500,
        width    = 8,
        height   = 7,
        bg       = "white"
      )
    } else {
      message("[", obj_name, "] Skipping missing gene: ", i)
    }
  }
}

# ---------------------------- #
# 5) Run for B cells
# ---------------------------- #
make_rna_extra_plots(
  seu             = OPIS_BCELL,
  feature_vec     = rna.features,
  out_vln_dir     = bcell_rna_vln_dir,
  out_vln_oud_dir = bcell_rna_vln_OUD_dir,
  out_feature_dir = bcell_rna_feature_dir,
  obj_name        = "OPIS_BCELL"
)

if (length(rna.features.2) > 0) {
  make_rna_extra_plots(
    seu             = OPIS_BCELL,
    feature_vec     = rna.features.2,
    out_vln_dir     = bcell_rna_vln_dir,
    out_vln_oud_dir = bcell_rna_vln_OUD_dir,
    out_feature_dir = bcell_rna_feature_dir,
    obj_name        = "OPIS_BCELL"
  )
}

message("Done.")

# =============================================================================
# Marker tables for annotation (RNA + ADT) -> CSV
# =============================================================================
marker.dir <- bcell.save.path
Idents(OPIS_BCELL) <- "Subcluster_ID"

# ---- RNA markers -------------------------------------------------------------
message("FindAllMarkers: RNA ...")
DefaultAssay(OPIS_BCELL) <- "RNA"
bcell.markers.rna <- FindAllMarkers(OPIS_BCELL,
                                    only.pos        = TRUE,
                                    min.pct         = 0.25,
                                    logfc.threshold = 0.25)
write.csv(bcell.markers.rna,
          file.path(marker.dir, "Bcell_markers_RNA.csv"),
          row.names = FALSE)

# ---- ADT markers -------------------------------------------------------------
message("FindAllMarkers: ADT ...")
DefaultAssay(OPIS_BCELL) <- "ADT"
bcell.markers.adt <- FindAllMarkers(OPIS_BCELL,
                                    only.pos        = TRUE,
                                    min.pct         = 0.25,
                                    logfc.threshold = 0.10)  # ADT effects are smaller
write.csv(bcell.markers.adt,
          file.path(marker.dir, "Bcell_markers_ADT.csv"),
          row.names = FALSE)

message("Marker tables written to: ", marker.dir)

# =============================================================================
# Average-expression export for annotation (RNA + ADT)
#   RNA : AverageExpression on the "data" slot. RNA data is log1p-normalized,
#         so AverageExpression's internal expm1() back-transform is CORRECT.
#   ADT : MANUAL rowMeans on the "data" slot -- NO expm1 back-transformation.
#         ADT here is DSB-normalized (isotype-corrected), NOT log1p counts, so
#         expm1() would warp the values. A plain arithmetic mean in DSB space
#         is the right per-cluster summary.
# =============================================================================
avg.dir  <- file.path(bcell.save.path, "avg_expression")
dir.create(avg.dir, recursive = TRUE, showWarnings = FALSE)

group_by <- "Subcluster_ID"

# min-max scale each row (gene/protein) to 0-1 across clusters (for heatmaps)
scale_01 <- function(m) {
  m <- as.matrix(m)
  out <- t(apply(m, 1, function(x) {
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) == 0) return(rep(0, length(x)))
    (x - rng[1]) / (rng[2] - rng[1])
  }))
  colnames(out) <- colnames(m)
  out
}

clusters <- OPIS_BCELL@meta.data[[group_by]]
cl.levels <- levels(factor(clusters))

# ---- RNA (AverageExpression is fine for log-normalized RNA) ------------------
message("AverageExpression: RNA ...")
DefaultAssay(OPIS_BCELL) <- "RNA"
rna_avail <- intersect(rna.features, rownames(OPIS_BCELL[["RNA"]]))

avg_rna <- AverageExpression(OPIS_BCELL, assays = "RNA", features = rna_avail,
                             group.by = group_by, slot = "data")$RNA
colnames(avg_rna) <- gsub("^g ", "", colnames(avg_rna))
write.csv(avg_rna, file.path(avg.dir, "AvgExpr_RNA.csv"))
write.csv(round(scale_01(avg_rna), 4), file.path(avg.dir, "AvgExpr_RNA_scaled.csv"))

avg_rna_all <- AverageExpression(OPIS_BCELL, assays = "RNA",
                                 group.by = group_by, slot = "data")$RNA
colnames(avg_rna_all) <- gsub("^g ", "", colnames(avg_rna_all))
write.csv(avg_rna_all, file.path(avg.dir, "AvgExpr_RNA_AllGenes.csv"))

# ---- ADT (manual mean on DSB data slot; NO expm1) ----------------------------
message("AverageExpression: ADT (manual, DSB-safe) ...")
DefaultAssay(OPIS_BCELL) <- "ADT"
adt_avail <- rownames(OPIS_BCELL[["ADT"]])
adt_mat   <- GetAssayData(OPIS_BCELL, assay = "ADT", slot = "data")[adt_avail, , drop = FALSE]

avg_adt <- sapply(cl.levels, function(cl) {
  cells <- which(clusters == cl)
  if (length(cells) == 0) return(rep(0, length(adt_avail)))
  rowMeans(adt_mat[, cells, drop = FALSE])
})
rownames(avg_adt) <- adt_avail
write.csv(avg_adt, file.path(avg.dir, "AvgExpr_ADT.csv"))
write.csv(round(scale_01(avg_adt), 4), file.path(avg.dir, "AvgExpr_ADT_scaled.csv"))

# ---- Percent expression ------------------------------------------------------
# RNA: % of cells per cluster with raw counts > 0
rna_counts <- GetAssayData(OPIS_BCELL, assay = "RNA", slot = "counts")[rna_avail, , drop = FALSE]
pct_rna <- sapply(cl.levels, function(cl) {
  cells <- which(clusters == cl)
  if (length(cells) == 0) return(rep(0, length(rna_avail)))
  rowMeans(rna_counts[, cells, drop = FALSE] > 0) * 100
})
rownames(pct_rna) <- rna_avail
write.csv(round(pct_rna, 2), file.path(avg.dir, "PctExpr_RNA.csv"))

# ADT: % of cells per cluster above background (DSB > 0)
pct_adt <- sapply(cl.levels, function(cl) {
  cells <- which(clusters == cl)
  if (length(cells) == 0) return(rep(0, length(adt_avail)))
  rowMeans(adt_mat[, cells, drop = FALSE] > 0) * 100
})
rownames(pct_adt) <- adt_avail
write.csv(round(pct_adt, 2), file.path(avg.dir, "PctExpr_ADT.csv"))

message("Average-expression tables written to: ", avg.dir)
message("Done.")

# =============================================================================
# CLEAN + ANNOTATE B-cell subclusters
#   - Removes contaminant / doublet subclusters (NO re-clustering; keeps the
#     existing PCA / WNN / wnn.umap embedding -- cheap, just drops cells)
#   - Assigns a manual identity to each surviving subcluster
#   - Justification for every call is documented inline
#   Run AFTER subclustering + marker/avg-expr export, with OPIS_BCELL in memory.
# =============================================================================

# ---- 1. Drop contaminant / doublet subclusters -------------------------------
# These three are NOT B-cell states -- they are doublets / ambient contamination,
# and account for much of the apparent over-fragmentation:
#   6  T cells / B-T doublets    : CD3D/CD3E/CD8A/CD8B/IL7R/GNLY/NKG7 (RNA);
#                                  CD3/CD8/TCR-AB/CD5/CD2 (ADT); lowest MS4A1.
#   7  Platelet doublets         : PPBP/PF4/ITGA2B (RNA); CD41(ITGA2B)/CD62P(SELP)/
#                                  CD36 (ADT); had 0 significant RNA markers.
#   9  Monocyte/myeloid doublets : CD14/LYZ/S100A8/FCN1/CST3 (RNA);
#                                  CD33/CD11b(ITGAM)/CD64(FCGR1A) (ADT).
drop.subclusters <- c("6","7","9")
keep.subclusters <- setdiff(levels(factor(OPIS_BCELL$Subcluster_ID)), drop.subclusters)

OPIS_BCELL <- subset(OPIS_BCELL, subset = Subcluster_ID %in% keep.subclusters)
OPIS_BCELL$Subcluster_ID <- droplevels(factor(OPIS_BCELL$Subcluster_ID))
Idents(OPIS_BCELL) <- "Subcluster_ID"
message("After cleaning: ", ncol(OPIS_BCELL), " cells across subclusters ",
        paste(levels(OPIS_BCELL$Subcluster_ID), collapse = ", "))

# ---- 2. Manual annotation map (RNA + ADT evidence + reasoning) ---------------
#
#  0  "Naive B (IEG-high)"
#       Naive surface phenotype (IGHD+, CD23/FCER2+, CD27-, ADT CD1D+) overlaid
#       with a strong immediate-early / heat-shock program (FOS, FOSB, DUSP1,
#       DNAJB1, HSP90AA1). Reads as naive B carrying an activation / dissociation-
#       stress signature rather than a separate lineage.
#       [VERIFY: is it donor- or batch-restricted? if so it's technical.]
#
#  1  "Naive B (resting)"
#       Textbook resting naive: IGHD = max, FCER2/CD23 = max, SELL+, CD27-,
#       TCL1A+, no activation program. The core naive population.
#
#  2  "Memory B"
#       CD27 = max with TNFRSF13B (TACI), CD80, CD70, FAS; ADT CD27 hi, CD24 hi.
#       Classic CD27+ (switched/activated) memory phenotype.
#
#  3  "Intermediate B (HLA-G+)"
#       Intermediate CD27, FCRL4 hi, distinctive HLA-G / HLA-DQA2 program, very
#       high surface BAFF-R (TNFRSF13C, ADT). Unswitched/intermediate B with a
#       possible regulatory lean.
#       [VERIFY: IL10 / CD5 to assess a Breg interpretation.]
#
#  4  "Transitional B"
#       TCL1A = max, MZB1 = max, PAX5 = max, IGHD+; ADT CD24 / CD38 elevated.
#       TCL1A^hi MZB1^hi CD24^hi CD38^hi = canonical transitional signature.
#
#  5  "Atypical B (ABC/DN2)"
#       T-bet (TBX21), CD11c (ITGAX), FCRL5, ZEB2, FGR, TOX2 all = max; ADT
#       CD11c / CD86 / CD22 hi. Unambiguous age-associated / atypical B (DN2).
#
#  8  "Transitional/Naive (T2-like)"
#       TCL1A+, IGHD+, PAX5+, ROR1+; ADT CD38 hi, CD21(CR2) hi, IGHD hi. Sits
#       between transitional and naive -- a later (T2/T3) transitional or
#       early-naive stage.
#       [VERIFY: CD24/CD38 co-level vs cluster 4 to place it on the T1->naive axis.]

annotation.map <- c(
  "0" = "Naive B (IEG-high)",
  "1" = "Naive B (resting)",
  "2" = "Memory B",
  "3" = "Intermediate B (HLA-G+)",
  "4" = "Transitional B",
  "5" = "Atypical B (ABC/DN2)",
  "8" = "Transitional/Naive (T2-like)"
)

OPIS_BCELL$Bcell_Annotation <- unname(annotation.map[as.character(OPIS_BCELL$Subcluster_ID)])

# Order labels along a rough developmental axis for tidy plots / legends
annotation.levels <- c(
  "Transitional B",
  "Transitional/Naive (T2-like)",
  "Naive B (resting)",
  "Naive B (IEG-high)",
  "Intermediate B (HLA-G+)",
  "Memory B",
  "Atypical B (ABC/DN2)"
)
OPIS_BCELL$Bcell_Annotation <- factor(OPIS_BCELL$Bcell_Annotation, levels = annotation.levels)
Idents(OPIS_BCELL) <- "Bcell_Annotation"

# sanity check: every subcluster maps to exactly one label
print(table(OPIS_BCELL$Subcluster_ID, OPIS_BCELL$Bcell_Annotation))

# ---- 3. Re-plot on the EXISTING embedding (no re-clustering) ------------------
umap.annot <- DimPlot2(OPIS_BCELL, reduction = "wnn.umap",
                       group.by = "Bcell_Annotation", label = TRUE,
                       repel = TRUE, box = TRUE, label.size = 5) +
  ggtitle("B-cell | Annotated (doublets removed)")

ggsave(file.path(bcell.save.path, "plots/UMAP", "Bcell_Annotated_UMAP.png"),
       plot = umap.annot, dpi = 500, width = 11, height = 8, bg = "white")

# ---- 4. Save annotated object (new file; leaves PreAnnotation intact) --------
qs_save(OPIS_BCELL, file = file.path(load.path, "OPIS_BCELL_Annotated.qs2"))
message("Annotated B-cell object saved: ", file.path(load.path, "OPIS_BCELL_Annotated.qs2"))
