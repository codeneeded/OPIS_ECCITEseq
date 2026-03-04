# =============================================================================
# OPIS ECCITEseq - Subclustering: CD4 & NK/CD8
# Pre-Annotation Pipeline
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

cd4.save.path    <- file.path(subclust.root, "CD4")
nkcd8.save.path  <- file.path(subclust.root, "NK_CD8")

# Create all output subdirectories
for (base in c(cd4.save.path, nkcd8.save.path)) {
  for (sub in c("plots/ADT/violin", "plots/ADT/violin_OUD",
                "plots/ADT/feature",
                "plots/RNA/violin", "plots/RNA/violin_OUD",
                "plots/RNA/feature",
                "plots/UMAP",
                "qs2_objects")) {
    dir.create(file.path(base, sub), recursive = TRUE, showWarnings = FALSE)
  }
}

# ---- Load full Seurat object -------------------------------------------------
message("Loading OPIS_ALL ...")
OPIS_ALL <- qs_read(file.path(load.path, "OPIS_ALL_PostAnnotation.qs2"))

# ---- Feature lists -----------------------------------------------------------
rna.features <- c(
  'ASCL2','BATF','BATF3','BCL6','BACH2','C1QBP','CCL2','CCL3','CCL4L2','CCL5','CCND3','CD14','CD19','CD1C',
  'CD200','CD27','CD3D','CD3E','CD36','CD4','CD40','CD40LG','CD70','CD7','CD79A','CD8A','CD8B', "CCR5",
  'CLEC9A','CR2','CTLA4','CTSW','CXCL8','CXCR3','CXCR5','EBI3','ENTPD1','FABP5','FCGR2B','FCGR3A',
  'FCRL5','FOXP3','GNLY','GP1BA','GP9','GATA3','GZMK','HAVCR2','HIF1A','HIST1H4C','HLA-DPA1',
  'HLA-DRA','HLA-DRB1','ICOS','IFI30','IFNG','IGFBP2','IGFBP4','IGHA1','IGHA2','IGHG1','IGHG2','IGHG3','IGHG4','IGHM','IKZF2','IL10',
  'IL17A','IL18BP','IL18RAP','IL1B','IL21','IL2RA','IL2RB','IRF4','IRF8','ITGAX','JCHAIN','KLRB1',
  'KLRD1','KLRC1','LAG3','LDHA','LGALS1','LTA','LTB','MAF','MAL','MALAT1','MIR155HG','MKI67',
  'MT-ND1','MT-ND5','MS4A1','NELL2','NCAM1','NKG7','NR4A1','PDCD1','PF4','PPBP','PRDM1','PRF1',
  'RORC','SELL','SERPINA1','SERPING1','SH2D1A','TCF4','TCF7','TIGIT','TNF','TNFAIP2','TNFRSF18',
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
run_subclustering <- function(obj, resolution = 0.4) {
  # Re-run PCA on RNA
  DefaultAssay(obj) <- "RNA"
  obj <- ScaleData(obj, verbose = TRUE)
  obj <- RunPCA(obj, verbose = TRUE)
  
  # Re-run LSI / SVD on ADT
  DefaultAssay(obj) <- "ADT"
  obj <- ScaleData(obj, verbose = TRUE)
  obj <- RunPCA(obj, reduction.name = "apca", verbose = TRUE)
  
  # WNN graph
  obj <- FindMultiModalNeighbors(
    obj,
    reduction.list = list("pca", "apca"),
    dims.list      = list(1:30, 1:18),
    modality.weight.name = "RNA.weight",
    verbose        = FALSE
  )
  
  # Re-run UMAP on WNN graph (keeps existing wnn.umap if you prefer — swap names)
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
# 1.  CD4 SUBCLUSTER  (clusters 5, 6, 7, 8)
# =============================================================================
message("===== CD4 Subclustering =====")

cd4.clusters <- c("5","6","7","8")
OPIS_CD4 <- subset(OPIS_ALL,
                   subset = Cluster_ID %in% cd4.clusters)

message(paste0("CD4 subset: ", ncol(OPIS_CD4), " cells"))

# Run subclustering
OPIS_CD4 <- run_subclustering(OPIS_CD4, resolution = 0.5)

# UMAP coloured by subcluster
umap.cd4 <- DimPlot(OPIS_CD4, reduction = "wnn.umap",
                    group.by = "Subcluster_ID", label = TRUE,
                    repel = TRUE, pt.size = 0.5) +
  ggtitle("CD4 Subclusters (res 0.5)") +
  theme_cowplot()

ggsave(file.path(cd4.save.path, "plots/UMAP", "CD4_Subclusters_UMAP.png"),
       plot = umap.cd4, dpi = 500, width = 10, height = 8, bg = "white")

# Also show original Cluster_ID for reference
umap.cd4.orig <- DimPlot(OPIS_CD4, reduction = "wnn.umap",
                         group.by = "Cluster_ID", label = TRUE,
                         repel = TRUE, pt.size = 0.5) +
  ggtitle("CD4 — Original Cluster IDs") +
  theme_cowplot()

ggsave(file.path(cd4.save.path, "plots/UMAP", "CD4_OriginalClusterIDs_UMAP.png"),
       plot = umap.cd4.orig, dpi = 500, width = 10, height = 8, bg = "white")

# Generate annotation plots
generate_annotation_plots(OPIS_CD4, cd4.save.path, "CD4")

# Save object
qs_save(OPIS_CD4,
        file = file.path(cd4.save.path, "qs2_objects", "OPIS_CD4_PreAnnotation.qs2"))
message("CD4 object saved.")

# =============================================================================
# 2.  NK / CD8 SUBCLUSTER  (clusters 0,2,4,9,10,11,13,14,16,18,19)
# =============================================================================
message("===== NK/CD8 Subclustering =====")

nkcd8.clusters <- c("19","11","14","13","18","10","2","0","9","4","16")
OPIS_NKCD8 <- subset(OPIS_ALL,
                     subset = Cluster_ID %in% nkcd8.clusters)

message(paste0("NK/CD8 subset: ", ncol(OPIS_NKCD8), " cells"))

# Run subclustering
OPIS_NKCD8 <- run_subclustering(OPIS_NKCD8, resolution = 0.5)

# UMAP coloured by subcluster
umap.nkcd8 <- DimPlot(OPIS_NKCD8, reduction = "wnn.umap",
                      group.by = "Subcluster_ID", label = TRUE,
                      repel = TRUE, pt.size = 0.5) +
  ggtitle("NK/CD8 Subclusters (res 0.5)") +
  theme_cowplot()

ggsave(file.path(nkcd8.save.path, "plots/UMAP", "NKCD8_Subclusters_UMAP.png"),
       plot = umap.nkcd8, dpi = 500, width = 10, height = 8, bg = "white")

umap.nkcd8.orig <- DimPlot(OPIS_NKCD8, reduction = "wnn.umap",
                           group.by = "Cluster_ID", label = TRUE,
                           repel = TRUE, pt.size = 0.5) +
  ggtitle("NK/CD8 — Original Cluster IDs") +
  theme_cowplot()

ggsave(file.path(nkcd8.save.path, "plots/UMAP", "NKCD8_OriginalClusterIDs_UMAP.png"),
       plot = umap.nkcd8.orig, dpi = 500, width = 10, height = 8, bg = "white")

# Generate annotation plots
generate_annotation_plots(OPIS_NKCD8, nkcd8.save.path, "NK_CD8")

# Save object
qs_save(OPIS_NKCD8,
        file = file.path(nkcd8.save.path, "qs2_objects", "OPIS_NKCD8_PreAnnotation.qs2"))
message("NK/CD8 object saved.")

# =============================================================================
# Resave Subclustering UMAPs using SeuratExtend DimPlot2
# =============================================================================

cd4.umap.dir   <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/subclustering/pre_annotation/CD4/plots/UMAP"
nkcd8.umap.dir <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/subclustering/pre_annotation/NK_CD8/plots/UMAP"

# =============================================================================
# CD4
# =============================================================================

# Subcluster IDs
p1 <- DimPlot2(
  OPIS_CD4,
  reduction  = "wnn.umap",
  group.by   = "Subcluster_ID",
  label      = TRUE,
  repel      = TRUE,
  box        = TRUE,
  label.size = 5
) + ggtitle("CD4 | Subclusters (res 0.5)")

ggsave(file.path(cd4.umap.dir, "CD4_Subclusters_UMAP.png"),
       plot = p1, dpi = 500, width = 10, height = 8, bg = "white")

# Original Cluster IDs
p2 <- DimPlot2(
  OPIS_CD4,
  reduction  = "wnn.umap",
  group.by   = "Cluster_ID",
  label      = TRUE,
  repel      = TRUE,
  box        = TRUE,
  label.size = 5
) + ggtitle("CD4 | Original Cluster IDs")

ggsave(file.path(cd4.umap.dir, "CD4_OriginalClusterIDs_UMAP.png"),
       plot = p2, dpi = 500, width = 10, height = 8, bg = "white")

# =============================================================================
# NK/CD8
# =============================================================================

# Subcluster IDs
p3 <- DimPlot2(
  OPIS_NKCD8,
  reduction  = "wnn.umap",
  group.by   = "Subcluster_ID",
  label      = TRUE,
  repel      = TRUE,
  box        = TRUE,
  label.size = 5
) + ggtitle("NK/CD8 | Subclusters (res 0.5)")

ggsave(file.path(nkcd8.umap.dir, "NKCD8_Subclusters_UMAP.png"),
       plot = p3, dpi = 500, width = 10, height = 8, bg = "white")

# Original Cluster IDs
p4 <- DimPlot2(
  OPIS_NKCD8,
  reduction  = "wnn.umap",
  group.by   = "Cluster_ID",
  label      = TRUE,
  repel      = TRUE,
  box        = TRUE,
  label.size = 5
) + ggtitle("NK/CD8 | Original Cluster IDs")

ggsave(file.path(nkcd8.umap.dir, "NKCD8_OriginalClusterIDs_UMAP.png"),
       plot = p4, dpi = 500, width = 10, height = 8, bg = "white")

message("DimPlot2 UMAPs saved.")

