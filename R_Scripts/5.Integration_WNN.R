## ============================================
##   OPIS_ALL RNA + ADT Integration + WNN
##   OPIS_ECCITEseq (qs2 + SeuratExtend)
## ============================================

## --------------------------
## 0) Libraries
## --------------------------
library(Seurat)
library(SingleCellExperiment)
library(scater)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(hdf5r)
library(dsb)
library(data.table)
library(ggplot2)
library(gridExtra)
library(SeuratWrappers)
library(Azimuth)
library(ggrepel)
library(patchwork)
library(SeuratExtend)   # DimPlot2
library(circlize)
library(ComplexHeatmap)
library(readxl)
library(batchelor)
library(harmony)
library(reticulate)
library(qs2)            # qs_read / qs_save
library(clustree)
## --------------------------
## 1) Paths & folders
## --------------------------
base_dir   <- "~/Documents/OPIS_ECCITEseq"
load.path  <- file.path(base_dir, "saved_R_data")

integration_dir <- file.path(base_dir, "Integration")
plot_dir        <- file.path(base_dir, "Plots", "Integration")

dir.create(integration_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir,        recursive = TRUE, showWarnings = FALSE)

## --------------------------
## 2) Load Seurat object (qs2)
## --------------------------
# File: Seurat_isotype.qs2
seurat_isotype <- qs_read(file.path(load.path, "Seurat_isotype.qs2"))

## --------------------------
## 3) Prep RNA & define OPIS_ALL
## --------------------------
# Join RNA layers
seurat_isotype[["RNA"]] <- JoinLayers(seurat_isotype[["RNA"]])

# If Cohort column exists and has "OPIS", subset; otherwise use full object
if ("Cohort" %in% colnames(seurat_isotype@meta.data) &&
    "OPIS" %in% seurat_isotype$Cohort) {
  OPIS_ALL <- subset(seurat_isotype, subset = Cohort == "OPIS")
} else {
  OPIS_ALL <- seurat_isotype
}

rm(seurat_isotype)
gc()

## Split RNA layers per batch (orig.ident)
OPIS_ALL[["RNA"]] <- split(OPIS_ALL[["RNA"]], f = OPIS_ALL$orig.ident)

## --------------------------
## 4) RNA — unintegrated + integrated (CCA + FastMNN)
## --------------------------
DefaultAssay(OPIS_ALL) <- "RNA"

OPIS_ALL <- NormalizeData(OPIS_ALL)
OPIS_ALL <- FindVariableFeatures(OPIS_ALL)
OPIS_ALL <- ScaleData(OPIS_ALL)
OPIS_ALL <- RunPCA(OPIS_ALL)

OPIS_ALL <- FindNeighbors(OPIS_ALL, dims = 1:30, reduction = "pca")
OPIS_ALL <- FindClusters(OPIS_ALL, resolution = 2, cluster.name = "unintegrated_clusters")
OPIS_ALL <- RunUMAP(
  OPIS_ALL,
  dims           = 1:30,
  reduction      = "pca",
  reduction.name = "umap.unintegrated"
)

# CCA integration
OPIS_ALL <- IntegrateLayers(
  OPIS_ALL,
  method         = CCAIntegration,
  orig.reduction = "pca",
  assay          = "RNA",
  new.reduction  = "integrated.cca.rna"
)

# FastMNN integration
OPIS_ALL <- IntegrateLayers(
  OPIS_ALL,
  method        = FastMNNIntegration,
  assay         = "RNA",
  new.reduction = "integrated.mnn.rna"
)

# UMAP / clustering on integrated reductions
# CCA
OPIS_ALL <- FindNeighbors(OPIS_ALL, reduction = "integrated.cca.rna", dims = 1:30)
OPIS_ALL <- FindClusters(OPIS_ALL, resolution = 2, cluster.name = "cca_clusters_rna")
OPIS_ALL <- RunUMAP(
  OPIS_ALL,
  reduction      = "integrated.cca.rna",
  dims           = 1:30,
  reduction.name = "umap.cca.rna"
)

# FastMNN
OPIS_ALL <- FindNeighbors(OPIS_ALL, reduction = "integrated.mnn.rna", dims = 1:30)
OPIS_ALL <- FindClusters(OPIS_ALL, resolution = 2, cluster.name = "mnn_clusters_rna")
OPIS_ALL <- RunUMAP(
  OPIS_ALL,
  reduction      = "integrated.mnn.rna",
  dims           = 1:30,
  reduction.name = "umap.mnn.rna"
)

## --------------------------
## 5) RNA Integration Plots (DimPlot2)
## --------------------------
p_rna_cca_celltype <- DimPlot2(
  OPIS_ALL,
  reduction  = "umap.cca.rna",
  group.by   = "predicted.celltype.l2",
  label      = TRUE,
  label.size = 2
) + ggtitle("OPIS_ALL RNA CCA – Azimuth L2")

p_rna_cca_cluster <- DimPlot2(
  OPIS_ALL,
  reduction  = "umap.cca.rna",
  group.by   = "cca_clusters_rna",
  label      = TRUE,
  label.size = 2
) + ggtitle("OPIS_ALL RNA CCA – CCA clusters")

p_rna_mnn_celltype <- DimPlot2(
  OPIS_ALL,
  reduction  = "umap.mnn.rna",
  group.by   = "predicted.celltype.l2",
  label      = TRUE,
  label.size = 2
) + ggtitle("OPIS_ALL RNA MNN – Azimuth L2")

p_rna_mnn_cluster <- DimPlot2(
  OPIS_ALL,
  reduction  = "umap.mnn.rna",
  group.by   = "mnn_clusters_rna",
  label      = TRUE,
  label.size = 2
) + ggtitle("OPIS_ALL RNA MNN – MNN clusters")

combined_rna <- wrap_plots(
  p_rna_cca_celltype,
  p_rna_cca_cluster,
  p_rna_mnn_celltype,
  p_rna_mnn_cluster,
  ncol = 2
)

ggsave(
  filename = file.path(plot_dir, "OPIS_ALL_RNA_Integration.png"),
  plot     = combined_rna,
  width    = 18,
  height   = 12,
  dpi      = 300
)

## --------------------------
## 6) ADT Integration
## --------------------------
DefaultAssay(OPIS_ALL) <- "ADT"

OPIS_ALL[["ADT"]] <- as(object = OPIS_ALL[["ADT"]], Class = "Assay5")
OPIS_ALL[["ADT"]] <- split(OPIS_ALL[["ADT"]], f = OPIS_ALL$orig.ident)

prots <- rownames(OPIS_ALL[["ADT"]]$data)
isotype_genes <- c("Mouse-IgG1", "Mouse-IgG2a", "Mouse-IgG2b",
                   "Rat-IgG2b", "Rat-IgG1", "Rat-IgG2a", "Hamster-IgG")
prots <- setdiff(prots, isotype_genes)

VariableFeatures(OPIS_ALL) <- prots

OPIS_ALL <- ScaleData(OPIS_ALL)
OPIS_ALL <- RunPCA(OPIS_ALL, reduction.name = "apca")
OPIS_ALL <- FindNeighbors(OPIS_ALL, dims = 1:30, reduction = "apca")
OPIS_ALL <- FindClusters(OPIS_ALL, resolution = 2, cluster.name = "unintegrated_clusters_adt")
OPIS_ALL <- RunUMAP(
  OPIS_ALL,
  dims           = 1:30,
  reduction      = "apca",
  reduction.name = "umap.unintegrated.adt"
)

# FastMNN integration on ADT
OPIS_ALL <- IntegrateLayers(
  OPIS_ALL,
  orig.reduction = "apca",
  features       = prots,
  method         = FastMNNIntegration,
  assay          = "ADT",
  new.reduction  = "integrated.mnn.adt"
)

# CCA integration on ADT
OPIS_ALL <- IntegrateLayers(
  OPIS_ALL,
  method         = CCAIntegration,
  orig.reduction = "apca",
  features       = prots,
  assay          = "ADT",
  new.reduction  = "integrated.cca.adt"
)

# UMAP / clustering on integrated ADT
# CCA
OPIS_ALL <- FindNeighbors(OPIS_ALL, reduction = "integrated.cca.adt", dims = 1:30)
OPIS_ALL <- FindClusters(OPIS_ALL, resolution = 2, cluster.name = "cca_clusters_adt")
OPIS_ALL <- RunUMAP(
  OPIS_ALL,
  reduction      = "integrated.cca.adt",
  dims           = 1:30,
  reduction.name = "umap.cca.adt"
)

# FastMNN
OPIS_ALL <- FindNeighbors(OPIS_ALL, reduction = "integrated.mnn.adt", dims = 1:30)
OPIS_ALL <- FindClusters(OPIS_ALL, resolution = 2, cluster.name = "mnn_clusters_adt")
OPIS_ALL <- RunUMAP(
  OPIS_ALL,
  reduction      = "integrated.mnn.adt",
  dims           = 1:30,
  reduction.name = "umap.mnn.adt"
)

## --------------------------
## 7) ADT Integration Plots (DimPlot2)
## --------------------------
p_adt_cca_celltype <- DimPlot2(
  OPIS_ALL,
  reduction  = "umap.cca.adt",
  group.by   = "predicted.celltype.l2",
  label      = TRUE,
  label.size = 2
) + ggtitle("OPIS_ALL ADT CCA – Azimuth L2")

p_adt_cca_cluster <- DimPlot2(
  OPIS_ALL,
  reduction  = "umap.cca.adt",
  group.by   = "cca_clusters_adt",
  label      = TRUE,
  label.size = 2
) + ggtitle("OPIS_ALL ADT CCA – CCA clusters")

p_adt_mnn_celltype <- DimPlot2(
  OPIS_ALL,
  reduction  = "umap.mnn.adt",
  group.by   = "predicted.celltype.l2",
  label      = TRUE,
  label.size = 2
) + ggtitle("OPIS_ALL ADT MNN – Azimuth L2")

p_adt_mnn_cluster <- DimPlot2(
  OPIS_ALL,
  reduction  = "umap.mnn.adt",
  group.by   = "mnn_clusters_adt",
  label      = TRUE,
  label.size = 2
) + ggtitle("OPIS_ALL ADT MNN – MNN clusters")

combined_adt <- wrap_plots(
  p_adt_cca_celltype,
  p_adt_cca_cluster,
  p_adt_mnn_celltype,
  p_adt_mnn_cluster,
  ncol = 2
)

ggsave(
  filename = file.path(plot_dir, "OPIS_ALL_ADT_Integration.png"),
  plot     = combined_adt,
  width    = 18,
  height   = 12,
  dpi      = 300
)

## --------------------------
## 8) Save integrated RNA+ADT object
## --------------------------
qs_save(
  OPIS_ALL,
  file = file.path(load.path, "OPIS_ALL_RNA_ADT.qs2")
)
OPIS_ALL <- qs_read(file.path(load.path, "OPIS_ALL_RNA_ADT.qs2"))
## --------------------------
## 9) WNN neighbors
## --------------------------
OPIS_ALL <- FindMultiModalNeighbors(
  OPIS_ALL,
  reduction.list       = list("integrated.mnn.rna", "integrated.mnn.adt"),
  dims.list            = list(1:30, 1:20),
  modality.weight.name = "wnn.weight"
)

OPIS_ALL <- RunUMAP(
  OPIS_ALL,
  nn.name        = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key  = "wnnUMAP_"
)

OPIS_ALL <- FindClusters(OPIS_ALL, graph.name = "wsnn", algorithm = 2, resolution = 0.5, cluster.name = "snn.louvianmlr_0.5")
OPIS_ALL <- FindClusters(OPIS_ALL, graph.name = "wsnn", algorithm = 2, resolution = 1,   cluster.name = "snn.louvianmlr_1")
OPIS_ALL <- FindClusters(OPIS_ALL, graph.name = "wsnn", algorithm = 2, resolution = 1.5, cluster.name = "snn.louvianmlr_1.5")

OPIS_ALL <- FindClusters(OPIS_ALL, graph.name = "wsnn", algorithm = 3, resolution = 0.5, cluster.name = "snn.slm_0.5")
OPIS_ALL <- FindClusters(OPIS_ALL, graph.name = "wsnn", algorithm = 3, resolution = 1,   cluster.name = "snn.slm_1")
OPIS_ALL <- FindClusters(OPIS_ALL, graph.name = "wsnn", algorithm = 3, resolution = 1.5, cluster.name = "snn.slm_1.5")

## --------------------------
## 10) Order cluster columns
## --------------------------
cluster_cols <- c(
  "snn.louvianmlr_0.5", "snn.louvianmlr_1", "snn.louvianmlr_1.5",
  "snn.slm_0.5", "snn.slm_1", "snn.slm_1.5"
)

for (cluster_col in cluster_cols) {
  if (cluster_col %in% colnames(OPIS_ALL@meta.data)) {
    vals <- OPIS_ALL[[cluster_col]][, 1]
    levs_num <- suppressWarnings(as.numeric(levels(vals)))
    if (!any(is.na(levs_num))) {
      OPIS_ALL[[cluster_col]][, 1] <- factor(vals, levels = sort(levs_num))
    } else {
      OPIS_ALL[[cluster_col]][, 1] <- factor(vals)
    }
  }
}

## --------------------------
## 11) WNN plots (DimPlot2)
## --------------------------
p1 <- DimPlot2(OPIS_ALL, reduction = "wnn.umap", group.by = "predicted.celltype.l1",
               label = TRUE, label.size = 5) + ggtitle("Azimuth L1")

p2 <- DimPlot2(OPIS_ALL, reduction = "wnn.umap", group.by = "predicted.celltype.l2",
               label = TRUE, label.size = 5) + ggtitle("Azimuth L2")

p3 <- DimPlot2(OPIS_ALL, reduction = "wnn.umap", group.by = "predicted.celltype.l3",
               label = TRUE, label.size = 5) + ggtitle("Azimuth L3")

p4 <- DimPlot2(OPIS_ALL, reduction = "wnn.umap", group.by = "snn.slm_0.5",
               label = TRUE, label.size = 5) + ggtitle("snn.slm_0.5")

p5 <- DimPlot2(OPIS_ALL, reduction = "wnn.umap", group.by = "snn.slm_1",
               label = TRUE, label.size = 5) + ggtitle("snn.slm_1")

p6 <- DimPlot2(OPIS_ALL, reduction = "wnn.umap", group.by = "snn.slm_1.5",
               label = TRUE, label.size = 5) + ggtitle("snn.slm_1.5")

p7 <- DimPlot2(OPIS_ALL, reduction = "wnn.umap", group.by = "snn.louvianmlr_0.5",
               label = TRUE, label.size = 5) + ggtitle("snn.louvianmlr_0.5")

p8 <- DimPlot2(OPIS_ALL, reduction = "wnn.umap", group.by = "snn.louvianmlr_1",
               label = TRUE, label.size = 5) + ggtitle("snn.louvianmlr_1")

p9 <- DimPlot2(OPIS_ALL, reduction = "wnn.umap", group.by = "snn.louvianmlr_1.5",
               label = TRUE, label.size = 5) + ggtitle("snn.louvianmlr_1.5")

combined_wnn <- wrap_plots(p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol = 3, nrow = 3)

ggsave(
  filename = file.path(plot_dir, "OPIS_ALL_WNN_Clusters.png"),
  plot     = combined_wnn,
  width    = 24,
  height   = 17,
  dpi      = 300
)
############ cLUSTREE cLUSTER SELECTION

# Define the range of resolutions you want to test
resolutions <- seq(0.2, 2.0, by = 0.2)

for (res in resolutions) {
  OPIS_ALL <- FindClusters(
    OPIS_ALL,
    graph.name = "wsnn",
    resolution = res,
    algorithm = 3  # Leiden (recommended)
  )
}

for (col in grep("^wsnn_res\\.", colnames(OPIS_ALL@meta.data), value = TRUE)) {
  lvls <- as.character(sort(as.numeric(levels(OPIS_ALL[[col]][, 1]))))
  OPIS_ALL[[col]] <- factor(OPIS_ALL[[col]][, 1], levels = lvls)
}

setwd('/home/akshay-iyer/Documents/OPIS_ECCITEseq/Integration/Clustree')



clustree(OPIS_ALL, prefix = "wsnn_res.")


ggsave('clustree.png',  width = 15,  # Adjust width as needed
       height = 9,  # Adjust height as needed
       dpi = 300,
       bg='white')


DimPlot2(
  OPIS_ALL,
  reduction = "wnn.umap",
  group.by = "wsnn_res.0.4",
  cols = 'light',
  label = TRUE, box = TRUE, repel = TRUE,
  label.color = "black", pt.size = 1
)

ggsave('cd8nk_subcluster_res0.4.png',
       dpi = 300,
       bg='white')



## --------------------------
## 12) Save final WNN object
## --------------------------
qs_save(
  OPIS_ALL,
  file = file.path(load.path, "OPIS_ALL_WNN.qs2")
)
