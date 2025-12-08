# ========================================= #
#   Cell Cycle Effect & Doublet Detection   #
# ========================================= #

# ------------ #
# Libraries    #
# ------------ #
library(Seurat)          # Main analysis
library(ggplot2)         # Plotting
library(DoubletFinder)   # Doublet detection (optional; scDblFinder used below)
library(scDblFinder)     # Doublet detection
library(Azimuth)         # Reference-based annotation
library(scCustomize)
library(qs2)             # qs_read / qs_save

# ---------------------------- #
# Paths & Output Directories   #
# ---------------------------- #

base_dir       <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq"
cell_cycle_dir <- file.path(base_dir, "QC", "Cell_Cycle")
doublet_dir    <- file.path(base_dir, "QC", "Doublets")
saved_dir      <- file.path(base_dir, "saved_R_data")

dir.create(cell_cycle_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(doublet_dir,    recursive = TRUE, showWarnings = FALSE)

# ---------------------------- #
# Load Filtered Seurat Object  #
# ---------------------------- #

filtered_seurat <- qs_read(file.path(saved_dir, "Seuratv5_filtered_seurat.qs2"))
seurat_clean <- filtered_seurat

# ---------------------------- #
# Cell Cycle Scoring           #
# ---------------------------- #
DefaultAssay(seurat_clean) <- "RNA"

# Azimuth annotation using PBMC reference
seurat_phase <- RunAzimuth(seurat_clean, reference = "pbmcref")

# Basic Transformations
seurat_phase <- NormalizeData(seurat_phase)
seurat_phase <- FindVariableFeatures(seurat_phase, selection.method = "vst")
seurat_phase <- ScaleData(seurat_phase)
seurat_phase <- RunPCA(
  seurat_phase,
  features        = VariableFeatures(seurat_phase),
  ndims.print     = 6:10,
  nfeatures.print = 10
)
seurat_phase <- FindNeighbors(seurat_phase, dims = 1:30, reduction = "pca")
seurat_phase <- FindClusters(seurat_phase, resolution = 2, cluster.name = "unintegrated_clusters")
seurat_phase <- RunUMAP(
  seurat_phase,
  dims           = 1:30,
  reduction      = "pca",
  reduction.name = "umap.unintegrated"
)

# Get Seurat’s built-in cell cycle gene lists
s.genes   <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

# Score cells for cell cycle phase
seurat_phase <- CellCycleScoring(
  seurat_phase,
  s.features  = s.genes,
  g2m.features = g2m.genes,
  set.ident   = TRUE
)

# ---------------------------- #
# Visualize Cell Cycle Effect  #
# ---------------------------- #

# Ridge plots for key markers
png(file.path(cell_cycle_dir, "RidgePlot_CellCycleMarkers.png"), width = 1800, height = 1200)
RidgePlot(seurat_phase, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
dev.off()

# Violin plot of S.Score & G2M.Score by sample (uses orig.ident here)
png(file.path(cell_cycle_dir, "CellCycle_Scores_bySample.png"), width = 1800, height = 1200)
VlnPlot(
  seurat_phase,
  features = c("S.Score", "G2M.Score"),
  group.by = "orig.ident",
  pt.size  = 0.1
)
dev.off()

# Violin plot of S.Score & G2M.Score by OUD status
if ("OUD_status" %in% colnames(seurat_phase@meta.data)) {
  png(file.path(cell_cycle_dir, "CellCycle_Scores_byOUD.png"), width = 1800, height = 1200)
  VlnPlot(
    seurat_phase,
    features = c("S.Score", "G2M.Score"),
    group.by = "OUD_status",
    pt.size  = 0.1
  )
  dev.off()
}

# PCA colored by Phase
p1 <- DimPlot_scCustom(
  seurat_phase,
  reduction = "pca",
  group.by  = "Phase",
  label     = TRUE,
  repel     = TRUE,
  label.box = TRUE,
  label.size = 3.5,
  pt.size    = 1
)
ggsave(
  filename = file.path(cell_cycle_dir, "PCA_by_CellCyclePhase.png"),
  plot     = p1,
  width    = 13,
  height   = 8,
  dpi      = 300
)

# Azimuth annotated over UMAP (unintegrated)
p2 <- DimPlot_scCustom(
  seurat_phase,
  reduction = "umap.unintegrated",
  group.by  = "predicted.celltype.l2",
  label     = TRUE,
  repel     = TRUE,
  label.box = TRUE,
  label.size = 3.5,
  pt.size    = 1
)
ggsave(
  filename = file.path(cell_cycle_dir, "Azimuth_Umap_Unintegrated.png"),
  plot     = p2,
  width    = 13,
  height   = 8,
  dpi      = 300
)

# Phase UMAP (unintegrated)
p3 <- DimPlot_scCustom(
  seurat_phase,
  reduction = "umap.unintegrated",
  group.by  = "Phase",
  label     = TRUE,
  repel     = TRUE,
  label.box = TRUE,
  label.size = 3.5,
  pt.size    = 1
)
ggsave(
  filename = file.path(cell_cycle_dir, "Cell_Cycle_Umap_Unintegrated.png"),
  plot     = p3,
  width    = 13,
  height   = 8,
  dpi      = 300
)

# Phase UMAP split by OUD status (instead of IGRA)
if ("OUD_status" %in% colnames(seurat_phase@meta.data)) {
  p4 <- DimPlot_scCustom(
    seurat_phase,
    reduction = "umap.unintegrated",
    group.by  = "Phase",
    split.by  = "OUD_status",
    label     = TRUE,
    repel     = TRUE,
    label.box = TRUE,
    label.size = 3.5,
    pt.size    = 1
  )
  ggsave(
    filename = file.path(cell_cycle_dir, "Cell_Cycle_OUD_Umap_Unintegrated.png"),
    plot     = p4,
    width    = 13,
    height   = 8,
    dpi      = 300
  )
}

# ---------------------------- #
# Doublet Detection (scDblFinder)
# ---------------------------- #

# Split object by sample for individual doublet calling
# For this dataset we'll use orig.ident as the "sample" key
split_seurat <- SplitObject(seurat_phase, split.by = "orig.ident")
samples      <- names(split_seurat)

for (i in samples) {
  # Run scDblFinder on raw counts matrix
  counts_mat <- GetAssayData(split_seurat[[i]], slot = "counts")
  sce        <- scDblFinder(counts_mat)
  
  # Store doublet info in Seurat object
  split_seurat[[i]]$scDblFinder.score <- sce$scDblFinder.score
  split_seurat[[i]]$scDblFinder.class <- sce$scDblFinder.class
  
  # UMAP by doublet classification
  p_doublet <- DimPlot_scCustom(
    split_seurat[[i]],
    reduction = "umap.unintegrated",
    group.by  = "scDblFinder.class",
    colors_use = c("#1f78b4", "#e31a1c"),
    label     = TRUE,
    repel     = TRUE,
    label.box = TRUE,
    label.size = 3.5,
    pt.size    = 1
  ) + ggtitle(paste0(i, " Doublets"))
  
  ggsave(
    filename = file.path(doublet_dir, paste0(i, "_Doublets.png")),
    plot     = p_doublet,
    width    = 13,
    height   = 8,
    dpi      = 300
  )
}

# ---------------------------- #
# Remove Doublets              #
# ---------------------------- #
for (i in samples) {
  split_seurat[[i]] <- subset(split_seurat[[i]], subset = scDblFinder.class == "singlet")
}

# Merge back post-doublet removal
seurat_phase_clean <- merge(split_seurat[[1]], y = split_seurat[-1])

# ---------------------------- #
# Save Output                  #
# ---------------------------- #
qs_save(
  seurat_phase_clean,
  file = file.path(saved_dir, "Seuratv5_filtered_CellCycle_DoubletClean.qs2")
)
