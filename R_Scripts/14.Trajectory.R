# =============================================================================
# OPIS ECCITEseq - Subclustering: CD4 & NK/CD8
# Post-Annotation Pipeline
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
load.path     <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/saved_R_data"
subclust.root <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/subclustering/post_annotation"

dir.create(subclust.root, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1) Load Objects
# =============================================================================

opis_save_dir <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/saved_R_data"

OPIS_NKCD8 <- qs_read(file.path(opis_save_dir, "OPIS_NKCD8_PreAnnotation.qs2"))
OPIS_CD4   <- qs_read(file.path(opis_save_dir, "OPIS_CD4_PreAnnotation.qs2"))

# =============================================================================
# 2) Annotate NK/CD8 Object
# =============================================================================
# Cluster-to-annotation map (clusters sharing a label will be merged)

nkcd8_annotation_map <- c(
  "0"  = "CD8+ Naïve",
  "1"  = "CD8+ TEMRA",
  "2"  = "CD8+ Tem (PRF1+)",
  "3"  = "CD8+ Tem",
  "4"  = "CD8+ Tem",
  "5"  = "CD8+ TEMRA",
  "6"  = "CD56+ NK",
  "7"  = "CD56dim NK",
  "8"  = "CD8+ Memory",
  "9"  = "CD56- NK",
  "10" = "CD56dim NK",
  "11" = "CD8+ TEMRA",
  "12" = "Innate-like T",
  "13" = "CD8+ Memory",
  "14" = "CD8+ Memory",
  "15" = "ɣδ T",
  "16" = "CD8+ Memory",
  "17" = "NK-like CD8+",
  "18" = "CD8+ Tem",
  "19" = "NK-like CD8+",
  "20" = "Activated cytotoxic CD8+",
  "21" = "CD8+ TEMRA",
  "22" = "NK-like CD8+",
  "23" = "Activated cytotoxic CD8+",
  "24" = "CD8+ Tem"
)

# Apply annotations
# Assumes the active Ident or a metadata column named 'seurat_clusters' holds cluster IDs
OPIS_NKCD8$celltype_annotation <- unname(nkcd8_annotation_map[
  as.character(OPIS_NKCD8$Subcluster_ID)
])

Idents(OPIS_NKCD8) <- "celltype_annotation"

message("NK/CD8 annotation summary:")
print(table(OPIS_NKCD8$celltype_annotation))

# =============================================================================
# 3) Annotate CD4 Object
# =============================================================================
# Clusters 9 and 10 are flagged for removal

cd4_annotation_map <- c(
  "0"  = "CD4+ Naïve",
  "1"  = "CD4+ Tcm",
  "2"  = "CD4+ Tcm",
  "3"  = "CD4+ Naïve",
  "4"  = "CD4+ Treg",
  "5"  = "CD4+ Tcm",
  "6"  = "Activated CD4+ Naïve",
  "7"  = "CD4+ KLRG1+ Tem",
  "8"  = "CD4+ Tcm",
  "9"  = "Remove",
  "10" = "Remove",
  "11" = "CD4+ term memory",
  "12" = "Remove",
  "13" = "Activated CD4+ memory",
  "14" = "CD4+ Tscm-like "
)

# Apply annotations
OPIS_CD4$celltype_annotation <- unname(cd4_annotation_map[
  as.character(OPIS_CD4$Subcluster_ID)
])

# Remove clusters labelled "Remove"
cells_to_keep <- WhichCells(OPIS_CD4, expression = celltype_annotation != "Remove")

n_removed <- ncol(OPIS_CD4) - length(cells_to_keep)
message(sprintf("Removing %d cells from CD4 clusters 9 & 10.", n_removed))

OPIS_CD4 <- subset(OPIS_CD4, cells = cells_to_keep)

Idents(OPIS_CD4) <- "celltype_annotation"

message("CD4 annotation summary (after removal):")
print(table(OPIS_CD4$celltype_annotation))

# =============================================================================
# 4) Save Annotated Objects
# =============================================================================

qs_save(OPIS_NKCD8, file.path(load.path, "OPIS_NKCD8_Annotated.qs2"))
qs_save(OPIS_CD4,   file.path(load.path, "OPIS_CD4_Annotated.qs2"))

message("Annotated objects saved to: ", load.path)

# =============================================================================
# 5) Quick QC DimPlots
# =============================================================================

nkcd8_plot <- DimPlot2(
  OPIS_NKCD8,
  reduction = "wnn.umap",
  label     = TRUE,
  repel     = TRUE,
  group.by  = "celltype_annotation",
  box=T,
  cols = 'default',
  pt.size = 1
) + ggtitle("NK/CD8 – Annotated") + NoLegend()

cd4_plot <- DimPlot2(
  OPIS_CD4,
  reduction = "wnn.umap",
  label     = TRUE,
  repel     = TRUE,
  group.by  = "celltype_annotation",
  box=T,
  cols = 'default',
  pt.size = 1
) + ggtitle("CD4 – Annotated") + NoLegend()
cd4_plot
combined_plot <- nkcd8_plot | cd4_plot

ggsave(
  filename = file.path(subclust.root, "Annotated_DimPlots.png"),
  plot     = combined_plot,
  width    = 18,
  height   = 8,
  bg='white'
)

################################# Trajectory Analysis ####################
# =============================================================================
# OPIS ECCITEseq - Trajectory Analysis (Monocle3)
# CD4 & CD8 (NK clusters excluded)
# =============================================================================

library(monocle3)
library(SeuratWrappers)
library(Seurat)
library(SingleCellExperiment)
library(tidyverse)

# ---- Output path -------------------------------------------------------------
subclust.root <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/subclustering/post_annotation"
traj.out      <- file.path(subclust.root, "trajectory")
dir.create(traj.out, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1) Fix Default Assay — must be RNA before anything else
# =============================================================================

DefaultAssay(OPIS_NKCD8) <- "RNA"
DefaultAssay(OPIS_CD4)   <- "RNA"

# =============================================================================
# 2) CD8 OBJECT — Exclude NK clusters, keep pure CD8 T cells
# =============================================================================

nk_labels <- c("CD56+ NK", "CD56dim NK", "CD56- NK")

OPIS_CD8 <- subset(
  OPIS_NKCD8,
  subset = celltype_annotation %in% nk_labels,
  invert = TRUE
)

message("CD8-only dims: ", paste(dim(OPIS_CD8), collapse = " x "))
message("CD8-only cell type breakdown:")
print(table(OPIS_CD8$celltype_annotation))

# ---- Build CDS from raw RNA counts ------------------------------------------
counts_cd8 <- GetAssayData(OPIS_CD8, assay = "RNA", layer = "counts")

cds_cd8 <- new_cell_data_set(
  expression_data = counts_cd8,
  cell_metadata   = OPIS_CD8@meta.data,
  gene_metadata   = data.frame(
    gene_short_name = rownames(counts_cd8),
    row.names       = rownames(counts_cd8)
  )
)

# ---- Inject WNN UMAP ---------------------------------------------------------
reducedDims(cds_cd8)[["UMAP"]] <- Embeddings(OPIS_CD8, reduction = "wnn.umap")

# ---- Cluster and learn graph -------------------------------------------------
cds_cd8 <- cluster_cells(cds_cd8, reduction_method = "UMAP")
cds_cd8 <- learn_graph(cds_cd8, use_partition = TRUE)

# ---- Root at CD8+ Naïve ------------------------------------------------------
root_cells_cd8 <- rownames(subset(
  colData(cds_cd8),
  celltype_annotation == "CD8+ Naïve"
))
cds_cd8 <- order_cells(cds_cd8, root_cells = root_cells_cd8)

message(sprintf("CD8 trajectory rooted at %d CD8+ Naïve cells.", length(root_cells_cd8)))

# =============================================================================
# 3) CD4 OBJECT
# =============================================================================

# ---- Build CDS from raw RNA counts ------------------------------------------
counts_cd4 <- GetAssayData(OPIS_CD4, assay = "RNA", layer = "counts")

cds_cd4 <- new_cell_data_set(
  expression_data = counts_cd4,
  cell_metadata   = OPIS_CD4@meta.data,
  gene_metadata   = data.frame(
    gene_short_name = rownames(counts_cd4),
    row.names       = rownames(counts_cd4)
  )
)

# ---- Inject WNN UMAP ---------------------------------------------------------
reducedDims(cds_cd4)[["UMAP"]] <- Embeddings(OPIS_CD4, reduction = "wnn.umap")

# ---- Cluster and learn graph -------------------------------------------------
cds_cd4 <- cluster_cells(cds_cd4, reduction_method = "UMAP")
cds_cd4 <- learn_graph(cds_cd4, use_partition = TRUE)

# ---- Root at Naïve -----------------------------------------------------------
root_cells_cd4 <- rownames(subset(
  colData(cds_cd4),
  celltype_annotation == "CD4+ Naïve"
))
cds_cd4 <- order_cells(cds_cd4, root_cells = root_cells_cd4)

message(sprintf("CD4 trajectory rooted at %d Naïve cells.", length(root_cells_cd4)))

# =============================================================================
# 4) Trajectory Plots
# =============================================================================

# ---- CD8 pseudotime ----------------------------------------------------------
p_cd8_pseudo <- plot_cells(
  cds_cd8,
  color_cells_by      = "pseudotime",
  label_cell_groups   = FALSE,
  label_leaves        = TRUE,
  label_branch_points = TRUE
) + ggtitle("CD8 T cells - Pseudotime")

# ---- CD8 annotation ----------------------------------------------------------
p_cd8_annot <- plot_cells(
  cds_cd8,
  color_cells_by      = "celltype_annotation",
  label_cell_groups   = TRUE,
  label_leaves        = FALSE,
  label_branch_points = FALSE
) + ggtitle("CD8 T cells - Cell Type")

# ---- CD4 pseudotime ----------------------------------------------------------
p_cd4_pseudo <- plot_cells(
  cds_cd4,
  color_cells_by      = "pseudotime",
  label_cell_groups   = FALSE,
  label_leaves        = TRUE,
  label_branch_points = TRUE
) + ggtitle("CD4 T cells - Pseudotime")

# ---- CD4 annotation ----------------------------------------------------------
p_cd4_annot <- plot_cells(
  cds_cd4,
  color_cells_by      = "celltype_annotation",
  label_cell_groups   = TRUE,
  label_leaves        = FALSE,
  label_branch_points = FALSE
) + ggtitle("CD4 T cells - Cell Type")

# ---- Save trajectory plots ---------------------------------------------------
ggsave(file.path(traj.out, "CD8_Trajectory_Pseudotime.png"),  p_cd8_pseudo, width = 8, height = 7, dpi = 300)
ggsave(file.path(traj.out, "CD8_Trajectory_Annotation.png"),  p_cd8_annot,  width = 8, height = 7, dpi = 300)
ggsave(file.path(traj.out, "CD4_Trajectory_Pseudotime.png"),  p_cd4_pseudo, width = 8, height = 7, dpi = 300)
ggsave(file.path(traj.out, "CD4_Trajectory_Annotation.png"),  p_cd4_annot,  width = 8, height = 7, dpi = 300)

message("Trajectory plots saved.")

# =============================================================================
# 5) Transfer Pseudotime Back into Seurat Objects
# =============================================================================

OPIS_CD8$pseudotime <- pseudotime(cds_cd8)[colnames(OPIS_CD8)]
OPIS_CD4$pseudotime <- pseudotime(cds_cd4)[colnames(OPIS_CD4)]

p_cd8_feat <- FeaturePlot(OPIS_CD8, features = "pseudotime", reduction = "wnn.umap") +
  scale_colour_viridis_c() +
  ggtitle("CD8 - Pseudotime (Seurat)")

p_cd4_feat <- FeaturePlot(OPIS_CD4, features = "pseudotime", reduction = "wnn.umap") +
  scale_colour_viridis_c() +
  ggtitle("CD4 - Pseudotime (Seurat)")

ggsave(file.path(traj.out, "CD8_Pseudotime_SeuratUMAP.png"), p_cd8_feat, width = 7, height = 6, dpi = 300)
ggsave(file.path(traj.out, "CD4_Pseudotime_SeuratUMAP.png"), p_cd4_feat, width = 7, height = 6, dpi = 300)

message("Pseudotime transferred back to Seurat objects and plotted.")

# =============================================================================
# 6) Genes That Change Along Pseudotime
# =============================================================================

# ---- CD8: graph test ---------------------------------------------------------
message("Running graph_test for CD8 (this may take a few minutes)...")
deg_cd8 <- graph_test(cds_cd8, neighbor_graph = "principal_graph", cores = 4)

top_cd8 <- deg_cd8 %>%
  rownames_to_column("gene") %>%
  filter(status == "OK", q_value < 0.05, morans_I > 0.3) %>%
  arrange(desc(morans_I))


message("Top 10 pseudotime-varying genes (CD8):")
print(top_cd8[, c("gene", "morans_I", "q_value")])


# ---- CD4: graph test ---------------------------------------------------------
message("Running graph_test for CD4 (this may take a few minutes)...")
deg_cd4 <- graph_test(cds_cd4, neighbor_graph = "principal_graph", cores = 4)

top_cd4 <- deg_cd4 %>%
  rownames_to_column("gene") %>%
  filter(status == "OK", q_value < 0.05, morans_I > 0.3) %>%
  arrange(desc(morans_I))


message("Top 10 pseudotime-varying genes (CD4):")
print(top_cd4[, c("gene", "morans_I", "q_value")])

# ---- Shared top genes --------------------------------------------------------
shared_genes <- intersect(top_cd8$gene, top_cd4$gene)
message(sprintf("Shared top pseudotime genes: %s",
                ifelse(length(shared_genes) == 0, "none in top 10",
                       paste(shared_genes, collapse = ", "))))
# ---- Save DEG tables ---------------------------------------------------------
write.csv(
  deg_cd8 %>% rownames_to_column("gene") %>% arrange(q_value) %>% filter(status == "OK"),
  file.path(traj.out, "CD8_PseudotimeDEG.csv"), row.names = FALSE
)
write.csv(
  deg_cd4 %>% rownames_to_column("gene") %>% arrange(q_value) %>% filter(status == "OK"),
  file.path(traj.out, "CD4_PseudotimeDEG.csv"), row.names = FALSE
)
# =============================================================================
# 7) Gene Expression Plots Along Pseudotime
# =============================================================================

cd8_genes <- top_cd8$gene
cd4_genes <- top_cd4$gene

# ---- CD8: UMAP colored by gene expression ------------------------------------
p_cd8_genes <- plot_cells(
  cds_cd8,
  genes               = cd8_genes,
  label_cell_groups   = FALSE,
  label_leaves        = FALSE,
  label_branch_points = FALSE
) + ggtitle("CD8 - Top Pseudotime-Varying Genes")

ggsave(
  filename = file.path(traj.out, "CD8_Pseudotime_TopGenes.png"),
  plot     = p_cd8_genes,
  width    = 14, height = 12, dpi = 300
)

# ---- CD4: UMAP colored by gene expression ------------------------------------
p_cd4_genes <- plot_cells(
  cds_cd4,
  genes               = cd4_genes,
  label_cell_groups   = FALSE,
  label_leaves        = FALSE,
  label_branch_points = FALSE
) + ggtitle("CD4 - Top Pseudotime-Varying Genes")

ggsave(
  filename = file.path(traj.out, "CD4_Pseudotime_TopGenes.png"),
  plot     = p_cd4_genes,
  width    = 14, height = 12, dpi = 300
)

# ---- CD8: plot_genes_in_pseudotime (50 genes) --------------------------------
cds_cd8_subset <- cds_cd8[top_cd8$gene, ]

p_cd8_curves <- plot_genes_in_pseudotime(
  cds_subset     = cds_cd8_subset,
  min_expr       = 0.1,
  cell_size      = 0.5,
  ncol           = 5,
  color_cells_by = "celltype_annotation",
  trend_formula  = "~ splines::ns(pseudotime, df=3)"
)

ggsave(
  filename = file.path(traj.out, "CD8_Pseudotime_GeneCurves.png"),
  plot     = p_cd8_curves,
  width    = 15, height = 10, dpi = 300
)

# ---- CD4: plot_genes_in_pseudotime (25 genes) --------------------------------
cds_cd4_subset <- cds_cd4[top_cd4$gene, ]

p_cd4_curves <- plot_genes_in_pseudotime(
  cds_subset     = cds_cd4_subset,
  min_expr       = 0.1,
  cell_size      = 0.5,
  ncol           = 5,
  color_cells_by = "celltype_annotation",
  trend_formula  = "~ splines::ns(pseudotime, df=3)"
)

ggsave(
  filename = file.path(traj.out, "CD4_Pseudotime_GeneCurves.png"),
  plot     = p_cd4_curves,
  width    = 15, height = 10, dpi = 300
)
# ---- CD8: plot_genes_in_pseudotime (50 genes) --------------------------------
cds_cd8_subset <- cds_cd8[top_cd8$gene, ]

p_cd8_curves <- plot_genes_in_pseudotime(
  cds_subset     = cds_cd8_subset,
  min_expr       = 0.1,
  cell_size      = 0.5,
  ncol           = 5,
  color_cells_by = "celltype_annotation",
  trend_formula  = "~ splines::ns(pseudotime, df=3)"
)

ggsave(
  filename = file.path(traj.out, "CD8_Pseudotime_GeneCurves.png"),
  plot     = p_cd8_curves,
  width    = 15, height = 10, dpi = 300
)

# ---- CD4: plot_genes_in_pseudotime (25 genes) --------------------------------
cds_cd4_subset <- cds_cd4[top_cd4$gene, ]

p_cd4_curves <- plot_genes_in_pseudotime(
  cds_subset     = cds_cd4_subset,
  min_expr       = 0.1,
  cell_size      = 0.5,
  ncol           = 5,
  color_cells_by = "celltype_annotation",
  trend_formula  = "~ splines::ns(pseudotime, df=3)"
)

ggsave(
  filename = file.path(traj.out, "CD4_Pseudotime_GeneCurves.png"),
  plot     = p_cd4_curves,
  width    = 8, height = 4, dpi = 300
)
# ---- CD8: colored by pseudotime ----------------------------------------------
p_cd8_curves_pt <- plot_genes_in_pseudotime(
  cds_subset     = cds_cd8_subset,
  min_expr       = 0.1,
  cell_size      = 0.5,
  ncol           = 5,
  color_cells_by = "pseudotime",
  trend_formula  = "~ splines::ns(pseudotime, df=3)"
)

ggsave(
  filename = file.path(traj.out, "CD8_Pseudotime_GeneCurves_byPseudotime.png"),
  plot     = p_cd8_curves_pt,
  width    = 15, height = 10, dpi = 300
)

# ---- CD4: colored by pseudotime ----------------------------------------------
p_cd4_curves_pt <- plot_genes_in_pseudotime(
  cds_subset     = cds_cd4_subset,
  min_expr       = 0.1,
  cell_size      = 0.5,
  ncol           = 5,
  color_cells_by = "pseudotime",
  trend_formula  = "~ splines::ns(pseudotime, df=3)"
)

ggsave(
  filename = file.path(traj.out, "CD4_Pseudotime_GeneCurves_byPseudotime.png"),
  plot     = p_cd4_curves_pt,
  width    = 15, height = 8, dpi = 300
)
