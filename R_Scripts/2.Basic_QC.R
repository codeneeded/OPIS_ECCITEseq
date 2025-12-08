# Quality Control Visualizations
##Ref -> https://www.singlecellcourse.org/scrna-seq-analysis-with-bioconductor.html
##Ref -> http://bioconductor.org/books/3.15/OSCA.intro/getting-scrna-seq-datasets.html
#Ref -> https://hbctraining.github.io/scRNA-seq_online/schedule/links-to-lessons.html

library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(hdf5r)
library(data.table)
library(ggplot2)
library(biomaRt)
library(scCustomize)
library(patchwork)
library(qs2)  

############################################ Global Variables #######################################################

# Project root
project_root <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq"

# QC root folder
qc_root <- file.path(project_root, "QC")
if (!dir.exists(qc_root)) dir.create(qc_root, recursive = TRUE)

# Set working directory to QC root
setwd(qc_root)

# Define input and output paths
in.path  <- file.path(project_root, "10x_Genomics")     # if needed later
out.path <- file.path(project_root, "saved_R_data")     # where .qs2 objects live

# Load merged Seurat object (saved previously with qs)
merged_seurat <- qs_read(file.path(out.path, "Seuratv5_CITEseq_dsbnorm_merged_Seurat.qs2"))

############################################## Cell Level QC #########################################################


## Number of Cells per Sample

### Create Pre-QC Folder as needed
preqc_dir <- file.path(qc_root, "Pre-QC")
if (!dir.exists(preqc_dir)) {
  dir.create(preqc_dir)
  message("Folder 'Pre-QC' created.")
} else {
  message("Folder 'Pre-QC' already exists.")
}
setwd(preqc_dir)

metadata <- merged_seurat@meta.data

png(file = "Cells_per_sample.png", width = 1800, height = 1200)
metadata %>%
  ggplot(aes(x = orig.ident, fill = orig.ident)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x  = element_text(angle = 45, vjust = 1, hjust = 1),
        plot.title  = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("NCells")
dev.off()

png(file = "Cells_per_OUD_status.png", width = 1800, height = 1200)
metadata %>%
  ggplot(aes(x = OUD_status, fill = OUD_status)) +
  geom_bar() +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    plot.title  = element_text(hjust = 0.5, face = "bold")
  ) +
  ggtitle("NCells per OUD Status")
dev.off()

# QC Features rough look
feats.1 <- c("nCount_RNA", "nFeature_RNA", "nCount_ADT", "nFeature_ADT",
             "percent_mito", "percent_ribo", "percent_hb", "percent_plat")

png(file = "Pre-QC_features_grouped.png", width = 1800, height = 1200)
VlnPlot(merged_seurat, group.by = "orig.ident", features = feats.1,
        pt.size = 0.1, ncol = 2) +
  NoLegend()
dev.off()

# UMI Count
png(file = "UMI_Count.png", width = 1800, height = 1200)
metadata %>%
  ggplot(aes(color = orig.ident, x = nCount_RNA, fill = orig.ident)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)
dev.off()

# nGenes
png(file = "nGenes.png", width = 1800, height = 1200)
metadata %>%
  ggplot(aes(color = orig.ident, x = nFeature_RNA, fill = orig.ident)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 250)
dev.off()

# Complexity Score
png(file = "Complexity_Score.png", width = 1800, height = 1200)
metadata %>%
  ggplot(aes(x = log10GenesPerUMI, color = orig.ident, fill = orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
dev.off()

# Mito Ratio
png(file = "Mito_Ratio.png", width = 1800, height = 1200)
metadata %>%
  ggplot(aes(color = orig.ident, x = percent_mito, fill = orig.ident)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 15)
dev.off()

# Ribo Ratio
png(file = "Ribo_Ratio.png", width = 1800, height = 1200)
metadata %>%
  ggplot(aes(color = orig.ident, x = percent_ribo, fill = orig.ident)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 5)
dev.off()

# Heme Ratio
png(file = "Heme_Ratio.png", width = 1800, height = 1200)
metadata %>%
  ggplot(aes(color = orig.ident, x = percent_hb, fill = orig.ident)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 20)
dev.off()

# Platelet Ratio
png(file = "Platlet_Ratio.png", width = 1800, height = 1200)
metadata %>%
  ggplot(aes(color = orig.ident, x = percent_plat, fill = orig.ident)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 2)
dev.off()

######## SC_Customize Output
p1 <- QC_Plots_Genes(seurat_object = merged_seurat, low_cutoff = 600, high_cutoff = 5500)
p2 <- QC_Plots_UMIs(seurat_object = merged_seurat, low_cutoff = 1200, high_cutoff = 45000)
p3 <- QC_Plots_Mito(seurat_object = merged_seurat, high_cutoff = 20)
p4 <- QC_Plots_Complexity(seurat_object = merged_seurat, high_cutoff = 0.8)

png(file = "Grouped_Cuttoff.png", width = 1800, height = 1200)
wrap_plots(p1, p2, p3, p4, ncol = 4)
dev.off()

### Scatter QC
png(file = "UMIvsGene.png", width = 1800, height = 1200)
QC_Plot_UMIvsGene(
  seurat_object   = merged_seurat,
  low_cutoff_gene = 600,
  high_cutoff_gene = 5500,
  low_cutoff_UMI   = 500,
  high_cutoff_UMI  = 50000,
  group.by = "orig.ident"
)
dev.off()

png(file = "MitovsGene.png", width = 1800, height = 1200)
QC_Plot_GenevsFeature(
  seurat_object    = merged_seurat,
  feature1         = "percent_mito",
  low_cutoff_gene  = 600,
  high_cutoff_gene = 5500,
  high_cutoff_feature = 20,
  group.by = "orig.ident"
)
dev.off()

png(file = "MitovsGene_gradient.png", width = 1800, height = 1200)
QC_Plot_UMIvsGene(
  seurat_object      = merged_seurat,
  meta_gradient_name = "percent_mito",
  low_cutoff_gene    = 600,
  high_cutoff_gene   = 5500,
  high_cutoff_UMI    = 45000
)
dev.off()

#################################################### FILTERING ######################################################

merged_seurat <- JoinLayers(merged_seurat)
DefaultAssay(merged_seurat) <- "RNA"

# Filter out low quality cells using selected thresholds
filtered_seurat <- subset(
  x = merged_seurat,
  subset = (nCount_RNA      >= 500)  &
    (nFeature_RNA    >= 600)  &
    (log10GenesPerUMI > 0.80) &
    (percent_mito    < 15)    &
    (percent_ribo    > 5)     &
    (percent_hb      < 20)    &
    (percent_plat    < 2)
)

# Keep genes expressed in at least 10 cells
genes_in_10_cells <- rowSums(filtered_seurat@assays$RNA@layers$counts > 0) >= 10
filtered_seurat <- subset(filtered_seurat, features = names(genes_in_10_cells[genes_in_10_cells]))

### Make Post-QC Folder
postqc_dir <- file.path(qc_root, "Post-QC")
if (!dir.exists(postqc_dir)) {
  dir.create(postqc_dir)
  message("Folder 'Post-QC' created.")
} else {
  message("Folder 'Post-QC' already exists.")
}
setwd(postqc_dir)

# Cells Post QC
png(file = "Post-QC_Cells_per_sample.png", width = 1800, height = 1200)
filtered_seurat@meta.data %>%
  ggplot(aes(x = orig.ident, fill = orig.ident)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x  = element_text(angle = 45, vjust = 1, hjust = 1),
        plot.title  = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("NCells")
dev.off()

# Cells Post QC per OUD status (NEW)
png(file = "Post-QC_Cells_per_OUD_status.png", width = 1800, height = 1200)
filtered_seurat@meta.data %>%
  ggplot(aes(x = OUD_status, fill = OUD_status)) +
  geom_bar() +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    plot.title  = element_text(hjust = 0.5, face = "bold")
  ) +
  ggtitle("Post-QC NCells per OUD Status")
dev.off()

png(file = "Post-QC_features_grouped.png", width = 1800, height = 1200)
VlnPlot(filtered_seurat, group.by = "orig.ident", features = feats.1,
        pt.size = 0.1, ncol = 2) +
  NoLegend()
dev.off()


# Save filtered object with qs2 into saved_R_data folder
qs_save(filtered_seurat,
      file = file.path(out.path, "Seuratv5_filtered_seurat.qs2"))
