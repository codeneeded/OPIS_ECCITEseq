#Load Required Libraries
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
library(cowplot)
library(ggplot2)
library(gridExtra)
library(Azimuth)
library(ggrepel)
library(patchwork)
library(reticulate)
library(circlize)
library(ComplexHeatmap)
library(readxl)
library(scCustomize)
library(SeuratWrappers)
library(qs2)
library(future)
library(future.apply)

# Set up a parallel plan with multiple workers
plan(multisession, workers = 30)  # Adjust the number of workers based on your CPU cores

##set path to load data


##set path to load data
setwd('~/Documents/OPIS_ECCITEseq/Integration/Isotype')

load.path <- "~/Documents/OPIS_ECCITEseq/saved_R_data/"


seurat_layer <- qs_read(paste0(load.path,'Seuratv5_filtered_CellCycle_DoubletClean.qs2'))

#### Delete all assays that are not RNA and ADT (including the sketched assay- this is to control size)
### we also get rid of the azimuth assays, however we keep the azimuth annotations

# Define assays to keep
assays_to_keep <- c("RNA", "ADT")

# Identify assays to drop
current_assays <- names(seurat_layer@assays)
assays_to_drop <- setdiff(current_assays, assays_to_keep)

# Drop unwanted assays (Azimuth score assays etc.), but keep meta.data annotations
for (assay in assays_to_drop) {
  seurat_layer[[assay]] <- NULL
}

# Sanity check
names(seurat_layer@assays)


########################Calculate thresholds for Isotype Controls ##################################

#### Define the isotype controls and extract data for plotting:
isotype_genes <- c('Mouse-IgG1', 'Mouse-IgG2a','Mouse-IgG2b', 'Rat-IgG2b'
                   ,'Rat-IgG1', 'Rat-IgG2a', 'Hamster-IgG')

# Extract the data for isotype controls
isotype_data <- seurat_layer[["ADT"]]$data[isotype_genes, ]
seurat_layer[["ADT"]]$data['Mouse-IgG1',]
rownames(seurat_layer[["ADT"]]$data)

# Convert the matrix to a data frame, retain cell barcodes, and reshape to long format
plot_data <- as.data.frame(t(isotype_data)) %>%
  rownames_to_column("CellBarcode") %>%
  gather(key = "Isotype", value = "Expression", -CellBarcode)

# Add the orig.ident information using the cell barcodes
plot_data$orig.ident <- seurat_layer@meta.data[plot_data$CellBarcode, "orig.ident", drop = TRUE]


#### Calculate the 99% threshold for each isotype control, split by sample:
threshold_data_list <- lapply(isotype_genes, function(isotype) {
  thresholds <- tapply(seurat_layer[["ADT"]]$data[isotype, ], seurat_layer@meta.data$orig.ident, function(x) quantile(x, 0.99))
  data.frame(Isotype = isotype, orig.ident = names(thresholds), Threshold = as.numeric(thresholds))
})

threshold_data <- do.call(rbind, threshold_data_list)

seurat_layer@assays$ADT@data@Dimnames[[1]]
#### Plot using ggplot2 with the added threshold data:

ggplot(plot_data, aes(x = Isotype, y = Expression)) +
  geom_violin(scale = "width", fill = "lightblue") +
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.5) + 
  geom_line(data = threshold_data, aes(y = Threshold, group = orig.ident), color = "red", linetype = "solid", size = 0.5) +
  facet_wrap(~ orig.ident) +  # Split the plot by sample
  theme_bw() +
  labs(title = "Isotype Control Expression with 99% Thresholds", y = "Expression Level") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave('Isotype_Threshold.png',dpi=500, width = 17)

########################################### Isotype Correction #####################################
isotype_genes <- c('Mouse-IgG1', 'Mouse-IgG2a','Mouse-IgG2b',
                   'Rat-IgG2b','Rat-IgG1', 'Rat-IgG2a', 'Hamster-IgG')

# Read the Excel file
isotype_df <- read_excel("/home/akshay-iyer/Documents/OPIS_ECCITEseq/Isotype.xlsx")

# Named vector: names = proteins, values = isotype
protein_to_isotype <- setNames(isotype_df$isotype, isotype_df$names)

# ADT proteins actually present
adt_proteins <- rownames(seurat_layer[["ADT"]]@data)

# Keep only proteins that are present in ADT
missing_proteins <- setdiff(names(protein_to_isotype), adt_proteins)
if (length(missing_proteins) > 0) {
  message(
    "Proteins in Isotype.xlsx but not in ADT (skipped in correction):\n",
    paste(missing_proteins, collapse = ", ")
  )
}
protein_to_isotype <- protein_to_isotype[names(protein_to_isotype) %in% adt_proteins]

# Calculate 99% threshold per isotype per sample
names(isotype_genes) <- isotype_genes
thresholds_list <- lapply(isotype_genes, function(isotype) {
  tapply(
    seurat_layer[["ADT"]]@data[isotype, ],
    seurat_layer@meta.data$orig.ident,
    function(x) quantile(x, 0.99)
  )
})
names(thresholds_list) <- isotype_genes

# Also ensure all mapped isotypes exist in thresholds_list
protein_to_isotype <- protein_to_isotype[protein_to_isotype %in% names(thresholds_list)]

# Copy object
seurat_isotype <- seurat_layer

# Apply thresholds
for (sample in unique(seurat_layer@meta.data$orig.ident)) {
  sample_cells <- which(seurat_layer@meta.data$orig.ident == sample)
  
  for (protein in names(protein_to_isotype)) {
    iso <- protein_to_isotype[[protein]]
    
    # safety: skip if something is off
    if (!protein %in% rownames(seurat_isotype[["ADT"]]@data)) next
    if (is.null(thresholds_list[[iso]][sample])) next
    
    thr <- thresholds_list[[iso]][sample]
    
    # subtract threshold
    seurat_isotype[["ADT"]]@data[protein, sample_cells] <-
      seurat_isotype[["ADT"]]@data[protein, sample_cells] - thr
    
    # zero negatives
    neg_idx <- which(seurat_isotype[["ADT"]]@data[protein, ] < 0)
    if (length(neg_idx) > 0) {
      seurat_isotype[["ADT"]]@data[protein, neg_idx] <- 0
    }
  }
}

qs_save(seurat_isotype, file=paste0(load.path,"Seurat_isotype.qs2"))

