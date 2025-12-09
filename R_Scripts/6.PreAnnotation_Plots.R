#Load Required Libraries
library(Seurat)
library(scater)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(hdf5r)
library(data.table)
library(cowplot)
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
library(readxl)
library(SeuratExtend)
library(qs2)  
### Read Files in

setwd('/home/akshay-iyer/Documents/OPIS_ECCITEseq/Annotation')

load.path <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/saved_R_data/"



### TARA
OPIS_ALL <- qs_read(file.path(load.path, "OPIS_ALL_WNN.qs2"))


######## Rejoin layers ####
OPIS_ALL[["RNA"]] <- JoinLayers(OPIS_ALL[["RNA"]])
OPIS_ALL[["ADT"]] <- JoinLayers(OPIS_ALL[["ADT"]])



#################################### RNA and Protein Features of interest ##############################################library(dplyr)


rna.features <- c(
  'ASCL2','BATF','BATF3','BCL6','BACH2','C1QBP','CCL2','CCL3','CCL4L2','CCL5','CCND3','CD14','CD19','CD1C',
  'CD200','CD27','CD3D','CD3E','CD36','CD4','CD40','CD40LG','CD70','CD7','CD79A','CD8A','CD8B',
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

Idents(OPIS_ALL) <- 'wsnn_res.0.4'

#Remove clusters with <50 cells
# Identify clusters with fewer than 20 cells
cluster_sizes <- table(Idents(OPIS_ALL))
small_clusters <- names(cluster_sizes[cluster_sizes < 100])
large_clusters <- names(cluster_sizes[cluster_sizes > 20])
OPIS_ALL <- subset(OPIS_ALL, idents = large_clusters)

############# Plots ###############################
setwd('/home/akshay-iyer/Documents/OPIS_ECCITEseq/Annotation/Pre-Annotation/Cluster_Plots')
p1 <- DimPlot2(
  OPIS_ALL,
  reduction = "wnn.umap",
  group.by = "wsnn_res.0.4",
  label=T,
  repel = T,
  box=T,
  label.size = 5
) + ggtitle("wsnn_res.0.4")

ggsave(
  filename = "OPIS_WNN_wsnn_res.0.4.png",
  plot=p1,
  width = 10,  # Adjust width as needed
  height = 8,  # Adjust height as needed
  dpi = 300,
  bg='white'
)
p2 <- DimPlot2(
  OPIS_ALL,
  reduction = "wnn.umap",
  group.by = "predicted.celltype.l2",
  label = TRUE,
  repel = TRUE,
  box=T,
  label.size = 5
) + ggtitle("Azimuth 2")


combined_wnn <- wrap_plots(p1, p2, ncol = 2, nrow = 1)

ggsave(
  filename = "OPIS_ALL_WNN_Clusters_Azimuth_0.4.png",
  plot     = combined_wnn,
  width    = 19,
  height   = 9,
  dpi      = 300,
  bg='white'
)


p3 <- ClusterDistrBar(OPIS_ALL$orig.ident, OPIS_ALL$wsnn_res.0.4, cols = "default", flip = FALSE, border = "black") +
  theme(axis.title.x = element_blank())


ggsave(
  filename = "OPIS_Cluster_Distribution_by_sample.png",
  plot = p3,
  width = 19,
  height = 11,
  dpi = 300,
  bg = "white"
)
OPIS_ALL$OUD_status
p4 <- ClusterDistrBar(OPIS_ALL$OUD_status, OPIS_ALL$wsnn_res.0.4, cols = "default", flip = FALSE, border = "black") +
  theme(axis.title.x = element_blank())

p4
ggsave(
  filename = "OPIS_Cluster_Distribution_by_OUD.png",
  plot = p4,
  width = 13,
  height = 9,
  dpi = 300,
  bg = "white"
)
qs_save(
  OPIS_ALL,
  file = file.path(load.path, "OPIS_ALL_PreAnnotation.qs2")
)
########################################## Feature Plots and VLN Plots ###############################################


# Make sure RNA features are present in the object
rna.features <- intersect(rna.features, rownames(OPIS_ALL[["RNA"]]))

# ----------------------------------------- #
#        Base directory + subfolders        #
# ----------------------------------------- #

base_dir <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/Annotation/Pre-Annotation"

# ADT subfolders
adt_vln_dir          <- file.path(base_dir, "Violin_Plot", "ADT")
adt_vln_OUD_dir      <- file.path(base_dir, "Violin_Plot", "ADT_splitBy_OUD_status")
adt_feature_dir      <- file.path(base_dir, "Feature_Plot", "ADT")

# RNA subfolders
rna_vln_dir          <- file.path(base_dir, "Violin_Plot", "RNA")
rna_vln_OUD_dir      <- file.path(base_dir, "Violin_Plot", "RNA_splitBy_OUD_status")
rna_feature_dir      <- file.path(base_dir, "Feature_Plot", "RNA")

dirs <- c(
  adt_vln_dir, adt_vln_OUD_dir, adt_feature_dir,
  rna_vln_dir, rna_vln_OUD_dir, rna_feature_dir
)
lapply(dirs, function(x) dir.create(x, recursive = TRUE, showWarnings = FALSE))

# ========================================= #
#                    ADT                    #
# ========================================= #

DefaultAssay(OPIS_ALL) <- "ADT"

for (i in prots) {
  if (i %in% rownames(OPIS_ALL[["ADT"]])) {
    
    # ------------------------------- #
    #   1) SPLIT-BY OUD STATUS FIRST  #
    # ------------------------------- #
    vln.pl.split <- VlnPlot2(
      OPIS_ALL,
      features   = i,
      cols       = "default",
      split.by   = "OUD_status",
      stat.method = "wilcox.test",
      show.mean  = TRUE,
      mean_colors = c("red", "blue")
    ) + ggtitle(paste("ADT |", i, "| Split by OUD_status"))
    
    ggsave(
      filename = file.path(adt_vln_OUD_dir, paste0(i, "_ADT_VLNplot_splitBy_OUD_status.png")),
      plot     = vln.pl.split,
      dpi      = 500,
      width    = 14,
      height   = 8,
      bg       = "white"
    )
    
    # ------------------------------- #
    #   2) NON-SPLIT VIOLIN PLOT      #
    # ------------------------------- #
    vln.pl <- VlnPlot2(
      OPIS_ALL,
      features   = i,
      cols       = "default",
      show.mean  = TRUE,
      mean_colors = c("red", "blue")
    ) + ggtitle(paste("ADT |", i))
    
    ggsave(
      filename = file.path(adt_vln_dir, paste0(i, "_ADT_VLNplot.png")),
      plot     = vln.pl,
      dpi      = 500,
      width    = 14,
      height   = 8,
      bg       = "white"
    )
    
    # ------------------------------- #
    #   3) FEATURE PLOT               #
    # ------------------------------- #
    pal <- viridis(n = 10, option = "A")
    fea.pl <- FeaturePlot_scCustom(
      OPIS_ALL,
      reduction  = "wnn.umap",
      features   = i,
      colors_use = pal,
      order      = TRUE
    )
    
    ggsave(
      filename = file.path(adt_feature_dir, paste0(i, "_ADT_Featureplot_Magma.png")),
      plot     = fea.pl,
      dpi      = 500,
      width    = 8
    )
  }
}

# ========================================= #
#                    RNA                    #
# ========================================= #

DefaultAssay(OPIS_ALL) <- "RNA"

for (i in rna.features) {
  if (i %in% rownames(OPIS_ALL[["RNA"]])) {
    
    # ------------------------------- #
    #   1) SPLIT-BY OUD STATUS FIRST  #
    # ------------------------------- #
    vln.pl.split <- VlnPlot2(
      OPIS_ALL,
      features   = i,
      cols       = "default",
      split.by   = "OUD_status",
      stat.method = "wilcox.test",
      show.mean  = TRUE,
      mean_colors = c("red", "blue")
    ) + ggtitle(paste("RNA |", i, "| Split by OUD_status"))
    
    ggsave(
      filename = file.path(rna_vln_OUD_dir, paste0(i, "_RNA_VLNplot_splitBy_OUD_status.png")),
      plot     = vln.pl.split,
      dpi      = 500,
      width    = 14,
      height   = 8,
      bg       = "white"
    )
    
    # ------------------------------- #
    #   2) NON-SPLIT VIOLIN PLOT      #
    # ------------------------------- #
    vln.pl <- VlnPlot2(
      OPIS_ALL,
      features   = i,
      cols       = "default",
      show.mean  = TRUE,
      mean_colors = c("red", "blue")
    ) + ggtitle(paste("RNA |", i))
    
    ggsave(
      filename = file.path(rna_vln_dir, paste0(i, "_RNA_VLNplot.png")),
      plot     = vln.pl,
      dpi      = 500,
      width    = 14,
      height   = 8,
      bg       = "white"
    )
    
    # ------------------------------- #
    #   3) FEATURE PLOT               #
    # ------------------------------- #
    pal <- viridis(n = 10, option = "A")
    
    fea.pl <- FeaturePlot_scCustom(
      OPIS_ALL,
      reduction  = "wnn.umap",
      features   = i,
      colors_use = pal,
      order      = TRUE
    )
    
    ggsave(
      filename = file.path(rna_feature_dir, paste0(i, "_RNA_Featureplot_Magma.png")),
      plot     = fea.pl,
      dpi      = 500,
      width    = 8
    )
  }
}


