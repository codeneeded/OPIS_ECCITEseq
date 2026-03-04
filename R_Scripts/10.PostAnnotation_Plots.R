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



### 
OPIS_ALL <- qs_read(file.path(load.path, "OPIS_ALL_PostAnnotation.qs2"))

levels(as.factor(OPIS_ALL$Cluster_ID))
###################### Features ######################
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
########################################## Feature Plots and VLN Plots ###############################################


# Make sure RNA features are present in the object
rna.features <- intersect(rna.features, rownames(OPIS_ALL[["RNA"]]))

# ----------------------------------------- #
#        Base directory + subfolders        #
# ----------------------------------------- #

base_dir <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/Annotation/Post-Annotation"

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

