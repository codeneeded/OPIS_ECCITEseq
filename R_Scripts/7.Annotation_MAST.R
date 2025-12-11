# ===============================
# OPIS_ALL: Annotation + DGE/DPE
# ===============================

library(Seurat)
library(dplyr)
library(ggplot2)
library(qs2)
library(SeuratExtend)
# -------------------------------------------------
# 0) Load OPIS_ALL from qs2
# -------------------------------------------------

load.path <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/saved_R_data"

OPIS_ALL <- qs_read(
  file.path(load.path, "OPIS_ALL_PreAnnotation.qs2")
)

DefaultAssay(OPIS_ALL) <- "RNA"

# =================
# ============================
# 1) Apply annotation from CSV
# ============================

annot_csv <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/Annotation/OPIS final annotation.csv"
annot_df <- read.csv(annot_csv, check.names = FALSE, stringsAsFactors = FALSE)

annot_df <- annot_df %>%
  dplyr::rename(
    Cluster_ID           = `Cluster ID`,
    Confirmed_Annotation = `Confirmed annotation`
  )

# Make sure cluster ID matches your clustering column
OPIS_ALL$Cluster_ID <- as.character(OPIS_ALL$wsnn_res.0.4)

# Map numeric cluster ID (0–25) to annotation
annot_map <- setNames(
  annot_df$Confirmed_Annotation,
  as.character(annot_df$Cluster_ID)
)

# Look up annotation per cell based on its Cluster_ID, and drop names
OPIS_ALL$Manual_Annotation <- unname(annot_map[OPIS_ALL$Cluster_ID])


# ============================
# 2) Remove clusters 23, 24, 25
# ============================

OPIS_ALL <- subset(OPIS_ALL, subset = !(wsnn_res.0.4 %in% c(23, 24, 25)))

# Drop unused factor levels just to keep things tidy
OPIS_ALL$wsnn_res.0.4      <- factor(OPIS_ALL$wsnn_res.0.4)
OPIS_ALL$Manual_Annotation <- factor(OPIS_ALL$Manual_Annotation)

OPIS_ALL$wsnn_res.0.4      <- droplevels(OPIS_ALL$wsnn_res.0.4)
OPIS_ALL$Manual_Annotation <- droplevels(OPIS_ALL$Manual_Annotation)

# ============================
# 3) Post-Annotation DimPlots
# ============================

post_annot_dir <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/Annotation/Post-Annotation"
dir.create(post_annot_dir, recursive = TRUE, showWarnings = FALSE)

# 3A) DimPlot by wsnn_res.0.4 (post-filter)
Idents(OPIS_ALL) <- "wsnn_res.0.4"

p_clusters <- DimPlot2(
  OPIS_ALL,
  reduction  = "wnn.umap",
  group.by   = "wsnn_res.0.4",
  label      = TRUE,
  repel      = TRUE,
  box        = TRUE,
  label.size = 5
) + ggtitle("OPIS_ALL | wsnn_res.0.4 (Post-Annotation)")

ggsave(
  filename = file.path(post_annot_dir, "OPIS_wsnn_res0.4_postAnnotation.png"),
  plot     = p_clusters,
  width    = 8,
  height   = 7,
  dpi      = 400,
  bg='white'
)

# 3B) DimPlot by Manual_Annotation
Idents(OPIS_ALL) <- "Manual_Annotation"

p_manual <- DimPlot2(
  OPIS_ALL,
  reduction  = "wnn.umap",
  group.by   = "Manual_Annotation",
  label      = TRUE,
  repel      = TRUE,
  box        = TRUE,
  label.size = 3
) + ggtitle("OPIS_ALL | Manual Annotation")

ggsave(
  filename = file.path(post_annot_dir, "OPIS_ManualAnnotation_postAnnotation.png"),
  plot     = p_manual,
  width    = 15,
  height   = 7,
  dpi      = 400,
  bg='white'
)

# ======================================
# 4) Helper: Cluster-wise DE function
# ======================================

safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)
count_sig <- function(df, pcol = "p_val_adj", thr = 0.05) sum(df[[pcol]] < thr, na.rm = TRUE)

run_clusterwise_de <- function(obj,
                               assay,
                               group_col,
                               ident1,
                               ident2,
                               out_dir,
                               latent  = c("nCount_RNA"),
                               min_cells_per_grp = 10,
                               title_prefix = "",
                               file_stub = NULL) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  DefaultAssay(obj) <- assay
  Idents(obj) <- "Manual_Annotation"
  
  de_list   <- list()
  count_tbl <- data.frame()
  
  cl_levels <- levels(obj)
  
  for (cl in cl_levels) {
    message("[", assay, "] Processing cluster: ", cl)
    
    cl_cells <- WhichCells(obj, idents = cl)
    sc <- subset(obj, cells = cl_cells)
    sc$.__grp <- sc[[group_col]]
    
    tab <- table(sc$.__grp)
    
    # Require both groups and min cells
    if (all(c(ident1, ident2) %in% names(tab)) &&
        all(tab[c(ident1, ident2)] >= min_cells_per_grp)) {
      
      de <- FindMarkers(
        sc,
        ident.1   = ident1,
        ident.2   = ident2,
        group.by  = ".__grp",
        test.use  = "MAST",
        latent.vars = latent
      )
      
      comp_name <- if (is.null(file_stub)) paste0(ident1, "_vs_", ident2) else file_stub
      out_csv   <- file.path(out_dir, paste0(safe_name(cl), "_", comp_name, ".csv"))
      
      write.csv(de, out_csv)
      de_list[[cl]] <- de
      
      count_tbl <- rbind(
        count_tbl,
        data.frame(Cluster = cl, DE_Genes = count_sig(de))
      )
    }
  }
  
  # Simple barplot of DE gene counts per cluster
  if (nrow(count_tbl) > 0) {
    p <- ggplot(count_tbl, aes(x = reorder(Cluster, -DE_Genes), y = DE_Genes, fill = Cluster)) +
      geom_bar(stat = "identity") +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none") +
      xlab("Cluster") +
      ylab("# of DE Features (adj. p < 0.05)") +
      ggtitle(paste0(title_prefix, ": Clusters Ranked by DE Features (",
                     ident1, " vs ", ident2, " | ", assay, ")"))
    
    ggsave(
      filename = file.path(
        out_dir,
        paste0(
          "Cluster_DE_Counts_",
          if (is.null(file_stub)) paste0(ident1, "vs", ident2) else file_stub,
          "_", assay, ".png"
        )
      ),
      plot   = p,
      width  = 10,
      height = 6,
      dpi    = 400
    )
  }
  
  invisible(list(de = de_list, counts = count_tbl))
}

# ======================================
# 5) DGE (RNA) & DPE (ADT) for OUD
# ======================================

de_base <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/Differential_Expression"
dge_dir <- file.path(de_base, "DGE", "OUD_Pos_vs_Neg")
dpe_dir <- file.path(de_base, "DPE", "OUD_Pos_vs_Neg")

dir.create(dge_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dpe_dir, recursive = TRUE, showWarnings = FALSE)

# Make sure OUD_status is a factor and remove NAs
OPIS_ALL$OUD_status <- factor(OPIS_ALL$OUD_status)
OPIS_OUD <- subset(OPIS_ALL, subset = !is.na(OUD_status))

# ---- Set the group labels here if needed ----
# Check levels(OPIS_OUD$OUD_status) to confirm
print(levels(OPIS_OUD$OUD_status))

# Adjust these if your actual level names differ:
ident1_OUD <- "OUD+"  # group 1
ident2_OUD <- "OUD-"  # group 2

# -------------------------
# 5A) DGE: RNA (MAST)
# -------------------------

dge_results <- run_clusterwise_de(
  obj         = OPIS_OUD,
  assay       = "RNA",
  group_col   = "OUD_status",
  ident1      = ident1_OUD,
  ident2      = ident2_OUD,
  out_dir     = dge_dir,
  latent      = c("nCount_RNA"),
  title_prefix = "OPIS_ALL OUD (RNA)",
  file_stub   = "OUD_Pos_vs_Neg_RNA"
)

# -------------------------
# 5B) DPE: ADT (MAST)
# -------------------------

dpe_results <- run_clusterwise_de(
  obj         = OPIS_OUD,
  assay       = "ADT",
  group_col   = "OUD_status",
  ident1      = ident1_OUD,
  ident2      = ident2_OUD,
  out_dir     = dpe_dir,
  latent      = c("nCount_RNA", "nCount_ADT"),
  title_prefix = "OPIS_ALL OUD (ADT)",
  file_stub   = "OUD_Pos_vs_Neg_ADT"
)

qs_save(
  OPIS_ALL,
  file = file.path(load.path, "OPIS_ALL_PostAnnotation.qs2")
)
