################################################################################
# OPIS ECCITE-seq — Manuscript Figure Generation & DGE Analysis
#
# Sections:
#   1. UMAP PBMC: Combine clusters → broad cell types, plot
#   2. Reannotate CD4 with updated labels
#   3. Annotation heatmaps (CD4 & CD8/NK)
#   4. Cluster-wise DGE (MAST): OUD+ vs OUD−
#
# Output directory:
#   /home/akshay-iyer/Documents/OPIS_ECCITEseq/Manuscript_Figures
################################################################################

library(Seurat)
library(SeuratObject)
library(qs2)
library(ggplot2)
library(pheatmap)
library(grid)
library(dplyr)
library(scCustomize)   # for DimPlot2
library(SeuratExtend)
# ── Global ggplot theme: white background ─────────────────────────────────────
theme_set(theme_classic(base_size = 14) + theme(
  panel.background = element_rect(fill = "white", colour = NA),
  plot.background  = element_rect(fill = "white", colour = NA)
))

# ── Paths ─────────────────────────────────────────────────────────────────────
load.path <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/saved_R_data/"
out_base  <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/Manuscript_Figures"

out_umap <- file.path(out_base, "UMAP")
out_heat <- file.path(out_base, "Heatmaps")
out_dge  <- file.path(out_base, "DGE")

dir.create(out_umap, recursive = TRUE, showWarnings = FALSE)
dir.create(out_heat, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dge,  recursive = TRUE, showWarnings = FALSE)

################################################################################
# 1) UMAP PBMC — Combine clusters to broad cell types
################################################################################

OPIS_ALL <- qs_read(file.path(load.path, "OPIS_ALL_PreAnnotation.qs2"))

# Remove clusters 26–30 (not assigned to any broad cell type)
OPIS_ALL <- subset(OPIS_ALL, wsnn_res.0.4 %in% as.character(0:25))

# Define cluster → cell type mapping (wsnn_res.0.4, clusters 0–25)
cluster_to_celltype <- c(
  # CD4+
  "5"  = "CD4+", "6"  = "CD4+", "7"  = "CD4+", "8"  = "CD4+",
  # CD8+
  "0"  = "CD8+", "2"  = "CD8+", "9"  = "CD8+", "10" = "CD8+",
  "11" = "CD8+", "13" = "CD8+", "14" = "CD8+", "18" = "CD8+",
  "19" = "CD8+",
  # NK
  "4"  = "NK", "16" = "NK",
  # B
  "3"  = "B", "21" = "B", "24" = "B", "12" = "B",
  # Monocytes
  "17" = "Monocytes", "1"  = "Monocytes", "22" = "Monocytes",
  "23" = "Monocytes", "15" = "Monocytes", "25" = "Monocytes",
  # pDC
  "20" = "pDC"
)

# Apply annotation — unname() is critical to avoid Seurat "No cell overlap" error
OPIS_ALL$broad_celltype <- unname(cluster_to_celltype[as.character(OPIS_ALL$wsnn_res.0.4)])
OPIS_ALL$broad_celltype <- factor(
  OPIS_ALL$broad_celltype,
  levels = c("CD4+", "CD8+", "NK", "B", "Monocytes", "pDC")
)

# Plot: original clusters
p1 <- DimPlot2(
  OPIS_ALL,
  reduction = "wnn.umap",
  group.by  = "wsnn_res.0.4",
  label      = TRUE,
  repel      = TRUE,
  box        = TRUE,
  label.size = 5
) + ggtitle("wsnn_res.0.4")

# Plot: broad cell types
p2 <- DimPlot2(
  OPIS_ALL,
  reduction   = "wnn.umap",
  group.by    = "broad_celltype",
  label       = TRUE,
  repel       = TRUE,
  box         = TRUE,
  label.size  = 8,
  label.color = "black",
  theme       = list(NoLegend(), NoAxes(), theme_umap_arrows())
) + ggtitle("Broad Cell Type Annotation")

# Save combined
png(file.path(out_umap, "PBMC_UMAP_Clusters_and_CellTypes.png"),
    width = 20, height = 8, units = "in", res = 300, bg = "white")
print(p1 | p2)
dev.off()

# Save individual
ggsave(file.path(out_umap, "PBMC_UMAP_Clusters.png"),
       plot = p1, width = 10, height = 8, dpi = 300)
ggsave(file.path(out_umap, "PBMC_UMAP_BroadCellTypes.png"),
       plot = p2, width = 10, height = 8, dpi = 300)

message("✓ Section 1: PBMC UMAP saved")

################################################################################
# 2) Load CD4 & NK/CD8 annotated objects + DimPlots
################################################################################

# Load objects
OPIS_CD4   <- qs_read(file.path(load.path, "OPIS_CD4_Annotated.qs2"))
OPIS_NKCD8 <- qs_read(file.path(load.path, "OPIS_NKCD8_Annotated.qs2"))

# Verify annotation columns
stopifnot("celltype_annotation" %in% colnames(OPIS_CD4@meta.data))
stopifnot("celltype_annotation" %in% colnames(OPIS_NKCD8@meta.data))

# CD4 UMAP
p_cd4_umap <- DimPlot2(
  OPIS_CD4,
  reduction   = "wnn.umap",
  group.by    = "celltype_annotation",
  label       = TRUE,
  repel       = TRUE,
  box         = TRUE,
  label.size  = 8,
  theme       = list(NoLegend(), NoAxes(), theme_umap_arrows())
) + ggtitle("CD4+ T Cell Subtypes")

ggsave(file.path(out_umap, "CD4_UMAP_Annotated.png"),
       plot = p_cd4_umap, width = 10, height = 8, dpi = 300)

# NK/CD8 UMAP
p_nkcd8_umap <- DimPlot2(
  OPIS_NKCD8,
  reduction   = "wnn.umap",
  group.by    = "celltype_annotation",
  label       = TRUE,
  repel       = TRUE,
  box         = TRUE,
  label.size  = 8,
  theme       = list(NoLegend(), NoAxes(), theme_umap_arrows())
) + ggtitle("NK & CD8+ T Cell Subtypes")

ggsave(file.path(out_umap, "NKCD8_UMAP_Annotated.png"),
       plot = p_nkcd8_umap, width = 10, height = 8, dpi = 300)

message("✓ Section 2: CD4 & NK/CD8 UMAPs saved")

################################################################################
# 3) Annotation Heatmaps — CD4, CD8, NK
################################################################################

# ── Shared utilities ──────────────────────────────────────────────────────────
scale_01 <- function(mat) {
  scaled <- t(apply(mat, 1, function(x) (x - min(x)) / (max(x) - min(x) + 1e-9)))
  colnames(scaled) <- colnames(mat)
  scaled
}

get_gaps <- function(annot_df) {
  grps <- as.character(annot_df$Group)
  which(grps[-length(grps)] != grps[-1])
}

cluster_within_groups <- function(mat, annot_df) {
  groups    <- levels(droplevels(annot_df$Group))
  row_order <- character(0)
  for (grp in groups) {
    rows <- rownames(annot_df)[annot_df$Group == grp]
    if (length(rows) <= 1) {
      row_order <- c(row_order, rows)
    } else {
      hc        <- hclust(dist(mat[rows, , drop = FALSE]), method = "ward.D2")
      row_order <- c(row_order, rows[hc$order])
    }
  }
  row_order
}

heatmap_colors <- colorRampPalette(c(
  "#F7FCF5", "#C7E9C0", "#74C476", "#31A354", "#006D2C"
))(100)

# ── Generic heatmap builder ───────────────────────────────────────────────────
# obj:           Seurat object
# assay_name:    "RNA" or "ADT"
# marker_groups: named vector (feature = group label)
# group_col:     metadata column for grouping columns (celltype_annotation)
# title:         plot title
# slot_name:     data slot to use
build_annotation_heatmap <- function(obj, assay_name, marker_groups,
                                     group_col = "celltype_annotation",
                                     title = "", slot_name = "data") {
  
  DefaultAssay(obj) <- assay_name
  
  # Filter to features actually present in the object
  available <- rownames(obj)
  keep      <- names(marker_groups) %in% available
  if (sum(keep) == 0) {
    message("  ⚠ No features found for ", assay_name, " in: ", title)
    return(NULL)
  }
  marker_groups <- marker_groups[keep]
  
  # Average expression per cluster
  avg <- AverageExpression(
    obj,
    assays   = assay_name,
    features = names(marker_groups),
    group.by = group_col,
    slot     = slot_name
  )[[assay_name]]
  
  # Scale
  if (assay_name == "ADT") {
    avg_scaled <- scale_01(log10(avg + 1))
  } else {
    avg_scaled <- scale_01(avg)
  }
  colnames(avg_scaled) <- gsub("^g ", "", colnames(avg_scaled))
  
  # Group colors — build dynamically
  all_groups <- unique(marker_groups)
  # Base palette
  base_palette <- c(
    "Naïve/Memory"                = "#74C2E1",
    "Naïve / Memory"              = "#74C2E1",
    "Naïve/Stemness"              = "#74C2E1",
    "Activation"                  = "#FFD54F",
    "Regulatory"                  = "#AB47BC",
    "Effector/Term Diff"          = "#EF5350",
    "Effector/Cytotoxic"          = "#EF5350",
    "Term Diff"                   = "#FF7043",
    "Lineage"                     = "#78909C",
    "NK Identity"                 = "#F06292",
    "NK-like Receptor Program"    = "#F06292",
    "CD56+ Less Mature"           = "#CE93D8",
    "Gamma Delta T Identity"      = "#9C6FD6",
    "Innate-like T Identity"      = "#FFAB91",
    "Memory/Stemness"             = "#66BB6A",
    "Co-stimulation"              = "#AED581",
    "Effector/TEMRA"              = "#EF5350",
    "Effector TFs"                = "#FFA726",
    "NK-like/Innate"              = "#F06292",
    "Exhaustion"                  = "#A1887F"
  )
  # Assign colors; fallback to grey for unmapped groups
  group_colors <- setNames(
    sapply(all_groups, function(g) {
      if (g %in% names(base_palette)) base_palette[g] else "#BDBDBD"
    }),
    all_groups
  )
  
  # Row annotation
  annot_df <- data.frame(
    Group = factor(marker_groups[rownames(avg_scaled)],
                   levels = unique(marker_groups)),
    row.names = rownames(avg_scaled)
  )
  
  # Cluster rows within groups
  row_order     <- cluster_within_groups(avg_scaled, annot_df)
  avg_ordered   <- avg_scaled[row_order, ]
  annot_ordered <- annot_df[row_order, , drop = FALSE]
  gaps          <- get_gaps(annot_ordered)
  
  annot_colors <- list(
    Group = group_colors[levels(droplevels(annot_ordered$Group))]
  )
  
  p <- pheatmap(
    avg_ordered,
    cluster_rows         = FALSE,
    cluster_cols         = FALSE,
    scale                = "none",
    color                = heatmap_colors,
    border_color         = "white",
    annotation_row       = annot_ordered,
    annotation_colors    = annot_colors,
    annotation_names_row = FALSE,
    gaps_row             = gaps,
    cellwidth            = 55,
    cellheight           = 16,
    fontsize             = 11,
    fontsize_row         = 10,
    fontsize_col         = 10,
    angle_col            = 45,
    main                 = title,
    silent               = TRUE
  )
  
  return(p)
}


# ══════════════════════════════════════════════════════════════════════════════
# 3A) CD4 Heatmaps
# ══════════════════════════════════════════════════════════════════════════════

cd4_rna_markers <- c(
  # Naïve / Memory
  "TCF7"    = "Naïve/Memory",
  "LEF1"    = "Naïve/Memory",
  "BACH2"   = "Naïve/Memory",
  "CCR7"    = "Naïve/Memory",
  "IL7R"    = "Naïve/Memory",
  "S100A4"  = "Naïve/Memory",
  # Activation
  "HLA-DRA" = "Activation",
  "IL2RA"   = "Activation",
  "FOS"     = "Activation",
  "JUN"     = "Activation",
  # Regulatory
  "FOXP3"   = "Regulatory",
  "TIGIT"   = "Regulatory",
  "CTLA4"   = "Regulatory",
  # Effector / Terminal Differentiation
  "KLRG1"   = "Effector/Term Diff",
  "PRF1"    = "Effector/Term Diff",
  "GZMK"    = "Effector/Term Diff",
  # Lineage
  "CD3D"    = "Lineage",
  "TRAC"    = "Lineage"
)

cd4_adt_markers <- c(
  # Naïve / Memory
  "CD45RA"  = "Naïve/Memory",
  "CD27"    = "Naïve/Memory",
  "SELL"    = "Naïve/Memory",
  "CD45RO"  = "Naïve/Memory",
  # Activation
  "CD38"    = "Activation",
  "ICOS"    = "Activation",
  "HLA_DR"  = "Activation",
  # Regulatory
  "CD25"    = "Regulatory",
  "CCR4"    = "Regulatory",
  # Effector / Terminal Differentiation
  "KLRG1"   = "Effector/Term Diff",
  # Lineage
  "CD4"     = "Lineage",
  "CD8"     = "Lineage"
)

p_cd4_rna <- build_annotation_heatmap(
  OPIS_CD4, "RNA", cd4_rna_markers,
  title = "CD4+ T Cells — mRNA (RNA)"
)

p_cd4_adt <- build_annotation_heatmap(
  OPIS_CD4, "ADT", cd4_adt_markers,
  title = "CD4+ T Cells — Protein (ADT)"
)

# Save CD4 combined
if (!is.null(p_cd4_adt) && !is.null(p_cd4_rna)) {
  png(file.path(out_heat, "CD4_ADT_RNA_Heatmap_Combined.png"),
      width = 24, height = 14, units = "in", res = 300, bg = "white")
  grid.newpage()
  grid.text("CD4+ T Cell Annotation Validation: Protein (ADT) & mRNA",
            x = 0.5, y = 0.975,
            gp = gpar(fontsize = 18, fontface = "bold"))
  pushViewport(viewport(layout = grid.layout(1, 2)))
  pushViewport(viewport(layout.pos.col = 1)); grid.draw(p_cd4_adt$gtable); popViewport()
  pushViewport(viewport(layout.pos.col = 2)); grid.draw(p_cd4_rna$gtable); popViewport()
  dev.off()
}

if (!is.null(p_cd4_adt)) {
  png(file.path(out_heat, "CD4_ADT_Heatmap.png"),
      width = 11, height = 10, units = "in", res = 300, bg = "white")
  grid.draw(p_cd4_adt$gtable)
  dev.off()
}

if (!is.null(p_cd4_rna)) {
  png(file.path(out_heat, "CD4_RNA_Heatmap.png"),
      width = 11, height = 12, units = "in", res = 300, bg = "white")
  grid.draw(p_cd4_rna$gtable)
  dev.off()
}

message("✓ CD4 heatmaps saved")

# ══════════════════════════════════════════════════════════════════════════════
# 3B) CD8 Heatmaps (from NKCD8 object — subset to CD8 clusters)
# ══════════════════════════════════════════════════════════════════════════════

# Subset to CD8 + gamma-delta + innate-like T cells (exclude NK)
cd8_labels <- c(
  "Activated cytotoxic CD8+", "CD8+ Memory", "CD8+ Naïve",
  "CD8+ Tem", "CD8+ Tem (PRF1+)", "CD8+ TEMRA",
  "NK-like CD8+", "ɣδ T", "Innate-like T"
)
OPIS_CD8 <- subset(OPIS_NKCD8, celltype_annotation %in% cd8_labels)

cd8_rna_markers <- c(
  # Naïve / Memory
  "TCF7"    = "Naïve/Memory",
  "CCR7"    = "Naïve/Memory",
  "IL7R"    = "Naïve/Memory",
  "S100A4"  = "Naïve/Memory",
  # Effector / Cytotoxic
  "PRF1"    = "Effector/Cytotoxic",
  "GNLY"    = "Effector/Cytotoxic",
  "NKG7"    = "Effector/Cytotoxic",
  "GZMK"    = "Effector/Cytotoxic",
  # Terminal Differentiation
  "KLRG1"   = "Term Diff",
  "FGFBP2"  = "Term Diff",
  # Activation
  "HLA-DRA" = "Activation",
  "IL2RA"   = "Activation",
  "FOS"     = "Activation",
  "JUN"     = "Activation",
  # NK-like Receptor Program
  "KLRD1"   = "NK-like Receptor Program",
  "KLRC1"   = "NK-like Receptor Program",
  # Gamma Delta T Identity
  "TRDC"    = "Gamma Delta T Identity",
  "TRDV1"   = "Gamma Delta T Identity",
  "TRGC1"   = "Gamma Delta T Identity",
  # Innate-like T Identity
  "KLRB1"   = "Innate-like T Identity",
  "ZBTB16"  = "Innate-like T Identity",
  "IL18RAP" = "Innate-like T Identity"
)

cd8_adt_markers <- c(
  # Naïve / Memory
  "CD45RA"  = "Naïve/Memory",
  "CD27"    = "Naïve/Memory",
  "CD45RO"  = "Naïve/Memory",
  # Terminal Differentiation
  "KLRG1"   = "Term Diff",
  # Activation
  "CD38"    = "Activation",
  "HLA_DR"  = "Activation",
  # NK-like Receptor Program
  "FCGR3A"  = "NK-like Receptor Program",
  "NCAM1"   = "NK-like Receptor Program",
  # Gamma Delta T Identity
  "TCR-AB"  = "Gamma Delta T Identity",
  "CD3"     = "Gamma Delta T Identity"
)

p_cd8_rna <- build_annotation_heatmap(
  OPIS_CD8, "RNA", cd8_rna_markers,
  title = "CD8+ T Cells — mRNA (RNA)"
)

p_cd8_adt <- build_annotation_heatmap(
  OPIS_CD8, "ADT", cd8_adt_markers,
  title = "CD8+ T Cells — Protein (ADT)"
)

# Save CD8 combined
if (!is.null(p_cd8_adt) && !is.null(p_cd8_rna)) {
  png(file.path(out_heat, "CD8_ADT_RNA_Heatmap_Combined.png"),
      width = 24, height = 16, units = "in", res = 300, bg = "white")
  grid.newpage()
  grid.text("CD8+ T Cell Annotation Validation: Protein (ADT) & mRNA",
            x = 0.5, y = 0.975,
            gp = gpar(fontsize = 18, fontface = "bold"))
  pushViewport(viewport(layout = grid.layout(1, 2)))
  pushViewport(viewport(layout.pos.col = 1)); grid.draw(p_cd8_adt$gtable); popViewport()
  pushViewport(viewport(layout.pos.col = 2)); grid.draw(p_cd8_rna$gtable); popViewport()
  dev.off()
}

if (!is.null(p_cd8_adt)) {
  png(file.path(out_heat, "CD8_ADT_Heatmap.png"),
      width = 11, height = 10, units = "in", res = 300, bg = "white")
  grid.draw(p_cd8_adt$gtable)
  dev.off()
}

if (!is.null(p_cd8_rna)) {
  png(file.path(out_heat, "CD8_RNA_Heatmap.png"),
      width = 11, height = 14, units = "in", res = 300, bg = "white")
  grid.draw(p_cd8_rna$gtable)
  dev.off()
}

message("✓ CD8 heatmaps saved")

# ══════════════════════════════════════════════════════════════════════════════
# 3C) NK Heatmaps (from NKCD8 object — subset to NK clusters)
# ══════════════════════════════════════════════════════════════════════════════

nk_labels <- grep("NK|Natural",
                  levels(factor(OPIS_NKCD8$celltype_annotation)),
                  value = TRUE, ignore.case = TRUE)
OPIS_NK <- subset(OPIS_NKCD8, celltype_annotation %in% nk_labels)

nk_rna_markers <- c(
  # NK Identity
  "NCAM1"   = "NK Identity",
  "FCGR3A"  = "NK Identity",
  "NKG7"    = "NK Identity",
  "GNLY"    = "NK Identity",
  "PRF1"    = "NK Identity",
  "KLRD1"   = "NK Identity",
  "FCER1G"  = "NK Identity",
  "TYROBP"  = "NK Identity",
  # CD56+ Less Mature
  "KLRC1"   = "CD56+ Less Mature",
  "XCL1"    = "CD56+ Less Mature",
  "XCL2"    = "CD56+ Less Mature"
)

nk_adt_markers <- c(
  # NK Identity
  "NCAM1"   = "NK Identity",
  "FCGR3A"  = "NK Identity"
)

p_nk_rna <- build_annotation_heatmap(
  OPIS_NK, "RNA", nk_rna_markers,
  title = "NK Cells — mRNA (RNA)"
)

p_nk_adt <- build_annotation_heatmap(
  OPIS_NK, "ADT", nk_adt_markers,
  title = "NK Cells — Protein (ADT)"
)

# Save NK combined
if (!is.null(p_nk_adt) && !is.null(p_nk_rna)) {
  png(file.path(out_heat, "NK_ADT_RNA_Heatmap_Combined.png"),
      width = 22, height = 12, units = "in", res = 300, bg = "white")
  grid.newpage()
  grid.text("NK Cell Annotation Validation: Protein (ADT) & mRNA",
            x = 0.5, y = 0.975,
            gp = gpar(fontsize = 18, fontface = "bold"))
  pushViewport(viewport(layout = grid.layout(1, 2)))
  pushViewport(viewport(layout.pos.col = 1)); grid.draw(p_nk_adt$gtable); popViewport()
  pushViewport(viewport(layout.pos.col = 2)); grid.draw(p_nk_rna$gtable); popViewport()
  dev.off()
}

if (!is.null(p_nk_rna)) {
  png(file.path(out_heat, "NK_RNA_Heatmap.png"),
      width = 11, height = 10, units = "in", res = 300, bg = "white")
  grid.draw(p_nk_rna$gtable)
  dev.off()
}

message("✓ NK heatmaps saved")

################################################################################
# 4) DGE — Cluster-wise MAST: OUD+ vs OUD−
################################################################################

safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)
count_sig <- function(df, pcol = "p_val_adj", thr = 0.05) sum(df[[pcol]] < thr, na.rm = TRUE)

run_clusterwise_de <- function(obj,
                               assay,
                               group_col,
                               ident1,
                               ident2,
                               cluster_col,
                               out_dir,
                               latent  = c("nCount_RNA"),
                               min_cells_per_grp = 10,
                               title_prefix = "",
                               file_stub = NULL) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  DefaultAssay(obj) <- assay
  Idents(obj) <- cluster_col
  
  de_list   <- list()
  count_tbl <- data.frame()
  
  cl_levels <- levels(Idents(obj))
  
  for (cl in cl_levels) {
    message("[", assay, "] Processing cluster: ", cl)
    
    cl_cells <- WhichCells(obj, idents = cl)
    sc <- subset(obj, cells = cl_cells)
    sc$.__grp <- sc[[group_col, drop = TRUE]]
    
    tab <- table(sc$.__grp)
    
    if (all(c(ident1, ident2) %in% names(tab)) &&
        all(tab[c(ident1, ident2)] >= min_cells_per_grp)) {
      
      de <- FindMarkers(
        sc,
        ident.1     = ident1,
        ident.2     = ident2,
        group.by    = ".__grp",
        test.use    = "MAST",
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
    } else {
      message("  ⚠ Skipping ", cl, ": insufficient cells in one or both groups")
    }
  }
  
  # Barplot of DE gene counts per cluster
  if (nrow(count_tbl) > 0) {
    p <- ggplot(count_tbl, aes(x = reorder(Cluster, -DE_Genes), y = DE_Genes, fill = Cluster)) +
      geom_bar(stat = "identity") +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none") +
      xlab("Cluster") +
      ylab("# DE Features (adj. p < 0.05)") +
      ggtitle(paste0(title_prefix, ": DE Features per Cluster (",
                     ident1, " vs ", ident2, " | ", assay, ")"))
    
    ggsave(
      filename = file.path(
        out_dir,
        paste0("Cluster_DE_Counts_",
               if (is.null(file_stub)) paste0(ident1, "vs", ident2) else file_stub,
               "_", assay, ".png")
      ),
      plot   = p,
      width  = 12,
      height = 6,
      dpi    = 400
    )
  }
  
  invisible(list(de = de_list, counts = count_tbl))
}

# ── Prepare object for DGE ────────────────────────────────────────────────────
# Use the broad_celltype-annotated OPIS_ALL
OPIS_ALL$OUD_status <- factor(OPIS_ALL$OUD_status)
OPIS_OUD <- subset(OPIS_ALL, subset = !is.na(OUD_status))

message("OUD_status levels: ", paste(levels(OPIS_OUD$OUD_status), collapse = ", "))

ident1_OUD <- "OUD+"
ident2_OUD <- "OUD-"


# ── 4B) DGE on CD4 subtypes ──────────────────────────────────────────────────
# Subset OUD object to CD4 cells and use fine annotations
OPIS_CD4_OUD <- subset(OPIS_CD4, subset = !is.na(OUD_status))
OPIS_CD4_OUD$OUD_status <- factor(OPIS_CD4_OUD$OUD_status)

dge_cd4_rna_dir <- file.path(out_dge, "CD4_Subtypes", "RNA")

dge_cd4_rna <- run_clusterwise_de(
  obj          = OPIS_CD4_OUD,
  assay        = "RNA",
  group_col    = "OUD_status",
  ident1       = ident1_OUD,
  ident2       = ident2_OUD,
  cluster_col  = "celltype_annotation",
  out_dir      = dge_cd4_rna_dir,
  latent       = c("nCount_RNA"),
  title_prefix = "OPIS CD4 Subtypes (RNA)",
  file_stub    = "OUD_Pos_vs_Neg_RNA"
)

dge_cd4_adt_dir <- file.path(out_dge, "CD4_Subtypes", "ADT")

dge_cd4_adt <- run_clusterwise_de(
  obj          = OPIS_CD4_OUD,
  assay        = "ADT",
  group_col    = "OUD_status",
  ident1       = ident1_OUD,
  ident2       = ident2_OUD,
  cluster_col  = "celltype_annotation",
  out_dir      = dge_cd4_adt_dir,
  latent       = c("nCount_ADT"),
  title_prefix = "OPIS CD4 Subtypes (ADT)",
  file_stub    = "OUD_Pos_vs_Neg_ADT"
)

message("✓ DGE (CD4 subtypes, RNA + ADT) complete")

# ── 4C) DGE on NK/CD8 subtypes ───────────────────────────────────────────────
OPIS_NKCD8_OUD <- subset(OPIS_NKCD8, subset = !is.na(OUD_status))
OPIS_NKCD8_OUD$OUD_status <- factor(OPIS_NKCD8_OUD$OUD_status)

dge_nkcd8_rna_dir <- file.path(out_dge, "NKCD8_Subtypes", "RNA")

dge_nkcd8_rna <- run_clusterwise_de(
  obj          = OPIS_NKCD8_OUD,
  assay        = "RNA",
  group_col    = "OUD_status",
  ident1       = ident1_OUD,
  ident2       = ident2_OUD,
  cluster_col  = "celltype_annotation",
  out_dir      = dge_nkcd8_rna_dir,
  latent       = c("nCount_RNA"),
  title_prefix = "OPIS NK/CD8 Subtypes (RNA)",
  file_stub    = "OUD_Pos_vs_Neg_RNA"
)

dge_nkcd8_adt_dir <- file.path(out_dge, "NKCD8_Subtypes", "ADT")

dge_nkcd8_adt <- run_clusterwise_de(
  obj          = OPIS_NKCD8_OUD,
  assay        = "ADT",
  group_col    = "OUD_status",
  ident1       = ident1_OUD,
  ident2       = ident2_OUD,
  cluster_col  = "celltype_annotation",
  out_dir      = dge_nkcd8_adt_dir,
  latent       = c("nCount_ADT"),
  title_prefix = "OPIS NK/CD8 Subtypes (ADT)",
  file_stub    = "OUD_Pos_vs_Neg_ADT"
)

message("✓ DGE (NK/CD8 subtypes, RNA + ADT) complete")

