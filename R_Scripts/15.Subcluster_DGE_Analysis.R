################################################################################
# OPIS ECCITE-seq — Manuscript Figure Generation & DGE/DPE Analysis (UPDATED)
#
# Changes from previous version:
#   1A) Broad-cluster DGE (RNA) & DPE (ADT) on OPIS_ALL: OUD+ vs OUD−
#       → Volcano plots, top-30 heatmaps (lavender/violet, z-scored, 2-col: OUD+/OUD−)
#       → All plots generated from saved CSVs (re-runnable without MAST)
#   1B) Subset-level DGE/DPE remains on CD4 & NK/CD8 annotated objects
#   2)  NK heatmap: exclude "NK-like CD8+" from NK subset
#   3)  CD8 ADT heatmap: add CD8A under "Lineage"
#       CD8 RNA heatmap: add CD3D, TRAC under "Lineage"
#   4)  Feature plots & split-by-OUD violin plots for CD4 & NK/CD8 subsets
#       → Organized folder structure under Manuscript_Figures
#   5)  EnrichR pathway analysis on subset-level DGE/DPE results
#       → Publication-quality plots
#
# Output directory tree:
#   Manuscript_Figures/
#   ├── UMAP/
#   ├── Heatmaps/
#   ├── DGE/
#   │   ├── Broad_Clusters/
#   │   │   ├── RNA/
#   │   │   └── ADT/
#   │   ├── CD4_Subtypes/
#   │   │   ├── RNA/
#   │   │   └── ADT/
#   │   └── NKCD8_Subtypes/
#   │       ├── RNA/
#   │       └── ADT/
#   ├── Feature_Plots/
#   │   ├── CD4/
#   │   │   ├── RNA/
#   ├── DGE/
#   │   ├── Broad_Clusters/
#   │   │   ├── RNA/
#   │   │   │   ├── CSVs/
#   │   │   │   └── Plots/
#   │   │   │       ├── Volcanos/
#   │   │   │       └── Heatmaps/
#   │   │   └── ADT/
#   │   │       ├── CSVs/
#   │   │       └── Plots/
#   │   │           ├── Volcanos/
#   │   │           └── Heatmaps/
#   │   ├── CD4_Subtypes/
#   │   │   ├── RNA/  (same CSVs/Plots structure)
#   │   │   └── ADT/
#   │   └── NKCD8_Subtypes/
#   │       ├── RNA/
#   │       └── ADT/
#   ├── Feature_Plots/
#   │   ├── CD4/
#   │   │   ├── RNA/
#   │   │   └── ADT/
#   │   └── NKCD8/
#   │       ├── RNA/
#   │       └── ADT/
#   ├── Violin_Plots/
#   │   ├── CD4/
#   │   │   ├── RNA_splitBy_OUD/
#   │   │   ├── RNA/
#   │   │   ├── ADT_splitBy_OUD/
#   │   │   └── ADT/
#   │   └── NKCD8/
#   │       ├── RNA_splitBy_OUD/
#   │       ├── RNA/
#   │       ├── ADT_splitBy_OUD/
#   │       └── ADT/
#   └── Pathway_Analysis/
#       ├── CD4_Subtypes/
#       │   ├── RNA/
#       │   └── ADT/
#       └── NKCD8_Subtypes/
#           ├── RNA/
#           └── ADT/
################################################################################

library(Seurat)
library(SeuratObject)
library(qs2)
library(ggplot2)
library(pheatmap)
library(grid)
library(dplyr)
library(scCustomize)
library(SeuratExtend)
library(viridis)
library(enrichR)
library(openxlsx)
library(ggrepel)
library(patchwork)

# ── Global ggplot theme: white background ─────────────────────────────────────
theme_set(theme_classic(base_size = 14) + theme(
  panel.background = element_rect(fill = "white", colour = NA),
  plot.background  = element_rect(fill = "white", colour = NA)
))

# ── Paths ─────────────────────────────────────────────────────────────────────
load.path <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/saved_R_data/"
out_base  <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/Manuscript_Figures"

# Create full directory tree
out_umap <- file.path(out_base, "UMAP")
out_heat <- file.path(out_base, "Heatmaps")

# DGE directories (each gets CSVs/ and Plots/Volcanos/, Plots/Heatmaps/ subfolders)
out_dge_broad_rna   <- file.path(out_base, "DGE", "Broad_Clusters", "RNA")
out_dge_broad_adt   <- file.path(out_base, "DGE", "Broad_Clusters", "ADT")
out_dge_cd4_rna     <- file.path(out_base, "DGE", "CD4_Subtypes", "RNA")
out_dge_cd4_adt     <- file.path(out_base, "DGE", "CD4_Subtypes", "ADT")
out_dge_nkcd8_rna   <- file.path(out_base, "DGE", "NKCD8_Subtypes", "RNA")
out_dge_nkcd8_adt   <- file.path(out_base, "DGE", "NKCD8_Subtypes", "ADT")

# Pre-create DGE subfolder structure
for (dge_dir in c(out_dge_broad_rna, out_dge_broad_adt,
                  out_dge_cd4_rna, out_dge_cd4_adt,
                  out_dge_nkcd8_rna, out_dge_nkcd8_adt)) {
  dir.create(file.path(dge_dir, "CSVs"),             recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dge_dir, "Plots", "Volcanos"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dge_dir, "Plots", "Heatmaps"), recursive = TRUE, showWarnings = FALSE)
}

# Feature plot directories
out_feat_cd4_rna    <- file.path(out_base, "Feature_Plots", "CD4", "RNA")
out_feat_cd4_adt    <- file.path(out_base, "Feature_Plots", "CD4", "ADT")
out_feat_nkcd8_rna  <- file.path(out_base, "Feature_Plots", "NKCD8", "RNA")
out_feat_nkcd8_adt  <- file.path(out_base, "Feature_Plots", "NKCD8", "ADT")

# Violin plot directories
out_vln_cd4_rna_oud     <- file.path(out_base, "Violin_Plots", "CD4", "RNA_splitBy_OUD")
out_vln_cd4_rna         <- file.path(out_base, "Violin_Plots", "CD4", "RNA")
out_vln_cd4_adt_oud     <- file.path(out_base, "Violin_Plots", "CD4", "ADT_splitBy_OUD")
out_vln_cd4_adt         <- file.path(out_base, "Violin_Plots", "CD4", "ADT")
out_vln_nkcd8_rna_oud   <- file.path(out_base, "Violin_Plots", "NKCD8", "RNA_splitBy_OUD")
out_vln_nkcd8_rna       <- file.path(out_base, "Violin_Plots", "NKCD8", "RNA")
out_vln_nkcd8_adt_oud   <- file.path(out_base, "Violin_Plots", "NKCD8", "ADT_splitBy_OUD")
out_vln_nkcd8_adt       <- file.path(out_base, "Violin_Plots", "NKCD8", "ADT")

# Pathway analysis directories
out_path_cd4_rna    <- file.path(out_base, "Pathway_Analysis", "CD4_Subtypes", "RNA")
out_path_cd4_adt    <- file.path(out_base, "Pathway_Analysis", "CD4_Subtypes", "ADT")
out_path_nkcd8_rna  <- file.path(out_base, "Pathway_Analysis", "NKCD8_Subtypes", "RNA")
out_path_nkcd8_adt  <- file.path(out_base, "Pathway_Analysis", "NKCD8_Subtypes", "ADT")

# Create all directories
all_dirs <- c(
  out_umap, out_heat,
  out_feat_cd4_rna, out_feat_cd4_adt,
  out_feat_nkcd8_rna, out_feat_nkcd8_adt,
  out_vln_cd4_rna_oud, out_vln_cd4_rna,
  out_vln_cd4_adt_oud, out_vln_cd4_adt,
  out_vln_nkcd8_rna_oud, out_vln_nkcd8_rna,
  out_vln_nkcd8_adt_oud, out_vln_nkcd8_adt,
  out_path_cd4_rna, out_path_cd4_adt,
  out_path_nkcd8_rna, out_path_nkcd8_adt
)
lapply(all_dirs, function(x) dir.create(x, recursive = TRUE, showWarnings = FALSE))


################################################################################
# SHARED UTILITIES
################################################################################

safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)
count_sig <- function(df, pcol = "p_val_adj", thr = 0.05) sum(df[[pcol]] < thr, na.rm = TRUE)

# ── Heatmap utilities ─────────────────────────────────────────────────────────
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
build_annotation_heatmap <- function(obj, assay_name, marker_groups,
                                     group_col = "celltype_annotation",
                                     title = "", slot_name = "data") {
  DefaultAssay(obj) <- assay_name
  available <- rownames(obj)
  keep      <- names(marker_groups) %in% available
  if (sum(keep) == 0) {
    message("  ⚠ No features found for ", assay_name, " in: ", title)
    return(NULL)
  }
  marker_groups <- marker_groups[keep]
  
  avg <- AverageExpression(
    obj, assays = assay_name, features = names(marker_groups),
    group.by = group_col, slot = slot_name
  )[[assay_name]]
  
  if (assay_name == "ADT") {
    avg_scaled <- scale_01(log10(avg + 1))
  } else {
    avg_scaled <- scale_01(avg)
  }
  colnames(avg_scaled) <- gsub("^g ", "", colnames(avg_scaled))
  
  # Group colors
  base_palette <- c(
    "Naïve/Memory" = "#74C2E1", "Naïve / Memory" = "#74C2E1",
    "Naïve/Stemness" = "#74C2E1", "Activation" = "#FFD54F",
    "Regulatory" = "#AB47BC", "Effector/Term Diff" = "#EF5350",
    "Effector/Cytotoxic" = "#EF5350", "Term Diff" = "#FF7043",
    "Lineage" = "#78909C", "NK Identity" = "#F06292",
    "NK-like Receptor Program" = "#F06292", "CD56+ Less Mature" = "#CE93D8",
    "Gamma Delta T Identity" = "#9C6FD6", "Innate-like T Identity" = "#FFAB91",
    "Memory/Stemness" = "#66BB6A", "Co-stimulation" = "#AED581",
    "Effector/TEMRA" = "#EF5350", "Effector TFs" = "#FFA726",
    "NK-like/Innate" = "#F06292", "Exhaustion" = "#A1887F"
  )
  all_groups   <- unique(marker_groups)
  group_colors <- setNames(
    sapply(all_groups, function(g) if (g %in% names(base_palette)) base_palette[g] else "#BDBDBD"),
    all_groups
  )
  
  annot_df <- data.frame(
    Group = factor(marker_groups[rownames(avg_scaled)], levels = unique(marker_groups)),
    row.names = rownames(avg_scaled)
  )
  
  row_order     <- cluster_within_groups(avg_scaled, annot_df)
  avg_ordered   <- avg_scaled[row_order, ]
  annot_ordered <- annot_df[row_order, , drop = FALSE]
  gaps          <- get_gaps(annot_ordered)
  annot_colors  <- list(Group = group_colors[levels(droplevels(annot_ordered$Group))])
  
  p <- pheatmap(
    avg_ordered, cluster_rows = FALSE, cluster_cols = FALSE,
    scale = "none", color = heatmap_colors, border_color = "white",
    annotation_row = annot_ordered, annotation_colors = annot_colors,
    annotation_names_row = FALSE, gaps_row = gaps,
    cellwidth = 55, cellheight = 16, fontsize = 11,
    fontsize_row = 10, fontsize_col = 10, angle_col = 45,
    main = title, silent = TRUE
  )
  return(p)
}

# ── Cluster-wise DE function ──────────────────────────────────────────────────
run_clusterwise_de <- function(obj, assay, group_col, ident1, ident2,
                               cluster_col, out_dir,
                               latent = c("nCount_RNA"),
                               min_cells_per_grp = 10,
                               title_prefix = "", file_stub = NULL) {
  # Subfolder structure
  csv_dir  <- file.path(out_dir, "CSVs")
  plot_dir <- file.path(out_dir, "Plots")
  dir.create(csv_dir,  recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  DefaultAssay(obj) <- assay
  Idents(obj)       <- cluster_col
  
  de_list   <- list()
  count_tbl <- data.frame()
  cl_levels <- levels(Idents(obj))
  
  for (cl in cl_levels) {
    message("[", assay, "] Processing cluster: ", cl)
    cl_cells <- WhichCells(obj, idents = cl)
    sc       <- subset(obj, cells = cl_cells)
    sc$.__grp <- sc[[group_col, drop = TRUE]]
    tab <- table(sc$.__grp)
    
    if (all(c(ident1, ident2) %in% names(tab)) &&
        all(tab[c(ident1, ident2)] >= min_cells_per_grp)) {
      
      de <- FindMarkers(
        sc, ident.1 = ident1, ident.2 = ident2,
        group.by = ".__grp", test.use = "MAST", latent.vars = latent
      )
      comp_name <- if (is.null(file_stub)) paste0(ident1, "_vs_", ident2) else file_stub
      out_csv   <- file.path(csv_dir, paste0(safe_name(cl), "_", comp_name, ".csv"))
      write.csv(de, out_csv)
      de_list[[cl]] <- de
      count_tbl <- rbind(count_tbl, data.frame(Cluster = cl, DE_Genes = count_sig(de)))
    } else {
      message("  ⚠ Skipping ", cl, ": insufficient cells in one or both groups")
    }
  }
  
  # ── Summary bar plot of DE gene counts (lavender palette) ───────────────────
  if (nrow(count_tbl) > 0) {
    p <- ggplot(count_tbl, aes(x = reorder(Cluster, -DE_Genes), y = DE_Genes)) +
      geom_bar(stat = "identity", fill = "#B39DDB", color = "#7E57C2", linewidth = 0.4) +
      theme_classic(base_size = 14) +
      theme(
        axis.text.x  = element_text(angle = 45, hjust = 1, size = 11, color = "black"),
        axis.text.y  = element_text(size = 11, color = "black"),
        axis.title   = element_text(size = 13),
        plot.title   = element_text(size = 15, face = "bold", hjust = 0.5),
        panel.background = element_rect(fill = "white"),
        plot.background  = element_rect(fill = "white")
      ) +
      xlab("Cluster") +
      ylab("# DE Features (adj. p < 0.05)") +
      ggtitle(paste0(title_prefix, ": DE Features per Cluster (",
                     ident1, " vs ", ident2, " | ", assay, ")"))
    
    ggsave(
      filename = file.path(plot_dir, paste0("Cluster_DE_Counts_",
                                            if (is.null(file_stub)) paste0(ident1, "vs", ident2) else file_stub,
                                            "_", assay, ".png")),
      plot = p, width = 12, height = 6, dpi = 400, bg = "white"
    )
  }
  
  invisible(list(de = de_list, counts = count_tbl))
}

# ── Volcano plot function ─────────────────────────────────────────────────────
make_volcano <- function(de_df, cluster_name, assay_label, top_n = 30, out_dir) {
  if (is.null(de_df) || nrow(de_df) == 0) return(invisible(NULL))
  
  de_df$gene <- rownames(de_df)
  de_df$neg_log10_padj <- -log10(de_df$p_val_adj + 1e-300)
  de_df$significance <- case_when(
    de_df$p_val_adj < 0.05 & de_df$avg_log2FC > 0.25  ~ "Up",
    de_df$p_val_adj < 0.05 & de_df$avg_log2FC < -0.25 ~ "Down",
    TRUE ~ "NS"
  )
  
  # Top genes by combined rank of FC and significance
  top_genes <- de_df %>%
    filter(significance != "NS") %>%
    arrange(p_val_adj) %>%
    slice_head(n = top_n)
  
  p <- ggplot(de_df, aes(x = avg_log2FC, y = neg_log10_padj, color = significance)) +
    geom_point(size = 1.2, alpha = 0.6) +
    scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1565C0", "NS" = "grey70"),
                       name = "Direction") +
    geom_text_repel(
      data = top_genes, aes(label = gene),
      size = 3, max.overlaps = 25, segment.size = 0.3, color = "black"
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey40") +
    theme_classic(base_size = 13) +
    theme(
      plot.title   = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = "right"
    ) +
    xlab("Average log2 Fold Change") +
    ylab("-log10(Adjusted P-value)") +
    ggtitle(paste0(cluster_name, " — ", assay_label, " (OUD+ vs OUD−)"))
  
  ggsave(
    file.path(out_dir, paste0(safe_name(cluster_name), "_Volcano_", assay_label, ".png")),
    plot = p, width = 10, height = 8, dpi = 400, bg = "white"
  )
}

# ── Lavender/violet palette for top-30 DE heatmaps ───────────────────────────
lavender_colors <- colorRampPalette(c(
  "#F3E5F5", "#CE93D8", "#AB47BC", "#8E24AA", "#4A148C"
))(100)

# ── Top-30 heatmap (2-column: OUD+ vs OUD−, z-scored, lavender palette) ──────
# Takes a DE data.frame (with avg_log2FC, pct.1, pct.2, p_val_adj columns).
# Constructs a 2-column matrix: z-scored expression proxy for OUD+ and OUD−.
# pct.1 = fraction expressing in ident1 (OUD+), pct.2 = in ident2 (OUD−).
# We use avg_log2FC to derive relative z-scored values per gene:
#   OUD+ column = +avg_log2FC/2, OUD− column = −avg_log2FC/2, then row-z-score.
make_top30_heatmap <- function(de_df, cluster_name, assay_label, out_dir) {
  if (is.null(de_df) || nrow(de_df) == 0) return(invisible(NULL))
  
  de_df$gene <- rownames(de_df)
  sig <- de_df %>%
    filter(p_val_adj < 0.05) %>%
    arrange(desc(abs(avg_log2FC))) %>%
    slice_head(n = 30)
  
  if (nrow(sig) == 0) return(invisible(NULL))
  
  # Build matrix: rows = genes, cols = OUD+ / OUD−
  # Use avg_log2FC to create a symmetric representation, then z-score rows
  mat <- matrix(NA, nrow = nrow(sig), ncol = 2,
                dimnames = list(sig$gene, c("OUD+", "OUD\u2212")))
  mat[, 1] <-  sig$avg_log2FC / 2
  mat[, 2] <- -sig$avg_log2FC / 2
  
  # Row-wise z-score
  mat_z <- t(apply(mat, 1, function(x) {
    s <- sd(x)
    if (is.na(s) || s == 0) return(rep(0, length(x)))
    (x - mean(x)) / s
  }))
  colnames(mat_z) <- colnames(mat)
  
  # Order by avg_log2FC (most upregulated in OUD+ at top)
  mat_z <- mat_z[order(sig$avg_log2FC, decreasing = TRUE), ]
  
  title_str <- paste0("Top 30 DE — ", cluster_name, " (", assay_label, ", OUD+ vs OUD\u2212)")
  
  p <- pheatmap(
    mat_z,
    cluster_rows  = FALSE,
    cluster_cols  = FALSE,
    scale         = "none",
    color         = lavender_colors,
    border_color  = "white",
    cellwidth     = 55,
    cellheight    = 16,
    fontsize      = 11,
    fontsize_row  = 10,
    fontsize_col  = 12,
    angle_col     = 0,
    main          = title_str,
    silent        = TRUE
  )
  
  out_file <- file.path(out_dir, paste0(safe_name(cluster_name), "_Top30_Heatmap_", assay_label, ".png"))
  png(out_file, width = 6, height = 12, units = "in", res = 400, bg = "white")
  grid.draw(p$gtable)
  dev.off()
  message("  Saved top-30 heatmap: ", out_file)
}

# ══════════════════════════════════════════════════════════════════════════════
# CSV-based plotting: read saved DGE/DPE CSVs and generate volcano + heatmap
# This lets you re-run plots without re-running MAST.
#
# base_dir should be e.g. .../DGE/Broad_Clusters/RNA
# Reads from:  base_dir/CSVs/
# Writes to:   base_dir/Plots/Volcanos/  and  base_dir/Plots/Heatmaps/
# ══════════════════════════════════════════════════════════════════════════════
make_plots_from_csvs <- function(base_dir, assay_label, top_n = 30) {
  csv_dir     <- file.path(base_dir, "CSVs")
  volcano_dir <- file.path(base_dir, "Plots", "Volcanos")
  heatmap_dir <- file.path(base_dir, "Plots", "Heatmaps")
  dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(heatmap_dir, recursive = TRUE, showWarnings = FALSE)
  
  csv_files <- list.files(csv_dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(csv_files) == 0) {
    message("  No CSVs found in: ", csv_dir)
    return(invisible(NULL))
  }
  
  for (csv_path in csv_files) {
    de_df <- tryCatch(
      read.csv(csv_path, row.names = 1, check.names = FALSE),
      error = function(e) { message("  Read failed: ", csv_path); return(NULL) }
    )
    if (is.null(de_df)) next
    if (!"p_val_adj" %in% colnames(de_df)) next
    
    # Derive cluster name from filename (strip file_stub suffix)
    fname <- tools::file_path_sans_ext(basename(csv_path))
    # Remove trailing _OUD_Pos_vs_Neg_RNA or _OUD_Pos_vs_Neg_ADT
    cl_name <- gsub("_OUD_Pos_vs_Neg_(RNA|ADT)$", "", fname)
    cl_name <- gsub("_", " ", cl_name)  # restore spaces
    
    message("  Plotting from CSV: ", basename(csv_path), " → cluster: ", cl_name)
    make_volcano(de_df, cl_name, assay_label, top_n = top_n, out_dir = volcano_dir)
    make_top30_heatmap(de_df, cl_name, assay_label, out_dir = heatmap_dir)
  }
}


################################################################################
# 1) UMAP PBMC — Combine clusters to broad cell types
################################################################################

OPIS_ALL <- qs_read(file.path(load.path, "OPIS_ALL_PreAnnotation.qs2"))
OPIS_ALL <- subset(OPIS_ALL, wsnn_res.0.4 %in% as.character(0:25))

cluster_to_celltype <- c(
  "5" = "CD4+", "6" = "CD4+", "7" = "CD4+", "8" = "CD4+",
  "0" = "CD8+", "2" = "CD8+", "9" = "CD8+", "10" = "CD8+",
  "11" = "CD8+", "13" = "CD8+", "14" = "CD8+", "18" = "CD8+", "19" = "CD8+",
  "4" = "NK", "16" = "NK",
  "3" = "B", "21" = "B", "24" = "B", "12" = "B",
  "17" = "Monocytes", "1" = "Monocytes", "22" = "Monocytes",
  "23" = "Monocytes", "15" = "Monocytes", "25" = "Monocytes",
  "20" = "pDC"
)

OPIS_ALL$broad_celltype <- unname(cluster_to_celltype[as.character(OPIS_ALL$wsnn_res.0.4)])
OPIS_ALL$broad_celltype <- factor(
  OPIS_ALL$broad_celltype,
  levels = c("CD4+", "CD8+", "NK", "B", "Monocytes", "pDC")
)

# Plot: original clusters
p1 <- DimPlot2(
  OPIS_ALL, reduction = "wnn.umap", group.by = "wsnn_res.0.4",
  label = TRUE, repel = TRUE, box = TRUE, label.size = 5
) + ggtitle("wsnn_res.0.4")

# Plot: broad cell types
p2 <- DimPlot2(
  OPIS_ALL, reduction = "wnn.umap", group.by = "broad_celltype",
  label = TRUE, repel = TRUE, box = TRUE, label.size = 8,
  label.color = "black",
  theme = list(NoLegend(), NoAxes(), theme_umap_arrows())
) + ggtitle("Broad Cell Type Annotation")

# Save combined
png(file.path(out_umap, "PBMC_UMAP_Clusters_and_CellTypes.png"),
    width = 20, height = 8, units = "in", res = 300, bg = "white")
print(p1 | p2)
dev.off()

ggsave(file.path(out_umap, "PBMC_UMAP_Clusters.png"),
       plot = p1, width = 10, height = 8, dpi = 300)
ggsave(file.path(out_umap, "PBMC_UMAP_BroadCellTypes.png"),
       plot = p2, width = 10, height = 8, dpi = 300)

ClusterDistrBar(origin = OPIS_ALL$OUD_status, cluster = OPIS_ALL$broad_celltype,
                flip = FALSE, reverse_order = FALSE)
ggsave(file.path(out_umap, "PBMC_Cluster_Prop_Stacked_Bar.png"),
       width = 6, height = 5, bg = "white")

ClusterDistrPlot(
  origin = OPIS_ALL$orig.ident, cluster = OPIS_ALL$broad_celltype,
  condition = OPIS_ALL$OUD_status, cols = c('#1180808B', '#F54927')
)
ggsave(file.path(out_umap, "PBMC_Clust_Prop_wilcoxins.png"),
       width = 10, height = 6, bg = "white")

message("✓ Section 1: PBMC UMAP saved")

################################################################################
# 1A) NEW: Broad-cluster DGE (RNA) & DPE (ADT) — OUD+ vs OUD−
#     Volcano + Top-30 z-scored heatmaps (lavender palette)
#     Plots are generated from saved CSVs — re-runnable without MAST
################################################################################

OPIS_ALL$OUD_status <- factor(OPIS_ALL$OUD_status)
OPIS_OUD <- subset(OPIS_ALL, subset = !is.na(OUD_status))

ident1_OUD <- "OUD+"
ident2_OUD <- "OUD-"

message("=== 1A: Broad-cluster DGE (RNA) on OPIS_ALL ===")
dge_broad_rna <- run_clusterwise_de(
  
  obj          = OPIS_OUD,
  assay        = "RNA",
  group_col    = "OUD_status",
  ident1       = ident1_OUD,
  ident2       = ident2_OUD,
  cluster_col  = "broad_celltype",
  out_dir      = out_dge_broad_rna,
  latent       = c("nCount_RNA"),
  title_prefix = "OPIS Broad Clusters (RNA)",
  file_stub    = "OUD_Pos_vs_Neg_RNA"
)

# Generate volcano + top-30 heatmaps from saved CSVs (RNA broad clusters)
make_plots_from_csvs(out_dge_broad_rna, "RNA", top_n = 30)

message("=== 1A: Broad-cluster DPE (ADT) on OPIS_ALL ===")
dge_broad_adt <- run_clusterwise_de(
  obj          = OPIS_OUD,
  assay        = "ADT",
  group_col    = "OUD_status",
  ident1       = ident1_OUD,
  ident2       = ident2_OUD,
  cluster_col  = "broad_celltype",
  out_dir      = out_dge_broad_adt,
  latent       = c("nCount_ADT"),
  title_prefix = "OPIS Broad Clusters (ADT)",
  file_stub    = "OUD_Pos_vs_Neg_ADT"
)

# Generate volcano + top-30 heatmaps from saved CSVs (ADT broad clusters)
make_plots_from_csvs(out_dge_broad_adt, "ADT", top_n = 30)

message("✓ Section 1A: Broad-cluster DGE/DPE complete")

################################################################################
# 2) Load CD4 & NK/CD8 annotated objects + DimPlots
################################################################################

OPIS_CD4   <- qs_read(file.path(load.path, "OPIS_CD4_Annotated.qs2"))
OPIS_NKCD8 <- qs_read(file.path(load.path, "OPIS_NKCD8_Annotated.qs2"))

stopifnot("celltype_annotation" %in% colnames(OPIS_CD4@meta.data))
stopifnot("celltype_annotation" %in% colnames(OPIS_NKCD8@meta.data))

# CD4 UMAP
p_cd4_umap <- DimPlot2(
  OPIS_CD4, reduction = "wnn.umap", group.by = "celltype_annotation",
  label = TRUE, repel = TRUE, box = TRUE, pt.size = 1.5, label.size = 7,
  theme = list(NoLegend(), NoAxes(), theme_umap_arrows())
) + ggtitle("CD4+ T Cell Subtypes")

ggsave(file.path(out_umap, "CD4_UMAP_Annotated.png"),
       plot = p_cd4_umap, width = 10, height = 8, dpi = 300)

# NK/CD8 UMAP
p_nkcd8_umap <- DimPlot2(
  OPIS_NKCD8, reduction = "wnn.umap", group.by = "celltype_annotation",
  label = TRUE, repel = TRUE, box = TRUE, pt.size = 1.5, label.size = 7,
  theme = list(NoLegend(), NoAxes(), theme_umap_arrows())
) + ggtitle("NK & CD8+ T Cell Subtypes")

ggsave(file.path(out_umap, "NKCD8_UMAP_Annotated.png"),
       plot = p_nkcd8_umap, width = 10, height = 8, dpi = 300)

message("✓ Section 2: CD4 & NK/CD8 UMAPs saved")

################################################################################
# 3) Annotation Heatmaps — CD4, CD8, NK
################################################################################

# ══════════════════════════════════════════════════════════════════════════════
# 3A) CD4 Heatmaps (unchanged)
# ══════════════════════════════════════════════════════════════════════════════

cd4_rna_markers <- c(
  "TCF7" = "Naïve/Memory", "LEF1" = "Naïve/Memory", "BACH2" = "Naïve/Memory",
  "CCR7" = "Naïve/Memory", "IL7R" = "Naïve/Memory", "S100A4" = "Naïve/Memory",
  "HLA-DRA" = "Activation", "IL2RA" = "Activation", "FOS" = "Activation", "JUN" = "Activation",
  "FOXP3" = "Regulatory", "TIGIT" = "Regulatory", "CTLA4" = "Regulatory",
  "KLRG1" = "Effector/Term Diff", "PRF1" = "Effector/Term Diff", "GZMK" = "Effector/Term Diff",
  "CD3D" = "Lineage", "TRAC" = "Lineage"
)

cd4_adt_markers <- c(
  "CD45RA" = "Naïve/Memory", "CD27" = "Naïve/Memory", "SELL" = "Naïve/Memory", "CD45RO" = "Naïve/Memory",
  "CD38" = "Activation", "ICOS" = "Activation", "HLA_DR" = "Activation",
  "CD25" = "Regulatory", "CCR4" = "Regulatory",
  "KLRG1" = "Effector/Term Diff",
  "CD4" = "Lineage", "CD8" = "Lineage"
)

p_cd4_rna <- build_annotation_heatmap(OPIS_CD4, "RNA", cd4_rna_markers, title = "CD4+ T Cells — mRNA (RNA)")
p_cd4_adt <- build_annotation_heatmap(OPIS_CD4, "ADT", cd4_adt_markers, title = "CD4+ T Cells — Protein (ADT)")

if (!is.null(p_cd4_adt) && !is.null(p_cd4_rna)) {
  png(file.path(out_heat, "CD4_ADT_RNA_Heatmap_Combined.png"),
      width = 24, height = 14, units = "in", res = 300, bg = "white")
  grid.newpage()
  grid.text("CD4+ T Cell Annotation Validation: Protein (ADT) & mRNA",
            x = 0.5, y = 0.975, gp = gpar(fontsize = 18, fontface = "bold"))
  pushViewport(viewport(layout = grid.layout(1, 2)))
  pushViewport(viewport(layout.pos.col = 1)); grid.draw(p_cd4_adt$gtable); popViewport()
  pushViewport(viewport(layout.pos.col = 2)); grid.draw(p_cd4_rna$gtable); popViewport()
  dev.off()
}
if (!is.null(p_cd4_adt)) {
  png(file.path(out_heat, "CD4_ADT_Heatmap.png"), width = 11, height = 10, units = "in", res = 300, bg = "white")
  grid.draw(p_cd4_adt$gtable); dev.off()
}
if (!is.null(p_cd4_rna)) {
  png(file.path(out_heat, "CD4_RNA_Heatmap.png"), width = 11, height = 12, units = "in", res = 300, bg = "white")
  grid.draw(p_cd4_rna$gtable); dev.off()
}
message("✓ CD4 heatmaps saved")

# ══════════════════════════════════════════════════════════════════════════════
# 3B) CD8 Heatmaps — UPDATED: added CD8A (ADT Lineage), CD3D+TRAC (RNA Lineage)
# ══════════════════════════════════════════════════════════════════════════════

cd8_labels <- c(
  "Activated cytotoxic CD8+", "CD8+ Memory", "CD8+ Naïve",
  "CD8+ Tem", "CD8+ Tem (PRF1+)", "CD8+ TEMRA",
  "NK-like CD8+", "ɣδ T", "Innate-like T"
)
OPIS_CD8 <- subset(OPIS_NKCD8, celltype_annotation %in% cd8_labels)

# UPDATED: added CD3D, TRAC under Lineage
cd8_rna_markers <- c(
  "TCF7" = "Naïve/Memory", "CCR7" = "Naïve/Memory", "IL7R" = "Naïve/Memory", "S100A4" = "Naïve/Memory",
  "PRF1" = "Effector/Cytotoxic", "GNLY" = "Effector/Cytotoxic", "NKG7" = "Effector/Cytotoxic", "GZMK" = "Effector/Cytotoxic",
  "KLRG1" = "Term Diff", "FGFBP2" = "Term Diff",
  "HLA-DRA" = "Activation", "IL2RA" = "Activation", "FOS" = "Activation", "JUN" = "Activation",
  "KLRD1" = "NK-like Receptor Program", "KLRC1" = "NK-like Receptor Program",
  "TRDC" = "Gamma Delta T Identity", "TRDV1" = "Gamma Delta T Identity", "TRGC1" = "Gamma Delta T Identity",
  "KLRB1" = "Innate-like T Identity", "ZBTB16" = "Innate-like T Identity", "IL18RAP" = "Innate-like T Identity",
  "CD3D" = "Lineage", "TRAC" = "Lineage"
)

# UPDATED: added CD8A under Lineage
cd8_adt_markers <- c(
  "CD45RA" = "Naïve/Memory", "CD27" = "Naïve/Memory", "CD45RO" = "Naïve/Memory",
  "KLRG1" = "Term Diff",
  "CD38" = "Activation", "HLA_DR" = "Activation",
  "FCGR3A" = "NK-like Receptor Program", "NCAM1" = "NK-like Receptor Program",
  "TCR-AB" = "Gamma Delta T Identity", "CD3" = "Gamma Delta T Identity",
  "CD8A" = "Lineage"
)

p_cd8_rna <- build_annotation_heatmap(OPIS_CD8, "RNA", cd8_rna_markers, title = "CD8+ T Cells — mRNA (RNA)")
p_cd8_adt <- build_annotation_heatmap(OPIS_CD8, "ADT", cd8_adt_markers, title = "CD8+ T Cells — Protein (ADT)")

if (!is.null(p_cd8_adt) && !is.null(p_cd8_rna)) {
  png(file.path(out_heat, "CD8_ADT_RNA_Heatmap_Combined.png"),
      width = 24, height = 16, units = "in", res = 300, bg = "white")
  grid.newpage()
  grid.text("CD8+ T Cell Annotation Validation: Protein (ADT) & mRNA",
            x = 0.5, y = 0.975, gp = gpar(fontsize = 18, fontface = "bold"))
  pushViewport(viewport(layout = grid.layout(1, 2)))
  pushViewport(viewport(layout.pos.col = 1)); grid.draw(p_cd8_adt$gtable); popViewport()
  pushViewport(viewport(layout.pos.col = 2)); grid.draw(p_cd8_rna$gtable); popViewport()
  dev.off()
}
if (!is.null(p_cd8_adt)) {
  png(file.path(out_heat, "CD8_ADT_Heatmap.png"), width = 11, height = 10, units = "in", res = 300, bg = "white")
  grid.draw(p_cd8_adt$gtable); dev.off()
}
if (!is.null(p_cd8_rna)) {
  png(file.path(out_heat, "CD8_RNA_Heatmap.png"), width = 11, height = 14, units = "in", res = 300, bg = "white")
  grid.draw(p_cd8_rna$gtable); dev.off()
}
message("✓ CD8 heatmaps saved")

# ══════════════════════════════════════════════════════════════════════════════
# 3C) NK Heatmaps — UPDATED: exclude "NK-like CD8+" from NK subset
# ══════════════════════════════════════════════════════════════════════════════

nk_labels <- grep("NK|Natural",
                  levels(factor(OPIS_NKCD8$celltype_annotation)),
                  value = TRUE, ignore.case = TRUE)
# UPDATED: remove "NK-like CD8+" from NK labels
nk_labels <- setdiff(nk_labels, "NK-like CD8+")
message("NK labels used (after excluding NK-like CD8+): ", paste(nk_labels, collapse = ", "))

OPIS_NK <- subset(OPIS_NKCD8, celltype_annotation %in% nk_labels)

nk_rna_markers <- c(
  "NCAM1" = "NK Identity", "FCGR3A" = "NK Identity", "NKG7" = "NK Identity",
  "GNLY" = "NK Identity", "PRF1" = "NK Identity", "KLRD1" = "NK Identity",
  "FCER1G" = "NK Identity", "TYROBP" = "NK Identity",
  "KLRC1" = "CD56+ Less Mature", "XCL1" = "CD56+ Less Mature", "XCL2" = "CD56+ Less Mature"
)

nk_adt_markers <- c(
  "NCAM1" = "NK Identity", "FCGR3A" = "NK Identity"
)

p_nk_rna <- build_annotation_heatmap(OPIS_NK, "RNA", nk_rna_markers, title = "NK Cells — mRNA (RNA)")
p_nk_adt <- build_annotation_heatmap(OPIS_NK, "ADT", nk_adt_markers, title = "NK Cells — Protein (ADT)")

if (!is.null(p_nk_adt) && !is.null(p_nk_rna)) {
  png(file.path(out_heat, "NK_ADT_RNA_Heatmap_Combined.png"),
      width = 22, height = 12, units = "in", res = 300, bg = "white")
  grid.newpage()
  grid.text("NK Cell Annotation Validation: Protein (ADT) & mRNA",
            x = 0.5, y = 0.975, gp = gpar(fontsize = 18, fontface = "bold"))
  pushViewport(viewport(layout = grid.layout(1, 2)))
  pushViewport(viewport(layout.pos.col = 1)); grid.draw(p_nk_adt$gtable); popViewport()
  pushViewport(viewport(layout.pos.col = 2)); grid.draw(p_nk_rna$gtable); popViewport()
  dev.off()
}
if (!is.null(p_nk_rna)) {
  png(file.path(out_heat, "NK_RNA_Heatmap.png"), width = 11, height = 10, units = "in", res = 300, bg = "white")
  grid.draw(p_nk_rna$gtable); dev.off()
}
message("✓ NK heatmaps saved (NK-like CD8+ excluded)")

################################################################################
# 4) Feature Plots & Violin Plots — CD4 and NK/CD8 subsets
#    Split-by-OUD violin plots, non-split violins, and feature plots
#    All organized under Manuscript_Figures/
################################################################################

# ── Feature lists ─────────────────────────────────────────────────────────────
rna.features <- c(
  'ASCL2','BATF','BATF3','BCL6','BACH2','C1QBP','CCL2','CCL3','CCL4L2','CCL5',
  'CCND3','CD14','CD19','CD1C','CD200','CD27','CD3D','CD3E','CD36','CD4',
  'CD40','CD40LG','CD70','CD7','CD79A','CD8A','CD8B','CCR5','CLEC9A','CR2',
  'CTLA4','CTSW','CXCL8','CXCR3','CXCR5','EBI3','ENTPD1','FABP5','FCGR2B',
  'FCGR3A','FCRL5','FOXP3','GNLY','GP1BA','GP9','GATA3','GZMK','HAVCR2',
  'HIF1A','HIST1H4C','HLA-DPA1','HLA-DRA','HLA-DRB1','ICOS','IFI30','IFNG',
  'IGFBP2','IGFBP4','IGHA1','IGHA2','IGHG1','IGHG2','IGHG3','IGHG4','IGHM',
  'IKZF2','IL10','IL17A','IL18BP','IL18RAP','IL1B','IL21','IL2RA','IL2RB',
  'IRF4','IRF8','ITGAX','JCHAIN','KLRB1','KLRD1','KLRC1','LAG3','LDHA',
  'LGALS1','LTA','LTB','MAF','MAL','MALAT1','MIR155HG','MKI67','MT-ND1',
  'MT-ND5','MS4A1','NELL2','NCAM1','NKG7','NR4A1','PDCD1','PF4','PPBP',
  'PRDM1','PRF1','RORC','SELL','SERPINA1','SERPING1','SH2D1A','TCF4','TCF7',
  'TIGIT','TNF','TNFAIP2','TNFRSF18','TNFRSF4','TNFRSF9','TOX','TBX21',
  'TRBC1','TRDC','TRDV1','TRDV2','TRGC1','TRGC2','TRGV9','XBP1','XCL1',
  'XCL2','ZBTB16','ZEB2'
)

rna.features.2 <- c(
  'TRAC','TRBC2','IL7R','FOS','FOSB','JUN','JUNB','JUND','EGR1','EGR2',
  'EGR3','NR4A2','NR4A3','DUSP1','DUSP2','DUSP4','DUSP5','IER2','IER3',
  'TNFAIP3','HSP90AA1','HSPA1A','HSPA1B','TNFRSF4','TNFRSF9','TNFRSF18',
  'CCR7','LEF1','KLF2','KLF3','CCL4','CXCR4','ITGA4','ITGB1','IL12RB2',
  'CCR6','IL23R','TOX2','IL21','GZMB','CTSW','FGFBP2','TOP2A','HMGB2',
  'TYMS','PCNA','STMN1','TUBB','UBE2C','CENPF','S100A4'
)

all_rna_features <- unique(c(rna.features, rna.features.2))

# ── Generic function for subset-level feature/violin plots ────────────────────
generate_subset_plots <- function(obj, subset_name,
                                  rna_features, rna_features_2,
                                  out_feat_rna, out_feat_adt,
                                  out_vln_rna_oud, out_vln_rna,
                                  out_vln_adt_oud, out_vln_adt) {
  
  all_rna <- unique(c(rna_features, rna_features_2))
  all_rna <- intersect(all_rna, rownames(obj[["RNA"]]))
  prots   <- rownames(obj[["ADT"]])
  
  # ── ADT plots ──
  DefaultAssay(obj) <- "ADT"
  for (feat in prots) {
    if (feat %in% rownames(obj[["ADT"]])) {
      # Split-by-OUD violin
      vln_split <- VlnPlot2(
        obj, features = feat, cols = "default",
        split.by = "OUD_status", stat.method = "wilcox.test",
        show.mean = TRUE, mean_colors = c("red", "blue")
      ) + ggtitle(paste("ADT |", feat, "| Split by OUD_status |", subset_name))
      ggsave(file.path(out_vln_adt_oud, paste0(feat, "_ADT_VLN_splitOUD.png")),
             plot = vln_split, dpi = 500, width = 14, height = 8, bg = "white")
      
      # Non-split violin
      vln <- VlnPlot2(
        obj, features = feat, cols = "default",
        show.mean = TRUE, mean_colors = c("red", "blue")
      ) + ggtitle(paste("ADT |", feat, "|", subset_name))
      ggsave(file.path(out_vln_adt, paste0(feat, "_ADT_VLN.png")),
             plot = vln, dpi = 500, width = 14, height = 8, bg = "white")
      
      # Feature plot
      pal <- viridis(n = 10, option = "A")
      fp <- FeaturePlot_scCustom(
        obj, reduction = "wnn.umap", features = feat,
        colors_use = pal, order = TRUE
      )
      ggsave(file.path(out_feat_adt, paste0(feat, "_ADT_FeaturePlot.png")),
             plot = fp, dpi = 500, width = 8, bg = "white")
    }
  }
  
  # ── RNA plots ──
  DefaultAssay(obj) <- "RNA"
  for (feat in all_rna) {
    if (feat %in% rownames(obj[["RNA"]])) {
      # Split-by-OUD violin
      vln_split <- VlnPlot2(
        obj, features = feat, cols = "default",
        split.by = "OUD_status", stat.method = "wilcox.test",
        show.mean = TRUE, mean_colors = c("red", "blue")
      ) + ggtitle(paste("RNA |", feat, "| Split by OUD_status |", subset_name))
      ggsave(file.path(out_vln_rna_oud, paste0(feat, "_RNA_VLN_splitOUD.png")),
             plot = vln_split, dpi = 500, width = 14, height = 8, bg = "white")
      
      # Non-split violin
      vln <- VlnPlot2(
        obj, features = feat, cols = "default",
        show.mean = TRUE, mean_colors = c("red", "blue")
      ) + ggtitle(paste("RNA |", feat, "|", subset_name))
      ggsave(file.path(out_vln_rna, paste0(feat, "_RNA_VLN.png")),
             plot = vln, dpi = 500, width = 14, height = 8, bg = "white")
      
      # Feature plot
      pal <- viridis(n = 10, option = "A")
      fp <- FeaturePlot_scCustom(
        obj, reduction = "wnn.umap", features = feat,
        colors_use = pal, order = TRUE
      )
      ggsave(file.path(out_feat_rna, paste0(feat, "_RNA_FeaturePlot.png")),
             plot = fp, dpi = 500, width = 8, bg = "white")
    }
  }
}

message("=== Generating CD4 subset feature & violin plots ===")
generate_subset_plots(
  obj = OPIS_CD4, subset_name = "CD4",
  rna_features = rna.features, rna_features_2 = rna.features.2,
  out_feat_rna = out_feat_cd4_rna, out_feat_adt = out_feat_cd4_adt,
  out_vln_rna_oud = out_vln_cd4_rna_oud, out_vln_rna = out_vln_cd4_rna,
  out_vln_adt_oud = out_vln_cd4_adt_oud, out_vln_adt = out_vln_cd4_adt
)

message("=== Generating NK/CD8 subset feature & violin plots ===")
generate_subset_plots(
  obj = OPIS_NKCD8, subset_name = "NKCD8",
  rna_features = rna.features, rna_features_2 = rna.features.2,
  out_feat_rna = out_feat_nkcd8_rna, out_feat_adt = out_feat_nkcd8_adt,
  out_vln_rna_oud = out_vln_nkcd8_rna_oud, out_vln_rna = out_vln_nkcd8_rna,
  out_vln_adt_oud = out_vln_nkcd8_adt_oud, out_vln_adt = out_vln_nkcd8_adt
)

message("✓ Section 4: Feature plots & violin plots for subsets complete")

################################################################################
# 4B) DGE on CD4 subtypes (subset clusters)
################################################################################

OPIS_CD4_OUD <- subset(OPIS_CD4, subset = !is.na(OUD_status))
OPIS_CD4_OUD$OUD_status <- factor(OPIS_CD4_OUD$OUD_status)

dge_cd4_rna <- run_clusterwise_de(
  obj = OPIS_CD4_OUD, assay = "RNA", group_col = "OUD_status",
  ident1 = ident1_OUD, ident2 = ident2_OUD,
  cluster_col = "celltype_annotation", out_dir = out_dge_cd4_rna,
  latent = c("nCount_RNA"), title_prefix = "OPIS CD4 Subtypes (RNA)",
  file_stub = "OUD_Pos_vs_Neg_RNA"
)

dge_cd4_adt <- run_clusterwise_de(
  obj = OPIS_CD4_OUD, assay = "ADT", group_col = "OUD_status",
  ident1 = ident1_OUD, ident2 = ident2_OUD,
  cluster_col = "celltype_annotation", out_dir = out_dge_cd4_adt,
  latent = c("nCount_ADT"), title_prefix = "OPIS CD4 Subtypes (ADT)",
  file_stub = "OUD_Pos_vs_Neg_ADT"
)

# Volcano + top-30 heatmaps from CSVs (CD4 subtypes)
make_plots_from_csvs(out_dge_cd4_rna, "RNA", top_n = 30)
make_plots_from_csvs(out_dge_cd4_adt, "ADT", top_n = 30)

message("✓ DGE (CD4 subtypes, RNA + ADT) complete")

# ── 4C) DGE on NK/CD8 subtypes ───────────────────────────────────────────────
OPIS_NKCD8_OUD <- subset(OPIS_NKCD8, subset = !is.na(OUD_status))
OPIS_NKCD8_OUD$OUD_status <- factor(OPIS_NKCD8_OUD$OUD_status)

dge_nkcd8_rna <- run_clusterwise_de(
  obj = OPIS_NKCD8_OUD, assay = "RNA", group_col = "OUD_status",
  ident1 = ident1_OUD, ident2 = ident2_OUD,
  cluster_col = "celltype_annotation", out_dir = out_dge_nkcd8_rna,
  latent = c("nCount_RNA"), title_prefix = "OPIS NK/CD8 Subtypes (RNA)",
  file_stub = "OUD_Pos_vs_Neg_RNA"
)

dge_nkcd8_adt <- run_clusterwise_de(
  obj = OPIS_NKCD8_OUD, assay = "ADT", group_col = "OUD_status",
  ident1 = ident1_OUD, ident2 = ident2_OUD,
  cluster_col = "celltype_annotation", out_dir = out_dge_nkcd8_adt,
  latent = c("nCount_ADT"), title_prefix = "OPIS NK/CD8 Subtypes (ADT)",
  file_stub = "OUD_Pos_vs_Neg_ADT"
)

# Volcano + top-30 heatmaps from CSVs (NK/CD8 subtypes)
make_plots_from_csvs(out_dge_nkcd8_rna, "RNA", top_n = 30)
make_plots_from_csvs(out_dge_nkcd8_adt, "ADT", top_n = 30)

message("✓ DGE (NK/CD8 subtypes, RNA + ADT) complete")

################################################################################
# 5) EnrichR Pathway Analysis — Publication-quality plots
#    Run on subset-level DGE/DPE results (CD4 & NK/CD8)
################################################################################

# ── Enrichr databases ─────────────────────────────────────────────────────────
databases <- c(
  "TRRUST_Transcription_Factors_2019", "ChEA_2022", "TRANSFAC_and_JASPAR_PWMs",
  "KEGG_2021_Human", "WikiPathways_2024_Human", "GO_Biological_Process_2023",
  "MSigDB_Hallmark_2020", "Panther_2016", "Reactome_2022", "BioPlanet_2019"
)

tf_databases <- c("TRRUST_Transcription_Factors_2019", "ChEA_2022", "TRANSFAC_and_JASPAR_PWMs")
pathway_databases <- setdiff(databases, tf_databases)

# ── Publication-quality enrichment theme ──────────────────────────────────────
theme_enrichr_pub <- function() {
  theme_classic(base_size = 13) %+replace%
    theme(
      plot.title       = element_text(size = 15, face = "bold", hjust = 0.5, margin = margin(b = 10)),
      plot.subtitle    = element_text(size = 11, hjust = 0.5, color = "grey40", margin = margin(b = 8)),
      axis.text.y      = element_text(size = 10, color = "black", lineheight = 0.9),
      axis.text.x      = element_text(size = 11, color = "black"),
      axis.title.x     = element_text(size = 12, face = "bold", margin = margin(t = 8)),
      axis.title.y     = element_blank(),
      legend.title     = element_text(size = 11, face = "bold"),
      legend.text      = element_text(size = 10),
      legend.position  = "bottom",
      legend.box       = "horizontal",
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
      plot.margin      = margin(15, 15, 10, 10)
    )
}

# ── Publication-quality color palettes for databases ──────────────────────────
db_colors_tf <- c(
  "TRRUST_Transcription_Factors_2019" = "#7B1FA2",
  "ChEA_2022"                         = "#E65100",
  "TRANSFAC_and_JASPAR_PWMs"          = "#00838F"
)

db_colors_pathway <- c(
  "KEGG_2021_Human"              = "#1565C0",
  "WikiPathways_2024_Human"      = "#2E7D32",
  "GO_Biological_Process_2023"   = "#C62828",
  "MSigDB_Hallmark_2020"         = "#F57F17",
  "Panther_2016"                 = "#6A1B9A",
  "Reactome_2022"                = "#00695C",
  "BioPlanet_2019"               = "#AD1457"
)

# ── Improved enrichment runner ────────────────────────────────────────────────
run_enrichment_pub <- function(gene_df, label, base_output) {
  gene_list <- rownames(gene_df)
  if (length(gene_list) == 0) {
    message("  No genes found for ", label)
    return(NULL)
  }
  
  dir.create(file.path(base_output, "CSVs"),  recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(base_output, "Plots"), recursive = TRUE, showWarnings = FALSE)
  
  message("  Submitting ", length(gene_list), " features to Enrichr: ", label)
  enrichment <- enrichr(gene_list, databases)
  
  # Save Excel
  excel_out <- file.path(base_output, "CSVs", paste0(label, "_Enrichment.xlsx"))
  wb <- createWorkbook()
  for (db in names(enrichment)) {
    addWorksheet(wb, substr(db, 1, 31))
    writeData(wb, substr(db, 1, 31), enrichment[[db]])
  }
  saveWorkbook(wb, excel_out, overwrite = TRUE)
  
  # Collect results
  top_tf_list      <- list()
  top_pathway_list <- list()
  
  for (db_name in names(enrichment)) {
    db_results <- enrichment[[db_name]]
    if ("Combined Score" %in% colnames(db_results)) {
      db_results <- db_results %>% rename(Combined.Score = `Combined Score`)
    }
    if (!"Combined.Score" %in% colnames(db_results)) next
    
    sig_results <- db_results %>% filter(Adjusted.P.value < 0.05)
    if (nrow(sig_results) > 0) {
      top_terms <- sig_results %>%
        arrange(desc(Combined.Score)) %>%
        slice_head(n = 10) %>%
        mutate(Database = db_name)
      
      if (db_name %in% tf_databases) {
        top_tf_list[[db_name]] <- top_terms
      } else if (db_name %in% pathway_databases) {
        top_pathway_list[[db_name]] <- top_terms
      }
    }
  }
  
  # ── TF Plot (publication quality) ──
  tf_df <- bind_rows(top_tf_list)
  if ("Combined.Score" %in% colnames(tf_df) && nrow(tf_df) > 0) {
    tf_df <- tf_df %>%
      arrange(desc(Combined.Score)) %>%
      slice_head(n = 20) %>%
      mutate(
        Term = gsub("_.*$", "", Term),    # Clean long term names
        Term = stringr::str_wrap(Term, width = 40),
        neg_log10_pval = -log10(Adjusted.P.value + 1e-300)
      )
    
    p_tf <- ggplot(tf_df, aes(x = reorder(Term, Combined.Score), y = Combined.Score)) +
      geom_segment(aes(xend = reorder(Term, Combined.Score), y = 0, yend = Combined.Score),
                   color = "grey60", linewidth = 0.5) +
      geom_point(aes(color = Database, size = neg_log10_pval)) +
      scale_color_manual(values = db_colors_tf, name = "Database") +
      scale_size_continuous(range = c(2, 7), name = "-log10(Adj. P)") +
      coord_flip() +
      theme_enrichr_pub() +
      labs(
        title    = paste("Transcription Factor Enrichment"),
        subtitle = label,
        x = "", y = "Combined Score"
      )
    
    ggsave(file.path(base_output, "Plots", paste0(label, "_TF_Enrichment.png")),
           plot = p_tf, width = 12, height = 10, dpi = 400, bg = "white")
    ggsave(file.path(base_output, "Plots", paste0(label, "_TF_Enrichment.pdf")),
           plot = p_tf, width = 12, height = 10)
  }
  
  # ── Pathway Plot (publication quality) ──
  pathway_df <- bind_rows(top_pathway_list)
  if ("Combined.Score" %in% colnames(pathway_df) && nrow(pathway_df) > 0) {
    pathway_df <- pathway_df %>%
      arrange(desc(Combined.Score)) %>%
      slice_head(n = 20) %>%
      mutate(
        Term = gsub(" \\(.*\\)$", "", Term),
        Term = stringr::str_wrap(Term, width = 45),
        neg_log10_pval = -log10(Adjusted.P.value + 1e-300)
      )
    
    p_path <- ggplot(pathway_df, aes(x = reorder(Term, Combined.Score), y = Combined.Score)) +
      geom_segment(aes(xend = reorder(Term, Combined.Score), y = 0, yend = Combined.Score),
                   color = "grey60", linewidth = 0.5) +
      geom_point(aes(color = Database, size = neg_log10_pval)) +
      scale_color_manual(values = db_colors_pathway, name = "Database") +
      scale_size_continuous(range = c(2, 7), name = "-log10(Adj. P)") +
      coord_flip() +
      theme_enrichr_pub() +
      labs(
        title    = paste("Pathway Enrichment"),
        subtitle = label,
        x = "", y = "Combined Score"
      )
    
    ggsave(file.path(base_output, "Plots", paste0(label, "_Pathway_Enrichment.png")),
           plot = p_path, width = 13, height = 10, dpi = 400, bg = "white")
    ggsave(file.path(base_output, "Plots", paste0(label, "_Pathway_Enrichment.pdf")),
           plot = p_path, width = 13, height = 10)
  }
}

# ── Run enrichment on all subset-level DGE results ───────────────────────────

# ── Run enrichment from saved DGE/DPE CSVs (no in-memory objects needed) ──────
# dge_base_dir: e.g. .../DGE/CD4_Subtypes/RNA  (reads from CSVs/ subfolder)
# enrichr_out_dir: e.g. .../Pathway_Analysis/CD4_Subtypes/RNA
run_enrichr_from_csvs <- function(dge_base_dir, enrichr_out_dir, assay_label) {
  csv_dir   <- file.path(dge_base_dir, "CSVs")
  csv_files <- list.files(csv_dir, pattern = "\\.csv$", full.names = TRUE)
  
  if (length(csv_files) == 0) {
    message("  No DGE CSVs found in: ", csv_dir)
    return(invisible(NULL))
  }
  
  for (csv_path in csv_files) {
    de_df <- tryCatch(
      read.csv(csv_path, row.names = 1, check.names = FALSE),
      error = function(e) { message("  Read failed: ", csv_path); return(NULL) }
    )
    if (is.null(de_df)) next
    if (!"p_val_adj" %in% colnames(de_df)) next
    
    sig_genes <- de_df %>% filter(!is.na(p_val_adj) & p_val_adj < 0.05)
    if (nrow(sig_genes) == 0) {
      message("  No significant features in ", basename(csv_path))
      next
    }
    
    # Derive label from filename
    fname <- tools::file_path_sans_ext(basename(csv_path))
    cl_name <- gsub("_OUD_Pos_vs_Neg_(RNA|ADT)$", "", fname)
    label   <- paste0(cl_name, "_OUD_Pos_vs_Neg_", assay_label)
    
    message("  EnrichR for: ", label)
    run_enrichment_pub(sig_genes, label, enrichr_out_dir)
  }
}

message("=== Running EnrichR on CD4 subtypes (RNA) ===")
run_enrichr_from_csvs(out_dge_cd4_rna, out_path_cd4_rna, "RNA")

message("=== Running EnrichR on CD4 subtypes (ADT) ===")
run_enrichr_from_csvs(out_dge_cd4_adt, out_path_cd4_adt, "ADT")

message("=== Running EnrichR on NK/CD8 subtypes (RNA) ===")
run_enrichr_from_csvs(out_dge_nkcd8_rna, out_path_nkcd8_rna, "RNA")

message("=== Running EnrichR on NK/CD8 subtypes (ADT) ===")
run_enrichr_from_csvs(out_dge_nkcd8_adt, out_path_nkcd8_adt, "ADT")

message("✓ Section 5: EnrichR pathway analysis complete")

################################################################################
message("\n══════════════════════════════════════════════════════════════")
message("✓ ALL SECTIONS COMPLETE")
message("  Output root: ", out_base)
message("══════════════════════════════════════════════════════════════")
################################################################################

