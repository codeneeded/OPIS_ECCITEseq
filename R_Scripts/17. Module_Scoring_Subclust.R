# ============================================================================
# OPIS — Module scores + Cohen's d (OUD+ vs OUD-)
# CD4 (Fig 3) and CD8 (Fig 4) scored & plotted independently.
# Heatmap layout: clusters = columns, modules = rows, fill = Cohen's d.
# ============================================================================

# -------- Libraries ---------------------------------------------------------
library(Seurat)
library(SeuratObject)
library(qs2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(scales)
library(scCustomize)
library(SeuratExtend)
library(viridis)
library(openxlsx)
library(patchwork)

# ============================================================================
# CONFIG — EDIT THIS BLOCK
# ============================================================================

load.path     <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/saved_R_data"
analysis_dir  <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/Module_Scoring_Subclust"
fig_dir_cd4   <- file.path(analysis_dir, "Fig3_CD4_Modules")
fig_dir_cd8   <- file.path(analysis_dir, "Fig4_CD8_Modules")
dir.create(fig_dir_cd4, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir_cd8, recursive = TRUE, showWarnings = FALSE)

# Metadata columns
oud_col          <- "OUD_status"              # OUD+ / OUD-
oud_pos_label    <- "OUD+"
oud_neg_label    <- "OUD-"
cd4_cluster_col  <- "celltype_annotation"
cd8_cluster_col  <- "celltype_annotation"

# Restrict to specific clusters for each panel
# CD4 (note: "CD4+ Naïve" uses ï — copy exactly from levels())
cd4_clusters_to_keep <- c("CD4+ Naïve",
                          "CD4+ Tcm",
                          "CD4+ Treg",
                          "CD4+ KLRG1+ Tem")

# CD8 (note: "CD8+ Tem" — there's also a "CD8+ Tem (PRF1+)" cluster you didn't list;
# only the plain "CD8+ Tem" is included. Add the other if you want both pooled separately.)
cd8_clusters_to_keep <- c("CD8+ Memory",
                          "CD8+ Tem",
                          "CD8+ TEMRA",
                          "Innate-like T")

# X-axis order on the heatmap / bar dodge order — naive → memory → effector
cd4_cluster_order <- c("CD4+ Naïve",
                       "CD4+ Tcm",
                       "CD4+ Treg",
                       "CD4+ KLRG1+ Tem")

cd8_cluster_order <- c("CD8+ Memory",
                       "CD8+ Tem",
                       "CD8+ TEMRA",
                       "Innate-like T")

# Effect-size color limits (symmetric). NULL = auto from data.
cd_color_limit <- NULL    # e.g. 1.0

# Per-cluster fill colors for the faceted bar plots
cd4_cluster_palette <- c("CD4+ Naïve"      = "#3B7FB8",   # steel blue
                         "CD4+ Tcm"        = "#4FA864",   # green
                         "CD4+ Treg"       = "#8E6BB0",   # purple
                         "CD4+ KLRG1+ Tem" = "#D95F02")   # burnt orange

cd8_cluster_palette <- c("CD8+ Memory"   = "#4FA864",     # green
                         "CD8+ Tem"      = "#D95F02",     # burnt orange
                         "CD8+ TEMRA"    = "#C4463A",     # brick red
                         "Innate-like T" = "#7570B3")     # violet

# ============================================================================
# MODULE GENE LISTS (Fig 3 — CD4)
# ============================================================================

cd4_modules <- list(
  Naive_Persistence     = c("CCR7","SELL","KLF2","LEF1","TCF7","SATB1","MAL",
                            "LTB","IL7R","BCL2"),
  IL2_STAT5_Cytokine    = c("IL2RA","STAT5A","STAT5B","SOCS1","SOCS3","BCL2",
                            "PIM1","CISH","JAK1","JAK3"),
  NFkB_Inflammatory     = c("NFKBIA","TNFAIP3","JUN","JUNB","FOS","RELB",
                            "BCL3","IL32","IER2"),
  Cytotoxic_Terminal    = c("CCL5","PRF1","GZMB","GZMK","GNLY","NKG7","KLRG1",
                            "ITGB1","CX3CR1","FGFBP2"),   # CXCR3R1 -> CX3CR1
  Interferon_Response   = c("IFITM1","IFITM2","IFI6","ISG15","MX1","STAT1",
                            "IRF7","OAS1","OASL","HLA-B","B2M"),
  Glycolysis            = c("SLC2A1","HK1","HK2","PFKP","PFKFB3","ALDOA",
                            "GAPDH","PGK1","ENO1","PKM","LDHA","PDK1"),
  OXPHOS                = c("NDUFA1","NDUFA2","NDUFB3","NDUFB8","SDHA","SDHB",
                            "UQCRC1","UQCRC2","COX5A","COX6C",
                            "ATP5F1A","ATP5F1B"),          # dropped lone "NDUFB"
  Mitochondrial_Fitness = c("PPARGC1A","TFAM","NRF1","ESRRA","CPT1A","SIRT1",
                            "PRKAA1","PRKAA2","BCL2"),
  Autophagy_Mitophagy   = c("MAP1LC3B","GABARAPL1","ATG5","ATG7","BECN1",
                            "SQSTM1","OPTN","PINK1","PRKN","BNIP3","BNIP3L")
)

cd4_module_labels <- c(
  Naive_Persistence     = "Naïve / persistence",
  IL2_STAT5_Cytokine    = "IL-2 / STAT5 / cytokine",
  NFkB_Inflammatory     = "NF-κB / inflammatory",
  Cytotoxic_Terminal    = "Cytotoxic / terminal",
  Interferon_Response   = "Interferon response",
  Glycolysis            = "Glycolysis",
  OXPHOS                = "OXPHOS",
  Mitochondrial_Fitness = "Mitochondrial fitness",
  Autophagy_Mitophagy   = "Autophagy / mitophagy"
)

cd4_module_order_phenotypic <- c(
  "Naïve / persistence",
  "IL-2 / STAT5 / cytokine",
  "NF-κB / inflammatory",
  "Cytotoxic / terminal",
  "Interferon response"
)

cd4_module_order_metabolic <- c(
  "Glycolysis",
  "OXPHOS",
  "Mitochondrial fitness",
  "Autophagy / mitophagy"
)

# ============================================================================
# MODULE GENE LISTS (Fig 4 — CD8)
# ============================================================================

cd8_modules <- list(
  Cytotoxic_Effector    = c("NKG7","PRF1","GZMB","GZMH","GZMA","GNLY","CCL5",
                            "CST7","FGFBP2"),
  TEMRA_Terminal        = c("KLRG1","CX3CR1","FGFBP2","PRF1","GZMB","TBX21",
                            "ZEB2","S1PR5"),
  Exhaustion            = c("PDCD1","LAG3","TIGIT","HAVCR2","TOX","CTLA4",
                            "ENTPD1","CXCL13"),
  Memory_StemLike       = c("TCF7","LEF1","CCR7","SELL","IL7R","BCL2","LTB"),
  Innate_NKlike         = c("KLRD1","KLRC1","KLRC2","KLRB1","NKG7","GNLY",
                            "TYROBP","FCGR3A","ZNF683","CXCR6"),
  Interferon_Response   = c("ISG15","IFITM1","IFITM2","IFI6","MX1","OAS1",
                            "OASL","STAT1","IRF7"),
  Glycolysis            = c("SLC2A1","HK1","HK2","PFKP","PFKFB3","ALDOA",
                            "GAPDH","PGK1","ENO1","PKM","LDHA","PDK1"),
  OXPHOS                = c("NDUFA1","NDUFA2","NDUFB3","NDUFB8","SDHA","SDHB",
                            "UQCRC1","UQCRC2","COX5A","COX6C",
                            "ATP5F1A","ATP5F1B"),
  Mitochondrial_Fitness = c("PPARGC1A","TFAM","NRF1","ESRRA","CPT1A","SIRT1",
                            "PRKAA1","PRKAA2","BCL2"),
  Autophagy_Mitophagy   = c("MAP1LC3B","GABARAPL1","ATG5","ATG7","BECN1",
                            "SQSTM1","OPTN","PINK1","PRKN","BNIP3","BNIP3L")
)

cd8_module_labels <- c(
  Cytotoxic_Effector    = "Cytotoxic effector",
  TEMRA_Terminal        = "TEMRA / terminal",
  Exhaustion            = "Exhaustion",
  Memory_StemLike       = "Memory / stem-like",
  Innate_NKlike         = "Innate-like / NK-like",
  Interferon_Response   = "Interferon response",
  Glycolysis            = "Glycolysis",
  OXPHOS                = "OXPHOS",
  Mitochondrial_Fitness = "Mitochondrial fitness",
  Autophagy_Mitophagy   = "Autophagy / mitophagy"
)

cd8_module_order_phenotypic <- c(
  "Memory / stem-like",
  "Cytotoxic effector",
  "TEMRA / terminal",
  "Exhaustion",
  "Innate-like / NK-like",
  "Interferon response"
)

cd8_module_order_metabolic <- c(
  "Glycolysis",
  "OXPHOS",
  "Mitochondrial fitness",
  "Autophagy / mitophagy"
)

# ============================================================================
# HELPERS
# ============================================================================

# --- Score modules and return a tidy per-cell data.frame
score_modules <- function(seu, modules, cluster_col, oud_col,
                          oud_pos, oud_neg, min_genes = 3) {
  
  set.seed(42)
  DefaultAssay(seu) <- "RNA"
  
  # Keep only genes present in the assay
  modules_filt <- lapply(modules, function(g) intersect(g, rownames(seu)))
  too_small    <- vapply(modules_filt, length, integer(1)) < min_genes
  
  if (any(too_small)) {
    message("  -> Dropping modules with <", min_genes, " detected genes: ",
            paste(names(modules_filt)[too_small], collapse = ", "))
    modules_filt <- modules_filt[!too_small]
  }
  if (length(modules_filt) == 0)
    stop("No usable modules remain after gene filtering.")
  
  # Report module sizes
  msg_sizes <- paste0(names(modules_filt), " (", lengths(modules_filt), ")",
                      collapse = ", ")
  message("  -> Scoring: ", msg_sizes)
  
  seu <- AddModuleScore(
    seu,
    features = modules_filt,
    name     = "MODULE_",
    nbin     = 24,
    ctrl     = 100,
    seed     = 42
  )
  
  score_cols <- paste0("MODULE_", seq_along(modules_filt))
  
  df <- data.frame(
    cell    = colnames(seu),
    Cluster = as.character(seu@meta.data[[cluster_col]]),
    OUD     = as.character(seu@meta.data[[oud_col]]),
    seu@meta.data[, score_cols, drop = FALSE],
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  colnames(df)[match(score_cols, colnames(df))] <- names(modules_filt)
  
  # Restrict to the two OUD groups of interest
  df <- df %>% filter(OUD %in% c(oud_pos, oud_neg))
  
  list(seurat = seu, scores = df, modules_used = names(modules_filt))
}

# --- Cohen's d per (module x cluster), OUD+ vs OUD-
compute_cohens_d <- function(scores_df, modules_used,
                             oud_pos, oud_neg, min_n = 3) {
  
  results <- list()
  clusters <- sort(unique(scores_df$Cluster))
  clusters <- clusters[!is.na(clusters) & clusters != ""]
  
  for (mod in modules_used) {
    for (cl in clusters) {
      sub <- scores_df %>% filter(Cluster == cl)
      x <- sub %>% filter(OUD == oud_pos) %>% pull(!!sym(mod))
      y <- sub %>% filter(OUD == oud_neg) %>% pull(!!sym(mod))
      
      if (length(x) >= min_n && length(y) >= min_n) {
        nx <- length(x); ny <- length(y)
        sx <- sd(x);     sy <- sd(y)
        pooled_sd <- sqrt(((nx - 1) * sx^2 + (ny - 1) * sy^2) /
                            (nx + ny - 2))
        cd <- if (!is.na(pooled_sd) && pooled_sd > 0)
          (mean(x) - mean(y)) / pooled_sd else NA_real_
        pv <- tryCatch(wilcox.test(x, y)$p.value,
                       error = function(e) NA_real_)
        results[[paste(mod, cl, sep = "__")]] <- data.frame(
          Module   = mod,
          Cluster  = cl,
          n_pos    = nx,
          n_neg    = ny,
          mean_pos = mean(x),
          mean_neg = mean(y),
          cohens_d = cd,
          p_value  = pv,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  out <- do.call(rbind, results)
  if (is.null(out) || nrow(out) == 0) return(NULL)
  
  out %>%
    mutate(p_adj = p.adjust(p_value, method = "BH"),
           star  = ifelse(is.na(p_adj), "",
                          ifelse(p_adj < 0.001, "***",
                                 ifelse(p_adj < 0.01,  "**",
                                        ifelse(p_adj < 0.05,  "*", "")))))
}

# --- Heatmap: clusters = columns, modules = rows
plot_cohens_d_heatmap <- function(cd_df, module_label_map, module_order,
                                  cluster_order = NULL,
                                  color_limit   = NULL,
                                  title) {
  
  cd_df$Module_Label <- module_label_map[cd_df$Module]
  cd_df$Module_Label <- factor(cd_df$Module_Label, levels = rev(module_order))
  
  if (!is.null(cluster_order)) {
    cd_df$Cluster <- factor(cd_df$Cluster, levels = cluster_order)
  } else {
    cd_df$Cluster <- factor(cd_df$Cluster, levels = sort(unique(cd_df$Cluster)))
  }
  
  cd_df <- cd_df %>% filter(!is.na(Cluster) & !is.na(Module_Label))
  
  lim <- if (is.null(color_limit))
    max(abs(cd_df$cohens_d), na.rm = TRUE) else color_limit
  lim <- max(lim, 0.3)  # avoid washed-out plots when effects are tiny
  
  ggplot(cd_df, aes(x = Cluster, y = Module_Label, fill = cohens_d)) +
    geom_tile(color = "white", linewidth = 1.2) +
    geom_text(aes(label = sprintf("%.2f", cohens_d)),
              color = "black", size = 3.6, vjust = -0.4) +
    geom_text(aes(label = star),
              color = "black", size = 5.5, fontface = "bold", vjust = 1.2) +
    scale_fill_gradient2(low = "#3A6FB0", mid = "white", high = "#C4463A",
                         midpoint = 0,
                         limits   = c(-lim, lim),
                         oob      = scales::squish,
                         name     = "Cohen's d\n(OUD+ vs OUD-)") +
    labs(x = NULL, y = NULL, title = title) +
    theme_cowplot(font_size = 16) +
    theme(axis.text.x     = element_text(size = 12, face = "bold",
                                         angle = 30, hjust = 1),
          axis.text.y     = element_text(size = 13, face = "bold"),
          plot.title      = element_text(size = 18, face = "bold", hjust = 0.5),
          plot.background = element_rect(fill = "white", color = NA),
          legend.position = "right")
}

# --- Faceted bar plot: one panel per cluster, modules on y, Cohen's d on x
plot_cohens_d_bars <- function(cd_df, module_label_map, module_order,
                               cluster_order = NULL,
                               cluster_palette = NULL, title) {
  
  cd_df$Module_Label <- module_label_map[cd_df$Module]
  cd_df$Module_Label <- factor(cd_df$Module_Label, levels = rev(module_order))
  
  if (!is.null(cluster_order)) {
    cd_df$Cluster <- factor(cd_df$Cluster, levels = cluster_order)
  } else {
    cd_df$Cluster <- factor(cd_df$Cluster, levels = sort(unique(cd_df$Cluster)))
  }
  
  cd_df <- cd_df %>% filter(!is.na(Cluster) & !is.na(Module_Label))
  
  # Symmetric x-axis with padding so stars don't get clipped
  x_max <- max(abs(cd_df$cohens_d), na.rm = TRUE)
  x_max <- max(x_max, 0.3)
  star_offset <- x_max * 0.06
  x_pad       <- x_max * 0.22
  xlim_vals   <- c(-x_max - x_pad, x_max + x_pad)
  
  cd_df$star_x     <- ifelse(cd_df$cohens_d > 0,
                             cd_df$cohens_d + star_offset,
                             cd_df$cohens_d - star_offset)
  cd_df$star_hjust <- ifelse(cd_df$cohens_d > 0, 0, 1)
  
  p <- ggplot(cd_df, aes(x = cohens_d, y = Module_Label, fill = Cluster)) +
    geom_col(width = 0.78, alpha = 0.95,
             color = "grey15", linewidth = 0.35) +
    geom_vline(xintercept = 0, linewidth = 0.5, color = "grey30") +
    geom_text(aes(x = star_x, label = star, hjust = star_hjust),
              size = 6, color = "black", vjust = 0.55, fontface = "bold") +
    facet_wrap(~ Cluster, nrow = 1) +
    scale_x_continuous(limits = xlim_vals, expand = expansion(mult = 0.02)) +
    coord_cartesian(clip = "off") +
    labs(x = "Cohen's d  (OUD+ vs OUD-)", y = NULL, title = title) +
    theme_cowplot(font_size = 16) +
    theme(strip.text       = element_text(size = 13, face = "bold"),
          strip.background = element_rect(fill = "grey92", color = NA),
          axis.text.y      = element_text(size = 12, face = "bold"),
          axis.text.x      = element_text(size = 11),
          axis.title.x     = element_text(size = 14, face = "bold",
                                          margin = margin(t = 8)),
          panel.spacing.x  = unit(1.4, "lines"),
          panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
          plot.title       = element_text(size = 18, face = "bold", hjust = 0.5),
          plot.margin      = margin(14, 24, 12, 14),
          plot.background  = element_rect(fill = "white", color = NA),
          legend.position  = "none")   # cluster shown in facet strip
  
  if (!is.null(cluster_palette))
    p <- p + scale_fill_manual(values = cluster_palette, drop = FALSE)
  
  p
}

# ============================================================================
# LOAD OBJECTS
# ============================================================================

cat("\nLoading Seurat objects...\n")
OPIS_CD4   <- qs_read(file.path(load.path, "OPIS_CD4_Annotated.qs2"))
OPIS_NKCD8 <- qs_read(file.path(load.path, "OPIS_NKCD8_Annotated.qs2"))

# Sanity checks
for (req in list(list(OPIS_CD4,   c(oud_col, cd4_cluster_col), "OPIS_CD4"),
                 list(OPIS_NKCD8, c(oud_col, cd8_cluster_col), "OPIS_NKCD8"))) {
  obj <- req[[1]]; cols <- req[[2]]; nm <- req[[3]]
  miss <- setdiff(cols, colnames(obj@meta.data))
  if (length(miss))
    stop(nm, " is missing required metadata columns: ",
         paste(miss, collapse = ", "))
}

# Show OUD level counts so you can confirm labels match the config
cat("\nOUD levels present (CD4):\n");   print(table(OPIS_CD4@meta.data[[oud_col]],   useNA = "ifany"))
cat("\nOUD levels present (NKCD8):\n"); print(table(OPIS_NKCD8@meta.data[[oud_col]], useNA = "ifany"))

# ============================================================================
# CD4 — Fig 3
# ============================================================================

cat("\n=== CD4 module scoring & Cohen's d ===\n")
if (!is.null(cd4_clusters_to_keep)) {
  miss <- setdiff(cd4_clusters_to_keep,
                  unique(OPIS_CD4@meta.data[[cd4_cluster_col]]))
  if (length(miss))
    warning("CD4 clusters not found in object: ", paste(miss, collapse = " | "))
  keep <- colnames(OPIS_CD4)[
    OPIS_CD4@meta.data[[cd4_cluster_col]] %in% cd4_clusters_to_keep
  ]
  OPIS_CD4 <- subset(OPIS_CD4, cells = keep)
  cat("  -> Subset CD4 to ", length(keep), " cells across ",
      length(cd4_clusters_to_keep), " cluster(s).\n", sep = "")
}

cd4_out <- score_modules(
  seu         = OPIS_CD4,
  modules     = cd4_modules,
  cluster_col = cd4_cluster_col,
  oud_col     = oud_col,
  oud_pos     = oud_pos_label,
  oud_neg     = oud_neg_label
)
write.csv(cd4_out$scores,
          file.path(fig_dir_cd4, "CD4_ModuleScores_per_cell.csv"),
          row.names = FALSE)

cd4_cd <- compute_cohens_d(
  scores_df    = cd4_out$scores,
  modules_used = cd4_out$modules_used,
  oud_pos      = oud_pos_label,
  oud_neg      = oud_neg_label
)
write.csv(cd4_cd,
          file.path(fig_dir_cd4, "CD4_CohenD_OUDpos_vs_OUDneg.csv"),
          row.names = FALSE)

p_cd4_hm_phen <- plot_cohens_d_heatmap(
  cd_df            = cd4_cd,
  module_label_map = cd4_module_labels,
  module_order     = cd4_module_order_phenotypic,
  cluster_order    = cd4_cluster_order,
  color_limit      = cd_color_limit,
  title            = "CD4: phenotypic modules — Cohen's d (OUD+ vs OUD-)"
)
ggsave(file.path(fig_dir_cd4, "Fig3_CD4_CohenD_Heatmap_Phenotypic.png"),
       p_cd4_hm_phen, width = 10, height = 6, dpi = 300, bg = "white")

p_cd4_hm_met <- plot_cohens_d_heatmap(
  cd_df            = cd4_cd,
  module_label_map = cd4_module_labels,
  module_order     = cd4_module_order_metabolic,
  cluster_order    = cd4_cluster_order,
  color_limit      = cd_color_limit,
  title            = "CD4: metabolic modules — Cohen's d (OUD+ vs OUD-)"
)
ggsave(file.path(fig_dir_cd4, "Fig3_CD4_CohenD_Heatmap_Metabolic.png"),
       p_cd4_hm_met, width = 10, height = 5, dpi = 300, bg = "white")

p_cd4_bar_phen <- plot_cohens_d_bars(
  cd_df            = cd4_cd,
  module_label_map = cd4_module_labels,
  module_order     = cd4_module_order_phenotypic,
  cluster_order    = cd4_cluster_order,
  cluster_palette  = cd4_cluster_palette,
  title            = "CD4: phenotypic modules — Cohen's d (OUD+ vs OUD-)"
)
ggsave(file.path(fig_dir_cd4, "Fig3_CD4_CohenD_Bars_Phenotypic.png"),
       p_cd4_bar_phen, width = 18, height = 6.5, dpi = 300, bg = "white")

p_cd4_bar_met <- plot_cohens_d_bars(
  cd_df            = cd4_cd,
  module_label_map = cd4_module_labels,
  module_order     = cd4_module_order_metabolic,
  cluster_order    = cd4_cluster_order,
  cluster_palette  = cd4_cluster_palette,
  title            = "CD4: metabolic modules — Cohen's d (OUD+ vs OUD-)"
)
ggsave(file.path(fig_dir_cd4, "Fig3_CD4_CohenD_Bars_Metabolic.png"),
       p_cd4_bar_met, width = 18, height = 5.5, dpi = 300, bg = "white")

# ============================================================================
# CD8 — Fig 4
# ============================================================================

cat("\n=== CD8 module scoring & Cohen's d ===\n")
OPIS_CD8 <- OPIS_NKCD8
if (!is.null(cd8_clusters_to_keep)) {
  miss <- setdiff(cd8_clusters_to_keep,
                  unique(OPIS_CD8@meta.data[[cd8_cluster_col]]))
  if (length(miss))
    warning("CD8 clusters not found in object: ", paste(miss, collapse = " | "))
  keep <- colnames(OPIS_CD8)[
    OPIS_CD8@meta.data[[cd8_cluster_col]] %in% cd8_clusters_to_keep
  ]
  OPIS_CD8 <- subset(OPIS_CD8, cells = keep)
  cat("  -> Subset CD8 to ", length(keep), " cells across ",
      length(cd8_clusters_to_keep), " cluster(s).\n", sep = "")
}

cd8_out <- score_modules(
  seu         = OPIS_CD8,
  modules     = cd8_modules,
  cluster_col = cd8_cluster_col,
  oud_col     = oud_col,
  oud_pos     = oud_pos_label,
  oud_neg     = oud_neg_label
)
write.csv(cd8_out$scores,
          file.path(fig_dir_cd8, "CD8_ModuleScores_per_cell.csv"),
          row.names = FALSE)

cd8_cd <- compute_cohens_d(
  scores_df    = cd8_out$scores,
  modules_used = cd8_out$modules_used,
  oud_pos      = oud_pos_label,
  oud_neg      = oud_neg_label
)
write.csv(cd8_cd,
          file.path(fig_dir_cd8, "CD8_CohenD_OUDpos_vs_OUDneg.csv"),
          row.names = FALSE)

p_cd8_hm_phen <- plot_cohens_d_heatmap(
  cd_df            = cd8_cd,
  module_label_map = cd8_module_labels,
  module_order     = cd8_module_order_phenotypic,
  cluster_order    = cd8_cluster_order,
  color_limit      = cd_color_limit,
  title            = "CD8: phenotypic modules — Cohen's d (OUD+ vs OUD-)"
)
ggsave(file.path(fig_dir_cd8, "Fig4_CD8_CohenD_Heatmap_Phenotypic.png"),
       p_cd8_hm_phen, width = 10, height = 7, dpi = 300, bg = "white")

p_cd8_hm_met <- plot_cohens_d_heatmap(
  cd_df            = cd8_cd,
  module_label_map = cd8_module_labels,
  module_order     = cd8_module_order_metabolic,
  cluster_order    = cd8_cluster_order,
  color_limit      = cd_color_limit,
  title            = "CD8: metabolic modules — Cohen's d (OUD+ vs OUD-)"
)
ggsave(file.path(fig_dir_cd8, "Fig4_CD8_CohenD_Heatmap_Metabolic.png"),
       p_cd8_hm_met, width = 10, height = 5, dpi = 300, bg = "white")

p_cd8_bar_phen <- plot_cohens_d_bars(
  cd_df            = cd8_cd,
  module_label_map = cd8_module_labels,
  module_order     = cd8_module_order_phenotypic,
  cluster_order    = cd8_cluster_order,
  cluster_palette  = cd8_cluster_palette,
  title            = "CD8: phenotypic modules — Cohen's d (OUD+ vs OUD-)"
)
ggsave(file.path(fig_dir_cd8, "Fig4_CD8_CohenD_Bars_Phenotypic.png"),
       p_cd8_bar_phen, width = 18, height = 7, dpi = 300, bg = "white")

p_cd8_bar_met <- plot_cohens_d_bars(
  cd_df            = cd8_cd,
  module_label_map = cd8_module_labels,
  module_order     = cd8_module_order_metabolic,
  cluster_order    = cd8_cluster_order,
  cluster_palette  = cd8_cluster_palette,
  title            = "CD8: metabolic modules — Cohen's d (OUD+ vs OUD-)"
)
ggsave(file.path(fig_dir_cd8, "Fig4_CD8_CohenD_Bars_Metabolic.png"),
       p_cd8_bar_met, width = 18, height = 5.5, dpi = 300, bg = "white")

cat("\nDone.\n  CD4 outputs: ", fig_dir_cd4,
    "\n  CD8 outputs: ", fig_dir_cd8, "\n", sep = "")