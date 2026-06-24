# =============================================================================
# OPIS ECCITEseq — B_Cell_Updates
# One drop-in script that patches the B-cell outputs WITHOUT rerunning the
# expensive steps (no re-subclustering, no MAST). It only:
#
#   (1) Green-palette heatmaps (CD8-manuscript style) added ALONGSIDE the
#       existing magma heatmaps — for Subcluster_ID (pre-annotation),
#       Bcell_Annotation, and Bcell_Final.
#   (2) % and raw-count frequency tables/plots for the final merged B clusters
#       -> final_merged/frequencies/
#   (3) DGE up + down: script 20's full table (DGE_MAST_full_allClusters.csv)
#       was written with only.pos = FALSE, so it already holds both directions.
#       Its significant export was up-only only because of avg_log2FC > 0.25.
#       Re-threshold the SAME full table with abs(avg_log2FC) > 0.25 into one
#       file. No MAST rerun.
#   (4) EnrichR pathway analysis on the B-cell merged-cluster DGE (all
#       significant genes per cluster, same approach as the broad EnrichR).
#   (5) EnrichR pathway analysis on the broad-cluster DGE CSVs in
#       Manuscript_Figures/DGE/Broad_Clusters/RNA/CSVs (B, CD4, ...).
#
# Reads : OPIS_BCELL_PreAnnotation.qs2, OPIS_BCELL_Annotated.qs2,
#         OPIS_BCELL_FinalMerged.qs2, DGE_MAST_full_allClusters.csv
# Writes: new files only; nothing existing is overwritten/deleted.
# =============================================================================

# ---- Toggle sections (run any subset) ---------------------------------------
RUN_GREEN_HEATMAPS <- TRUE
RUN_FREQUENCIES    <- TRUE
RUN_DGE_UPDOWN     <- TRUE
RUN_ENRICHR_BCELL  <- TRUE   # needs internet + enrichR + openxlsx
RUN_ENRICHR_BROAD  <- TRUE   # needs internet + enrichR + openxlsx

# ---- Libraries ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratExtend)
  library(qs2)
  library(tidyverse)
  library(pheatmap)
})

# ---- Paths (match the existing pipeline exactly) -----------------------------
load.path       <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/saved_R_data"
subclust.root   <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/subclustering/pre_annotation"
bcell.save.path <- file.path(subclust.root, "Bcells")
final.dir       <- file.path(bcell.save.path, "final_merged")
dge.dir         <- file.path(final.dir, "DGE_MAST")

# Broad-cluster DGE CSVs (request 5)
broad.csv.dir   <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/Manuscript_Figures/DGE/Broad_Clusters/RNA/CSVs"
broad.enr.out   <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/Manuscript_Figures/DGE/Broad_Clusters/RNA/Pathway_Analysis_EnrichR"

# =============================================================================
# SHARED HELPERS
# =============================================================================
.norm <- function(x) toupper(gsub("[-_. ]", "", x))

# min-max scale each row 0-1 (DSB ADT can be negative; min->0 is fine)
scale_01 <- function(m) {
  m <- as.matrix(m)
  out <- t(apply(m, 1, function(x) {
    rng <- range(x, na.rm = TRUE)
    if (!is.finite(diff(rng)) || diff(rng) == 0) return(rep(0, length(x)))
    (x - rng[1]) / (rng[2] - rng[1])
  }))
  dimnames(out) <- dimnames(m); out
}
get_gaps <- function(annot_df) {
  g <- as.character(annot_df$Group); which(g[-length(g)] != g[-1])
}
cluster_within_groups <- function(mat, annot_df) {
  ro <- character(0)
  for (grp in levels(droplevels(annot_df$Group))) {
    rows <- rownames(annot_df)[annot_df$Group == grp]
    if (length(rows) <= 1) { ro <- c(ro, rows); next }
    hc <- hclust(dist(mat[rows, , drop = FALSE]), method = "ward.D2")
    ro <- c(ro, rows[hc$order])
  }
  ro
}
resolve_feature <- function(label, aliases, available) {
  cand <- c(label, if (nzchar(aliases)) strsplit(aliases, ",")[[1]] else character(0))
  cand <- trimws(cand)
  for (c0 in cand) if (c0 %in% available) return(c0)
  av_norm <- .norm(available)
  for (c0 in cand) { hit <- available[av_norm == .norm(c0)]; if (length(hit)) return(hit[1]) }
  NA_character_
}
prepare_panel <- function(markers, cat_levels, available, modality) {
  res <- markers %>%
    mutate(feature = mapply(resolve_feature, label, aliases,
                            MoreArgs = list(available = available)))
  miss <- res %>% filter(is.na(feature))
  if (nrow(miss)) message("[", modality, "] not found, skipped: ",
                          paste(miss$label, collapse = ", "))
  res <- res %>% filter(!is.na(feature)) %>%
    mutate(category = factor(category, levels = cat_levels)) %>% arrange(category)
  res[!duplicated(res$feature), , drop = FALSE]
}

# ---- Green palette (CD8 manuscript) + per-category row-annotation colours ----
heatmap_colors_green <- colorRampPalette(
  c("#F7FCF5", "#C7E9C0", "#74C476", "#31A354", "#006D2C"))(100)

bcell_cat_palette <- c(
  "Lineage"               = "#78909C",
  "Naive/Memory"          = "#74C2E1",
  "Naive B (resting)"     = "#74C2E1",
  "Transitional B"        = "#4FC3F7",
  "IEG-High"              = "#FFD54F",
  "Intermediate HLA-G+"   = "#AB47BC",
  "Atypical B"            = "#EF5350",
  "Atypical"              = "#EF5350",
  "Memory B"              = "#66BB6A",
  "Activation"            = "#FFA726",
  "Plasmablast exclusion" = "#A1887F"
)

# ---- Curated marker panels (identical to script 20) --------------------------
adt_markers <- tibble::tribble(
  ~category,      ~label,   ~aliases,
  "Lineage",      "CD19",   "",
  "Lineage",      "CD20",   "MS4A1",
  "Lineage",      "CD22",   "",
  "Lineage",      "CD79B",  "",
  "Lineage",      "HLA-DR", "HLA.DR,HLADR,HLA-DRA",
  "Naive/Memory", "IgD",    "IGHD",
  "Naive/Memory", "IgM",    "IGHM",
  "Naive/Memory", "CD27",   "",
  "Naive/Memory", "CD38",   "",
  "Naive/Memory", "CD24",   "",
  "Naive/Memory", "CD86",   "",
  "Atypical",     "CXCR5",  "",
  "Atypical",     "LILRB1", "CD85J",
  "Activation",   "CD40",   "",
  "Activation",   "CD69",   ""
)
adt_cat_levels <- c("Lineage", "Naive/Memory", "Atypical", "Activation")

rna_markers <- tibble::tribble(
  ~category,               ~label,      ~aliases,
  "Lineage",               "MS4A1",     "",
  "Lineage",               "CD79A",     "",
  "Lineage",               "CD79B",     "",
  "Lineage",               "CD74",      "",
  "Lineage",               "BANK1",     "",
  "Lineage",               "CD19",      "",
  "Lineage",               "BLK",       "",
  "Naive B (resting)",     "IGHD",      "",
  "Naive B (resting)",     "IGHM",      "",
  "Naive B (resting)",     "TCL1A",     "",
  "Naive B (resting)",     "IL4R",      "",
  "Naive B (resting)",     "FCER2",     "",
  "Naive B (resting)",     "CCR7",      "",
  "Naive B (resting)",     "CXCR4",     "",
  "Naive B (resting)",     "BACH2",     "",
  "Transitional B",        "CD38",      "",
  "Transitional B",        "CD10",      "MME",
  "Transitional B",        "SOX4",      "",
  "Transitional B",        "VPREB3",    "",
  "Transitional B",        "IGLL5",     "",
  "Transitional B",        "TCL1A",     "",
  "Transitional B",        "CD24",      "",
  "Transitional B",        "IL7R",      "",
  "Transitional B",        "CD21",      "CR2",
  "IEG-High",              "FOS",       "",
  "IEG-High",              "JUN",       "",
  "IEG-High",              "JUNB",      "",
  "IEG-High",              "DUSP1",     "",
  "IEG-High",              "EGR1",      "",
  "IEG-High",              "IER2",      "",
  "IEG-High",              "NR4A1",     "",
  "Intermediate HLA-G+",   "HLA-G",     "HLA.G",
  "Intermediate HLA-G+",   "HLA-DRA",   "HLA.DRA",
  "Intermediate HLA-G+",   "HLA-DPA1",  "HLA.DPA1",
  "Intermediate HLA-G+",   "HLA-DPB1",  "HLA.DPB1",
  "Intermediate HLA-G+",   "HLA-DRB1",  "HLA.DRB1",
  "Intermediate HLA-G+",   "AIM2",      "",
  "Atypical B",            "TBX21",     "",
  "Atypical B",            "ITGAX",     "",
  "Atypical B",            "FCRL5",     "",
  "Atypical B",            "FCRL3",     "",
  "Atypical B",            "ZEB2",      "",
  "Atypical B",            "LILRB1",    "",
  "Atypical B",            "LILRB2",    "",
  "Atypical B",            "CXCR3",     "",
  "Memory B",              "CD27",      "",
  "Memory B",              "COCH",      "",
  "Memory B",              "LTB",       "",
  "Memory B",              "CD40",      "",
  "Memory B",              "CD83",      "",
  "Memory B",              "TNFRSF13B", "",
  "Plasmablast exclusion", "JCHAIN",    "",
  "Plasmablast exclusion", "XBP1",      "",
  "Plasmablast exclusion", "MZB1",      "",
  "Plasmablast exclusion", "PRDM1",     ""
)
rna_cat_levels <- c("Lineage", "Naive B (resting)", "Transitional B", "IEG-High",
                    "Intermediate HLA-G+", "Atypical B", "Memory B",
                    "Plasmablast exclusion")

# =============================================================================
# (1) GREEN HEATMAP BUILDER  (CD8-manuscript palette, pheatmap, grouped rows)
#     ADT averaging is DSB-safe (manual rowMeans, NO expm1) per project rule;
#     RNA uses AverageExpression on log1p "data".
# =============================================================================
build_green_heatmap <- function(obj, assay_name, panel_df, group_col,
                                clust_levels = NULL, title = "", out_png,
                                cellwidth = 55, cellheight = 16) {
  DefaultAssay(obj) <- assay_name
  available <- rownames(obj[[assay_name]])
  panel_df  <- panel_df[panel_df$feature %in% available, , drop = FALSE]
  if (nrow(panel_df) < 2) { message("  <2 features for ", assay_name, " / ", title, " — skipped"); return(invisible()) }
  
  grp <- droplevels(factor(obj@meta.data[[group_col]]))
  if (is.null(clust_levels)) clust_levels <- levels(grp)
  clust_levels <- clust_levels[clust_levels %in% levels(grp)]
  
  if (assay_name == "ADT") {
    dat <- GetAssayData(obj, assay = "ADT", slot = "data")[panel_df$feature, , drop = FALSE]
    avg <- sapply(clust_levels, function(cl) rowMeans(dat[, which(grp == cl), drop = FALSE]))
  } else {
    avgL <- AverageExpression(obj, assays = assay_name, features = panel_df$feature,
                              group.by = group_col, slot = "data")[[assay_name]]
    ci  <- match(.norm(clust_levels), .norm(colnames(avgL)))
    if (anyNA(ci)) ci <- seq_along(clust_levels)
    avg <- avgL[panel_df$feature, ci, drop = FALSE]
  }
  rownames(avg) <- panel_df$label; colnames(avg) <- clust_levels
  avg_scaled <- scale_01(avg)
  
  annot_df  <- data.frame(Group = factor(panel_df$category, levels = unique(panel_df$category)),
                          row.names = panel_df$label)
  ro        <- cluster_within_groups(avg_scaled, annot_df)
  avg_ord   <- avg_scaled[ro, , drop = FALSE]
  annot_ord <- annot_df[ro, , drop = FALSE]
  gaps      <- get_gaps(annot_ord)
  
  glv <- levels(droplevels(annot_ord$Group))
  gc  <- bcell_cat_palette[glv]; gc[is.na(gc)] <- "#BDBDBD"; names(gc) <- glv
  annot_colors <- list(Group = gc)
  
  cnts <- as.integer(table(grp)[clust_levels])
  col_labels <- paste0(clust_levels, " (n=", cnts, ")")
  
  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
  pheatmap(
    avg_ord, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none",
    color = heatmap_colors_green, border_color = "white",
    annotation_row = annot_ord, annotation_colors = annot_colors,
    annotation_names_row = FALSE, gaps_row = gaps,
    labels_col = col_labels,
    cellwidth = cellwidth, cellheight = cellheight,
    fontsize = 11, fontsize_row = 10, fontsize_col = 10, angle_col = 45,
    main = title, filename = out_png
  )
  message("  saved: ", out_png)
}

if (RUN_GREEN_HEATMAPS) {
  message("\n===== (1) Green heatmaps =====")
  
  # ---- 1a. Final merged (Bcell_Final): alongside the existing magma ones -----
  if (file.exists(file.path(load.path, "OPIS_BCELL_FinalMerged.qs2"))) {
    OBF <- qs_read(file.path(load.path, "OPIS_BCELL_FinalMerged.qs2"))
    final.levels <- levels(droplevels(factor(OBF$Bcell_Final)))
    out.dir <- file.path(final.dir, "heatmaps", "green")
    build_green_heatmap(OBF, "ADT", prepare_panel(adt_markers, adt_cat_levels, rownames(OBF[["ADT"]]), "ADT"),
                        "Bcell_Final", final.levels, "B cells | ADT panel (final merged)",
                        file.path(out.dir, "Bcell_FINAL_ADT_heatmap_GREEN.png"))
    build_green_heatmap(OBF, "RNA", prepare_panel(rna_markers, rna_cat_levels, rownames(OBF[["RNA"]]), "RNA"),
                        "Bcell_Final", final.levels, "B cells | RNA panel (final merged)",
                        file.path(out.dir, "Bcell_FINAL_RNA_heatmap_GREEN.png"))
  } else message("  OPIS_BCELL_FinalMerged.qs2 not found — skipping final-merged green heatmaps.")
  
  # ---- 1b. Annotated (Bcell_Annotation, 7 labels) ----------------------------
  if (file.exists(file.path(load.path, "OPIS_BCELL_Annotated.qs2"))) {
    OBA <- qs_read(file.path(load.path, "OPIS_BCELL_Annotated.qs2"))
    ann.levels <- levels(droplevels(factor(OBA$Bcell_Annotation)))
    out.dir <- file.path(bcell.save.path, "heatmaps_green_annotated")
    build_green_heatmap(OBA, "ADT", prepare_panel(adt_markers, adt_cat_levels, rownames(OBA[["ADT"]]), "ADT"),
                        "Bcell_Annotation", ann.levels, "B cells | ADT panel (7-label annotation)",
                        file.path(out.dir, "Bcell_ANNOT_ADT_heatmap_GREEN.png"))
    build_green_heatmap(OBA, "RNA", prepare_panel(rna_markers, rna_cat_levels, rownames(OBA[["RNA"]]), "RNA"),
                        "Bcell_Annotation", ann.levels, "B cells | RNA panel (7-label annotation)",
                        file.path(out.dir, "Bcell_ANNOT_RNA_heatmap_GREEN.png"))
  } else message("  OPIS_BCELL_Annotated.qs2 not found — skipping annotated green heatmaps.")
  
  # ---- 1c. Pre-annotation subclusters (Subcluster_ID) ------------------------
  #   This is the "subcluster pre-annotation folder" heatmap, in green.
  if (file.exists(file.path(load.path, "OPIS_BCELL_PreAnnotation.qs2"))) {
    OBP <- qs_read(file.path(load.path, "OPIS_BCELL_PreAnnotation.qs2"))
    sub.levels <- as.character(sort(as.integer(levels(factor(OBP$Subcluster_ID)))))
    out.dir <- file.path(bcell.save.path, "heatmaps_green_subcluster")
    build_green_heatmap(OBP, "ADT", prepare_panel(adt_markers, adt_cat_levels, rownames(OBP[["ADT"]]), "ADT"),
                        "Subcluster_ID", sub.levels, "B cells | ADT panel (subclusters, pre-annotation)",
                        file.path(out.dir, "Bcell_SUBCLUSTER_ADT_heatmap_GREEN.png"))
    build_green_heatmap(OBP, "RNA", prepare_panel(rna_markers, rna_cat_levels, rownames(OBP[["RNA"]]), "RNA"),
                        "Subcluster_ID", sub.levels, "B cells | RNA panel (subclusters, pre-annotation)",
                        file.path(out.dir, "Bcell_SUBCLUSTER_RNA_heatmap_GREEN.png"))
  } else message("  OPIS_BCELL_PreAnnotation.qs2 not found — skipping subcluster green heatmaps.")
}

# =============================================================================
# (2) FREQUENCIES — % and raw counts for the final merged B clusters
#     -> final_merged/frequencies/ (overall, by OUD_status, by donor) + plots
# =============================================================================
if (RUN_FREQUENCIES) {
  message("\n===== (2) Cluster frequencies =====")
  freq.dir <- file.path(final.dir, "frequencies")
  dir.create(freq.dir, recursive = TRUE, showWarnings = FALSE)
  
  OBF <- qs_read(file.path(load.path, "OPIS_BCELL_FinalMerged.qs2"))
  grp <- droplevels(factor(OBF$Bcell_Final))
  lvls <- levels(grp)
  clust_cols <- setNames(c("#4C72B0", "#55A868", "#C44E52", "#8172B3", "#CCB974")[seq_along(lvls)], lvls)
  
  # --- overall ---
  overall <- tibble(Cluster = lvls,
                    n_cells = as.integer(table(grp)[lvls])) %>%
    mutate(percent = round(100 * n_cells / sum(n_cells), 2))
  write.csv(overall, file.path(freq.dir, "Bcell_freq_overall.csv"), row.names = FALSE)
  
  # --- by OUD_status (counts + % within each OUD group) ---
  if ("OUD_status" %in% colnames(OBF@meta.data)) {
    oud <- droplevels(factor(OBF$OUD_status))
    tab <- table(Cluster = grp, OUD_status = oud)
    cnt_df <- as.data.frame.matrix(tab) %>% rownames_to_column("Cluster")
    write.csv(cnt_df, file.path(freq.dir, "Bcell_freq_byOUD_counts.csv"), row.names = FALSE)
    pct_df <- as.data.frame.matrix(round(prop.table(tab, margin = 2) * 100, 2)) %>%
      rownames_to_column("Cluster")
    write.csv(pct_df, file.path(freq.dir, "Bcell_freq_byOUD_percentWithinOUD.csv"), row.names = FALSE)
    
    # stacked composition bar: cluster makeup within each OUD group
    comp <- as.data.frame(tab) %>% group_by(OUD_status) %>%
      mutate(percent = 100 * Freq / sum(Freq)) %>% ungroup() %>%
      mutate(Cluster = factor(Cluster, levels = lvls))
    p_comp <- ggplot(comp, aes(OUD_status, percent, fill = Cluster)) +
      geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
      scale_fill_manual(values = clust_cols) +
      labs(title = "B-cell composition by OUD status", x = NULL, y = "% of B cells") +
      theme_minimal(base_size = 13)
    ggsave(file.path(freq.dir, "Bcell_composition_byOUD_stacked.png"),
           p_comp, width = 7, height = 6, dpi = 500, bg = "white")
  } else message("  OUD_status not found — skipping OUD breakdown.")
  
  # --- by donor (orig.ident) ---
  if ("orig.ident" %in% colnames(OBF@meta.data)) {
    don <- droplevels(factor(OBF$orig.ident))
    tab_d <- table(Cluster = grp, Donor = don)
    write.csv(as.data.frame.matrix(tab_d) %>% rownames_to_column("Cluster"),
              file.path(freq.dir, "Bcell_freq_byDonor_counts.csv"), row.names = FALSE)
    write.csv(as.data.frame.matrix(round(prop.table(tab_d, margin = 2) * 100, 2)) %>%
                rownames_to_column("Cluster"),
              file.path(freq.dir, "Bcell_freq_byDonor_percentWithinDonor.csv"), row.names = FALSE)
  } else message("  orig.ident not found — skipping donor breakdown.")
  
  # --- overall counts bar plot ---
  p_cnt <- overall %>% mutate(Cluster = factor(Cluster, levels = lvls)) %>%
    ggplot(aes(Cluster, n_cells, fill = Cluster)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = paste0(n_cells, "\n(", percent, "%)")), vjust = -0.25, size = 3.5) +
    scale_fill_manual(values = clust_cols) +
    labs(title = "B-cell cluster sizes (final merged)", x = NULL, y = "Cells") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1)) +
    expand_limits(y = max(overall$n_cells) * 1.15)
  ggsave(file.path(freq.dir, "Bcell_cluster_counts_bar.png"),
         p_cnt, width = 8, height = 6, dpi = 500, bg = "white")
  
  message("  frequency tables + plots -> ", freq.dir)
}

# =============================================================================
# (3) DGE up AND down in ONE significant file (your standard MAST format).
#     Script 20 wrote DGE_MAST_full_allClusters.csv with only.pos = FALSE, so it
#     already contains BOTH directions. Its significant export was up-only purely
#     because of the avg_log2FC > 0.25 filter. We just re-threshold the same full
#     table with abs(avg_log2FC) > 0.25 — no MAST rerun.
# =============================================================================
if (RUN_DGE_UPDOWN) {
  message("\n===== (3) DGE up + down =====")
  full.csv <- file.path(dge.dir, "DGE_MAST_full_allClusters.csv")
  if (!file.exists(full.csv)) {
    message("  ", full.csv, " not found — run script 20 section 7 first; skipping.")
  } else {
    full <- read.csv(full.csv, check.names = FALSE, stringsAsFactors = FALSE)
    
    # Cluster order: known final order first, then any extras present
    final.order <- c("Transitional B", "Naive B", "Intermediate B (HLA-G+)",
                     "Memory B", "Atypical B (ABC/DN2)")
    present <- unique(as.character(full$cluster))
    lvls    <- c(final.order[final.order %in% present], setdiff(present, final.order))
    full$cluster <- factor(full$cluster, levels = lvls)
    
    # ONE significant table, both directions. Only change vs script 20:
    # abs(avg_log2FC) > 0.25 instead of > 0.25, so down genes are kept.
    sig <- full %>%
      filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25) %>%
      arrange(cluster, desc(avg_log2FC))
    write.csv(sig, file.path(dge.dir, "DGE_MAST_significant.csv"), row.names = FALSE)
    for (cl in lvls) {
      tag <- gsub("[^A-Za-z0-9]+", "_", cl)
      write.csv(sig %>% filter(cluster == cl),
                file.path(dge.dir, paste0("DGE_MAST_", tag, ".csv")), row.names = FALSE)
    }
    message("  significant=", nrow(sig),
            " (up=", sum(sig$avg_log2FC > 0), ", down=", sum(sig$avg_log2FC < 0), ")",
            "  -> DGE_MAST_significant.csv (single file, both directions)")
  }
}

# =============================================================================
# ENRICHR ENGINE (adapted from 8_EnrichR.R; takes a gene VECTOR)
# =============================================================================
databases <- c(
  "TRRUST_Transcription_Factors_2019", "ChEA_2022", "TRANSFAC_and_JASPAR_PWMs",
  "KEGG_2021_Human", "WikiPathways_2024_Human", "GO_Biological_Process_2023",
  "MSigDB_Hallmark_2020", "Panther_2016", "Reactome_2022", "BioPlanet_2019"
)
tf_databases      <- c("TRRUST_Transcription_Factors_2019", "ChEA_2022", "TRANSFAC_and_JASPAR_PWMs")
pathway_databases <- setdiff(databases, tf_databases)

# enrichR must be ATTACHED (not just ::-namespaced): library() runs the hook that
# sets the base URL. Without it, enrichr() errors "Must specify ... url or handle".
# setEnrichrSite() then registers the live Enrichr endpoint.
init_enrichr <- function() {
  if (!requireNamespace("enrichR", quietly = TRUE) ||
      !requireNamespace("openxlsx", quietly = TRUE)) {
    message("  enrichR/openxlsx not installed — install.packages(c('enrichR','openxlsx'))")
    return(FALSE)
  }
  suppressPackageStartupMessages({ library(enrichR); library(openxlsx) })
  ok <- tryCatch({
    enrichR::setEnrichrSite("Enrichr")
    isTRUE(getOption("enrichR.live"))
  }, error = function(e) { message("  Enrichr site error: ", conditionMessage(e)); FALSE })
  if (!ok) message("  Enrichr unreachable (internet/proxy?) — pathway analysis will be skipped.")
  ok
}

run_enrichr_on_genes <- function(genes, label, base_output) {
  genes <- unique(genes[!is.na(genes) & nzchar(genes)])
  if (length(genes) < 5) { message("  <5 genes for ", label, " — skipped"); return(invisible()) }
  message("  - Enrichr: ", length(genes), " genes | ", label)
  enrichment <- enrichR::enrichr(genes, databases)
  
  dir.create(file.path(base_output, "CSVs"),  recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(base_output, "Plots"), recursive = TRUE, showWarnings = FALSE)
  
  wb <- openxlsx::createWorkbook()
  for (db in names(enrichment)) {
    openxlsx::addWorksheet(wb, substr(db, 1, 31))
    openxlsx::writeData(wb, substr(db, 1, 31), enrichment[[db]])
  }
  openxlsx::saveWorkbook(wb, file.path(base_output, "CSVs", paste0(label, "_Enrichment.xlsx")),
                         overwrite = TRUE)
  
  top_tf <- list(); top_path <- list()
  for (db_name in names(enrichment)) {
    dbr <- enrichment[[db_name]]
    if ("Combined Score" %in% colnames(dbr)) dbr <- dbr %>% rename(Combined.Score = `Combined Score`)
    if (!"Combined.Score" %in% colnames(dbr)) next
    sigr <- dbr %>% filter(Adjusted.P.value < 0.05)
    if (nrow(sigr) == 0) next
    tt <- sigr %>% arrange(desc(Combined.Score)) %>% slice_head(n = 10) %>% mutate(Database = db_name)
    if (db_name %in% tf_databases) top_tf[[db_name]] <- tt else top_path[[db_name]] <- tt
  }
  
  mk_bar <- function(dfl, ttl, ylab, png) {
    df <- bind_rows(dfl)
    if (!"Combined.Score" %in% colnames(df) || nrow(df) == 0) {
      message("  - no significant terms: ", ttl); return(invisible())
    }
    df <- df %>% arrange(desc(Combined.Score)) %>% slice_head(n = 20)
    p <- ggplot(df, aes(reorder(Term, Combined.Score), Combined.Score, fill = Database)) +
      geom_bar(stat = "identity") + scale_y_log10() + coord_flip() +
      labs(title = ttl, x = NULL, y = ylab) +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 11, color = "black"),
            axis.title = element_text(size = 14), plot.title = element_text(size = 15, hjust = 0.5),
            legend.title = element_text(size = 11), legend.text = element_text(size = 10))
    ggsave(png, p, width = 12, height = 10, dpi = 300, bg = "white")
    message("  - saved: ", png)
  }
  mk_bar(top_tf,   paste("Top Transcription Factors -", label), "log10(Combined Score)",
         file.path(base_output, "Plots", paste0(label, "_Transcription_Factors.png")))
  mk_bar(top_path, paste("Top Pathways -", label),             "log10(Combined Score)",
         file.path(base_output, "Plots", paste0(label, "_Pathways.png")))
}

# =============================================================================
# (4) ENRICHR on B-cell merged-cluster DGE (up AND down, per cluster)
# =============================================================================
if (RUN_ENRICHR_BCELL) {
  message("\n===== (4) EnrichR on B-cell DGE =====")
  if (!init_enrichr()) {
    message("  EnrichR unavailable — skipping B-cell pathway analysis.")
  } else {
    sig.csv <- file.path(dge.dir, "DGE_MAST_significant.csv")
    if (!file.exists(sig.csv)) {
      message("  ", sig.csv, " not found — run section (3) first; skipping.")
    } else {
      bcell.enr.out <- file.path(final.dir, "pathway", "EnrichR")
      sig <- read.csv(sig.csv, check.names = FALSE, stringsAsFactors = FALSE)
      # All significant genes per cluster (up + down together), same as the
      # broad-cluster EnrichR pattern in 8_EnrichR.R.
      for (cl in unique(sig$cluster)) {
        g   <- sig$gene[sig$cluster == cl]
        lab <- gsub("[^A-Za-z0-9]+", "_", cl)
        run_enrichr_on_genes(g, lab, bcell.enr.out)
      }
      message("  EnrichR (B cells) -> ", bcell.enr.out)
    }
  }
}

# =============================================================================
# (5) ENRICHR on broad-cluster DGE CSVs (B, CD4, ...) in Manuscript_Figures
#     Robust column detection: gene col or rownames; adj-p col; logFC col.
#     Runs up and down separately when a logFC column is present.
# =============================================================================
read_dge_csv <- function(csv) {
  df <- read.csv(csv, check.names = FALSE, stringsAsFactors = FALSE)
  # write.csv(..., row.names = TRUE) leaves the gene column with a blank header,
  # which dplyr refuses to operate on. Name blanks "gene" and de-duplicate.
  cn <- colnames(df)
  cn[is.na(cn) | cn == ""] <- "gene"
  colnames(df) <- make.unique(cn)
  
  gcol <- intersect(c("gene", "Gene", "genes", "Symbol", "symbol"), colnames(df))
  if (length(gcol)) {
    df$.gene <- as.character(df[[gcol[1]]])
  } else if (is.character(df[[1]]) && !all(grepl("^[0-9.]+$", df[[1]]))) {
    df$.gene <- as.character(df[[1]])                   # first col is gene symbols
  } else {
    df$.gene <- rownames(df)
  }
  pcol <- intersect(c("p_val_adj", "padj", "adj.P.Val", "FDR", "q_value", "qvalue", "Adjusted.P.value"),
                    colnames(df))
  lcol <- intersect(c("avg_log2FC", "avg_logFC", "log2FoldChange", "logFC", "LogFC"),
                    colnames(df))
  df$.padj <- if (length(pcol)) suppressWarnings(as.numeric(df[[pcol[1]]])) else NA_real_
  df$.lfc  <- if (length(lcol)) suppressWarnings(as.numeric(df[[lcol[1]]])) else NA_real_
  list(df = df, has_p = length(pcol) > 0, has_lfc = length(lcol) > 0)
}

if (RUN_ENRICHR_BROAD) {
  message("\n===== (5) EnrichR on broad-cluster DGE CSVs =====")
  if (!init_enrichr()) {
    message("  EnrichR unavailable — skipping broad-cluster pathway analysis.")
  } else if (!dir.exists(broad.csv.dir)) {
    message("  broad CSV dir not found: ", broad.csv.dir)
  } else {
    csvs <- list.files(broad.csv.dir, pattern = "\\.csv$", full.names = TRUE)
    if (length(csvs) == 0) message("  no CSVs in ", broad.csv.dir)
    for (csv in csvs) {
      base <- tools::file_path_sans_ext(basename(csv))     # e.g. "B", "CD4"
      message("\n  -> ", basename(csv))
      parsed <- tryCatch(read_dge_csv(csv), error = function(e) { message("     read failed: ", e$message); NULL })
      if (is.null(parsed)) next
      df <- parsed$df
      if (parsed$has_p) df <- df %>% filter(!is.na(.padj) & .padj < 0.05)
      if (nrow(df) == 0) { message("     no significant rows"); next }
      
      if (parsed$has_lfc) {
        run_enrichr_on_genes(df$.gene[df$.lfc >  0.25], paste0(base, "_UP"),   broad.enr.out)
        run_enrichr_on_genes(df$.gene[df$.lfc < -0.25], paste0(base, "_DOWN"), broad.enr.out)
      } else {
        run_enrichr_on_genes(df$.gene, base, broad.enr.out)   # no direction info
      }
    }
    message("  EnrichR (broad clusters) -> ", broad.enr.out)
  }
}

message("\nB_Cell_Updates complete.")

