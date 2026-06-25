# =============================================================================
# OPIS ECCITEseq — Monocytes: FINALISE A SINGLE RESOLUTION
#
# Runs markers + avg-expression + dendrogram/merge + signatures + curated dot
# plots for ONE chosen resolution only (NOT the whole sweep). Loads the
# pre-annotation object from the subclustering run (which already carries every
# wsnn_res.<r> column), so nothing is re-embedded or re-clustered.
#
# Pick the resolution by editing `chosen_res` below. Recommended: 0.4
#   - res 0.5 over-splits classical monocytes (merge_candidates: 0/1/3/5/8/10
#     differ by 13-44 DE genes; cluster 10 was mito-high low-quality).
#   - res >=0.6 is unstable on the clustree (cross-parent edges).
#   - res 0.4 resolves the CD16+ (FCGR3A/MS4A7) and inflammatory populations
#     while sitting just below the classical over-split; 0.3 is the conservative
#     alternative.
#
# Reads : OPIS_MONO_PreAnnotation.qs2
# Writes: subclustering/pre_annotation/Monocytes/selected_res<chosen_res>/
#         + OPIS_MONO_res<chosen_res>.qs2 (Subcluster_ID committed to that res)
# =============================================================================

library(Seurat)
library(SeuratExtend)
library(scCustomize)
library(qs2)
library(dplyr); library(tidyr); library(tibble)
library(ggplot2); library(patchwork); library(viridis)
library(pheatmap)

# ============================ CHOOSE RESOLUTION ==============================
chosen_res <- 0.4            # <-- the only knob. Must be one of the swept values.

# ============================ Paths ==========================================
load.path  <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/saved_R_data"
mono.base  <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/subclustering/pre_annotation/Monocytes"
out.base   <- file.path(mono.base, sprintf("selected_res%.1f", chosen_res))
mk.dir     <- file.path(out.base, "Markers")
ann.dir    <- file.path(out.base, "annotation_support")
ann.tab    <- file.path(ann.dir, "tables")
for (d in c(mk.dir, ann.tab)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ============================ Config (same as the sweep script) ==============
merge_corr_threshold <- 0.90
de_padj <- 0.05; de_lfc <- 0.5
RUN_SINGLER <- FALSE

myeloid_lineage_genes <- c("LYZ","S100A8","S100A9","LST1","CTSS","FCN1","TYROBP","FCER1G")
lineage_check_adt <- c("CD14","CD16","CD11b","CD33","CD64","HLA-DR","CX3CR1","CD11c","CD86")

mono_signatures <- list(
  Classical_Mono    = c("CD14","FCN1","S100A8","S100A9","VCAN","CCR2","LRP1","SELL"),
  Intermediate_Mono = c("CD14","FCGR3A","MS4A7","HLA-DRA","HLA-DRB1","CD74","LST1"),
  NonClassical_Mono = c("FCGR3A","MS4A7","CX3CR1","LST1","IFITM2","LILRB1","RHOC","C1QA","C1QB"),
  Inflammatory      = c("IL1B","CXCL8","NLRP3","CCL3","CCL4","PTGS2","TNF","NFKBIA","S100A8","S100A9"),
  IFN_ISG           = c("ISG15","IFIT1","IFIT2","IFIT3","IFI6","MX1","OAS1","RSAD2","IRF7"),
  Antigen_Presenting= c("HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","CD74","CIITA","CD86"),
  cDC1              = c("CLEC9A","XCR1","BATF3","IRF8","CADM1"),
  cDC2              = c("CD1C","FCER1A","CLEC10A","FCGR2B"),
  pDC               = c("LILRA4","IL3RA","GZMB","CLEC4C","TCF4","SERPINF1"),
  Proliferation     = c("MKI67","TOP2A","BIRC5","PCNA","CENPF"),
  Platelet_contam   = c("PPBP","PF4","GP9","ITGA2B","TUBB1")
)
mono_name_lookup <- c(
  Classical_Mono    = "Classical monocyte (CD14+CD16-)",
  Intermediate_Mono = "Intermediate monocyte (CD14+CD16+)",
  NonClassical_Mono = "Non-classical monocyte (CD16+)",
  Inflammatory      = "Inflammatory monocyte (IL1B/CXCL8)",
  IFN_ISG           = "IFN / ISG-responding monocyte",
  Antigen_Presenting= "Antigen-presenting (HLA-DR/CIITA high)",
  cDC1 = "cDC1 (CLEC9A+)", cDC2 = "cDC2 (CD1C+)", pDC = "pDC (LILRA4+)",
  Proliferation = "Proliferating myeloid", Platelet_contam = "EXCLUDE - platelet contamination"
)

.norm <- function(x) toupper(gsub("[-_. ]", "", x))
resolve_feature <- function(label, available) {
  if (label %in% available) return(label)
  hit <- available[.norm(available) == .norm(label)]; if (length(hit)) hit[1] else NA_character_
}

# ============================ Load + commit chosen resolution ================
message("Loading OPIS_MONO_PreAnnotation.qs2 ...")
OPIS_MONO <- qs_read(file.path(load.path, "OPIS_MONO_PreAnnotation.qs2"))
if (inherits(OPIS_MONO[["RNA"]], "Assay5"))
  OPIS_MONO[["RNA"]] <- tryCatch(JoinLayers(OPIS_MONO[["RNA"]]), error = function(e) OPIS_MONO[["RNA"]])

rc <- paste0("wsnn_res.", chosen_res)
if (!rc %in% colnames(OPIS_MONO@meta.data))
  stop("Column ", rc, " not in object. Available: ",
       paste(grep("^wsnn_res\\.", colnames(OPIS_MONO@meta.data), value = TRUE), collapse = ", "))

OPIS_MONO$Subcluster_ID <- factor(OPIS_MONO[[rc]][, 1],
                                  levels = sort(as.integer(levels(factor(OPIS_MONO[[rc]][, 1])))))
Idents(OPIS_MONO) <- "Subcluster_ID"
clusters <- levels(OPIS_MONO$Subcluster_ID)
group_var <- "Subcluster_ID"
message("Resolution ", chosen_res, " -> ", length(clusters), " clusters: ",
        paste(clusters, collapse = ", "))

# Chosen-resolution UMAP
DefaultAssay(OPIS_MONO) <- "RNA"
ggsave(file.path(out.base, sprintf("UMAP_clusters_res%.1f.png", chosen_res)),
       DimPlot2(OPIS_MONO, features = "Subcluster_ID", reduction = "wnn.umap",
                label = TRUE, box = TRUE, theme = theme_umap_arrows()) +
         ggtitle(sprintf("Monocytes — clusters @ res %.1f", chosen_res)),
       width = 9, height = 8, dpi = 300, bg = "white")

# ============================ 1. Markers + avg expression ====================
message("\n=== markers (RNA + ADT) ===")
DefaultAssay(OPIS_MONO) <- "RNA"
mk <- FindAllMarkers(OPIS_MONO, assay = "RNA", only.pos = TRUE,
                     min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
write.csv(mk, file.path(mk.dir, sprintf("FindAllMarkers_RNA_res%.2f.csv", chosen_res)), row.names = FALSE)
if (nrow(mk) > 0) {
  top_mk <- mk %>% group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 20, with_ties = FALSE) %>% ungroup()
  write.csv(top_mk, file.path(mk.dir, sprintf("top20_per_cluster_RNA_res%.2f.csv", chosen_res)), row.names = FALSE)
}
avg_rna <- AverageExpression(OPIS_MONO, assays = "RNA", layer = "data", group.by = group_var)$RNA
colnames(avg_rna) <- sub("^g", "", colnames(avg_rna))
write.csv(as.data.frame(avg_rna), file.path(mk.dir, sprintf("avg_expr_RNA_res%.2f.csv", chosen_res)))

DefaultAssay(OPIS_MONO) <- "ADT"
mk_adt <- FindAllMarkers(OPIS_MONO, assay = "ADT", only.pos = TRUE,
                         min.pct = 0.25, logfc.threshold = 0.10, verbose = FALSE)
write.csv(mk_adt, file.path(mk.dir, sprintf("FindAllMarkers_ADT_res%.2f.csv", chosen_res)), row.names = FALSE)
clv     <- OPIS_MONO$Subcluster_ID
adt_dat <- GetAssayData(OPIS_MONO, assay = "ADT", slot = "data")
avg_adt <- sapply(clusters, function(cl) rowMeans(adt_dat[, which(clv == cl), drop = FALSE]))   # DSB-safe
rownames(avg_adt) <- rownames(adt_dat)
write.csv(as.data.frame(avg_adt), file.path(mk.dir, sprintf("avg_expr_ADT_res%.2f.csv", chosen_res)))
DefaultAssay(OPIS_MONO) <- "RNA"
message("  RNA markers=", nrow(mk), ", ADT markers=", nrow(mk_adt))

# ============================ 2. Correlation + dendrogram + merge ============
message("\n=== merge candidates ===")
hvg  <- intersect(VariableFeatures(OPIS_MONO), rownames(OPIS_MONO))
avgH <- AverageExpression(OPIS_MONO, assays = "RNA", layer = "data",
                          features = hvg, group.by = group_var)$RNA
colnames(avgH) <- sub("^g", "", colnames(avgH))
avgH <- avgH[, order(as.integer(colnames(avgH))), drop = FALSE]
cormat <- cor(log1p(as.matrix(avgH)), method = "pearson")
pheatmap(cormat, main = sprintf("Monocytes res %.1f: cluster correlation", chosen_res),
         cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, number_format = "%.2f",
         filename = file.path(ann.dir, "cluster_correlation_heatmap.png"), width = 8, height = 7)
hc <- hclust(as.dist(1 - cormat), method = "average")
png(file.path(ann.dir, "cluster_dendrogram.png"), width = 1600, height = 1000, res = 200)
plot(hc, main = sprintf("Monocytes res %.1f: dendrogram (1 - correlation)", chosen_res), xlab = "", sub = "")
dev.off()

pairs <- which(upper.tri(cormat), arr.ind = TRUE)
cand  <- data.frame()
for (k in seq_len(nrow(pairs))) {
  i <- rownames(cormat)[pairs[k, 1]]; j <- colnames(cormat)[pairs[k, 2]]
  rr <- cormat[pairs[k, 1], pairs[k, 2]]
  if (rr < merge_corr_threshold) next
  dem  <- FindMarkers(OPIS_MONO, ident.1 = i, ident.2 = j, only.pos = FALSE,
                      min.pct = 0.1, logfc.threshold = 0.1, verbose = FALSE)
  n_de <- sum(dem$p_val_adj < de_padj & abs(dem$avg_log2FC) > de_lfc, na.rm = TRUE)
  cand <- rbind(cand, data.frame(cluster_a = i, cluster_b = j,
                                 correlation = round(rr, 3), n_DE_genes = n_de))
}
if (nrow(cand) > 0) cand <- cand[order(-cand$correlation, cand$n_DE_genes), ]
write.csv(cand, file.path(ann.tab, "merge_candidates.csv"), row.names = FALSE)
message("  ", nrow(cand), " merge-candidate pair(s) (corr >= ", merge_corr_threshold, ")")

# per-cluster QC
mtcol <- intersect(c("percent.mt","percent_mt"), colnames(OPIS_MONO@meta.data))[1]
if (is.na(mtcol)) { OPIS_MONO$percent_mt <- PercentageFeatureSet(OPIS_MONO, pattern = "^MT-"); mtcol <- "percent_mt" }
qc_tab <- OPIS_MONO@meta.data %>% group_by(Subcluster_ID) %>%
  summarise(n_cells = n(), med_nFeature = median(nFeature_RNA), med_nCount = median(nCount_RNA),
            med_pct_mt = round(median(.data[[mtcol]]), 2), .groups = "drop")
write.csv(qc_tab, file.path(ann.tab, "qc_per_cluster.csv"), row.names = FALSE)
print(qc_tab)

# ============================ 3. Signatures + heatmap ========================
message("\n=== signature scores ===")
sig_cols <- character(0)
for (sn in names(mono_signatures)) {
  feats <- intersect(mono_signatures[[sn]], rownames(OPIS_MONO))
  if (length(feats) < 2) { message("  skip ", sn, " (<2 genes)"); next }
  OPIS_MONO <- AddModuleScore(OPIS_MONO, features = list(feats), name = paste0("sig_", sn, "_"), seed = 42)
  sig_cols <- c(sig_cols, paste0("sig_", sn, "_1"))
}
sigmat <- OPIS_MONO@meta.data %>% group_by(Subcluster_ID) %>%
  summarise(across(all_of(sig_cols), mean), .groups = "drop") %>%
  column_to_rownames("Subcluster_ID") %>% as.matrix()
colnames(sigmat) <- sub("^sig_", "", sub("_1$", "", colnames(sigmat)))
write.csv(round(as.data.frame(sigmat), 4), file.path(ann.tab, "signature_scores_per_cluster.csv"))
pheatmap(t(sigmat), scale = "row", cluster_cols = FALSE, cluster_rows = TRUE,
         display_numbers = TRUE, number_format = "%.2f",
         main = sprintf("Monocytes res %.1f: signature scores (row-scaled)", chosen_res),
         filename = file.path(ann.dir, "signature_heatmap.png"), width = 9, height = 6)

# ============================ 4. Curated dot plots + heatmap =================
message("\n=== curated dot plots ===")
rna_long <- rbind(
  data.frame(category = "Lineage", gene = myeloid_lineage_genes, stringsAsFactors = FALSE),
  do.call(rbind, lapply(names(mono_signatures), function(ct)
    data.frame(category = ct, gene = mono_signatures[[ct]], stringsAsFactors = FALSE))))
rna_long <- rna_long[!duplicated(rna_long$gene), ]
rna_long <- rna_long[rna_long$gene %in% rownames(OPIS_MONO[["RNA"]]), ]
rna_long$category <- factor(rna_long$category, levels = c("Lineage", names(mono_signatures)))
rna_long <- rna_long[order(rna_long$category), ]

DefaultAssay(OPIS_MONO) <- "RNA"
p_dot_rna <- DotPlot(OPIS_MONO, features = rna_long$gene, cluster.idents = FALSE) + coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(face = "italic", size = 8)) +
  labs(title = sprintf("Monocytes res %.1f: curated RNA markers", chosen_res), x = NULL, y = "Cluster")
ggsave(file.path(ann.dir, "marker_dotplot_RNA.png"), p_dot_rna,
       width = 10, height = 0.18 * nrow(rna_long) + 3, dpi = 300, bg = "white", limitsize = FALSE)

avg_panel <- AverageExpression(OPIS_MONO, assays = "RNA", layer = "data",
                               features = rna_long$gene, group.by = group_var)$RNA
colnames(avg_panel) <- sub("^g", "", colnames(avg_panel))
avg_panel <- avg_panel[, order(as.integer(colnames(avg_panel))), drop = FALSE]
ann_row <- data.frame(Category = rna_long$category, row.names = rna_long$gene)
pheatmap(log1p(as.matrix(avg_panel)), scale = "row", cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_row = ann_row, fontsize_row = 7,
         main = sprintf("Monocytes res %.1f: curated RNA panel (row-scaled log avg)", chosen_res),
         filename = file.path(ann.dir, "marker_avgexpr_heatmap_RNA.png"),
         width = 9, height = 0.16 * nrow(rna_long) + 3)

DefaultAssay(OPIS_MONO) <- "ADT"
adt_dot <- unname(na.omit(vapply(lineage_check_adt, resolve_feature, character(1),
                                 available = rownames(OPIS_MONO))))
if (length(adt_dot)) {
  p_dot_adt <- DotPlot(OPIS_MONO, features = adt_dot, cluster.idents = FALSE) + coord_flip() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = sprintf("Monocytes res %.1f: curated ADT markers", chosen_res), x = NULL, y = "Cluster")
  ggsave(file.path(ann.dir, "marker_dotplot_ADT.png"), p_dot_adt,
         width = 9, height = 0.3 * length(adt_dot) + 3, dpi = 300, bg = "white", limitsize = FALSE)
}
DefaultAssay(OPIS_MONO) <- "RNA"

# ============================ 5. Suggested-names worksheet ===================
top20 <- tryCatch(read.csv(file.path(mk.dir, sprintf("top20_per_cluster_RNA_res%.2f.csv", chosen_res))),
                  error = function(e) NULL)
top_markers <- function(cl) {
  if (is.null(top20)) return("")
  paste(head(top20$gene[as.character(top20$cluster) == as.character(cl)], 6), collapse = ", ")
}
worksheet <- lapply(clusters, function(cl) {
  ord <- names(sort(sigmat[as.character(cl), ], decreasing = TRUE))
  data.frame(cluster = cl, n_cells = sum(OPIS_MONO$Subcluster_ID == cl),
             top_signature = ord[1], second_signature = ord[2],
             top_markers = top_markers(cl),
             suggested_name = unname(mono_name_lookup[ord[1]]), stringsAsFactors = FALSE)
}) %>% bind_rows()
write.csv(worksheet, file.path(ann.dir, "suggested_names.csv"), row.names = FALSE)
print(worksheet[, c("cluster","top_signature","suggested_name")])

# ============================ 5b. Proposed annotation + labeled UMAPs =========
# Apply the res-0.4 annotation worked out from the marker/ADT/signature review
# (see Monocyte_annotation_README). Cluster IDs are res-0.4-specific, so this
# block only applies when those exact clusters are present.
umap.pt.size <- 1.6        # <-- bump up for bigger points

mono.fine <- c(
  "0" = "Classical monocyte (HSP/stress)",
  "1" = "Classical monocyte (VCAN+/platelet-assoc)",
  "2" = "Inflammatory monocyte",
  "3" = "Classical monocyte (hypoxic/glycolytic)",
  "4" = "MHC-II-high antigen-presenting monocyte",
  "5" = "CD16+ monocyte (ISG-high)",
  "6" = "Classical monocyte (CCR2/IRF8-high) [REVIEW]",
  "7" = "HLA-G+ monocyte",
  "8" = "Non-classical monocyte (CD16++)",
  "9" = "EXCLUDE - T-cell doublet/contamination"
)
mono.group <- c(
  "0" = "Classical monocyte", "1" = "Classical monocyte", "3" = "Classical monocyte",
  "2" = "Inflammatory monocyte",
  "4" = "Antigen-presenting monocyte",
  "7" = "HLA-G+ monocyte",
  "5" = "Non-classical/CD16+ monocyte", "8" = "Non-classical/CD16+ monocyte",
  "6" = "REVIEW", "9" = "EXCLUDE"
)
fine.levels  <- c("Classical monocyte (HSP/stress)",
                  "Classical monocyte (VCAN+/platelet-assoc)",
                  "Classical monocyte (hypoxic/glycolytic)",
                  "Inflammatory monocyte", "MHC-II-high antigen-presenting monocyte",
                  "HLA-G+ monocyte", "CD16+ monocyte (ISG-high)",
                  "Non-classical monocyte (CD16++)",
                  "Classical monocyte (CCR2/IRF8-high) [REVIEW]",
                  "EXCLUDE - T-cell doublet/contamination")
group.levels <- c("Classical monocyte", "Inflammatory monocyte",
                  "Antigen-presenting monocyte", "HLA-G+ monocyte",
                  "Non-classical/CD16+ monocyte", "REVIEW", "EXCLUDE")

if (all(clusters %in% names(mono.fine))) {
  message("\n=== proposed annotation UMAPs ===")
  OPIS_MONO$Mono_fine  <- factor(unname(mono.fine [as.character(OPIS_MONO$Subcluster_ID)]),
                                 levels = intersect(fine.levels,  unname(mono.fine)))
  OPIS_MONO$Mono_group <- factor(unname(mono.group[as.character(OPIS_MONO$Subcluster_ID)]),
                                 levels = intersect(group.levels, unname(mono.group)))
  
  # Fine labels: too long to print on-plot -> legend only, bigger points
  ggsave(file.path(out.base, sprintf("UMAP_annotation_fine_res%.1f.png", chosen_res)),
         DimPlot2(OPIS_MONO, features = "Mono_fine", reduction = "wnn.umap",
                  pt.size = umap.pt.size, label = FALSE, theme = theme_umap_arrows()) +
           ggtitle(sprintf("Monocytes res %.1f — proposed annotation (fine)", chosen_res)),
         width = 12, height = 8, dpi = 300, bg = "white")
  
  # Merged groups: labelled on-plot, bigger points
  ggsave(file.path(out.base, sprintf("UMAP_annotation_group_res%.1f.png", chosen_res)),
         DimPlot2(OPIS_MONO, features = "Mono_group", reduction = "wnn.umap",
                  pt.size = umap.pt.size, label = TRUE, box = TRUE,
                  theme = theme_umap_arrows()) +
           ggtitle(sprintf("Monocytes res %.1f — proposed annotation (merged groups)", chosen_res)),
         width = 11, height = 8, dpi = 300, bg = "white")
  message("  saved annotated UMAPs (fine + merged groups).")
} else {
  message("\nSkipping proposed-annotation UMAPs: cluster IDs (",
          paste(clusters, collapse = ", "), ") don't match the res-0.4 map. ",
          "Edit mono.fine / mono.group for this resolution.")
}

# ============================ 5c. GREEN marker heatmaps (annotation confirm) ==
# CD8/B-cell manuscript style: grouped curated panels, green palette, row group
# annotation + gaps, columns = clusters (top bar = proposed group). RNA averaged
# on log1p "data"; ADT averaged DSB-safe (manual rowMeans, no expm1).
message("\n=== green marker heatmaps ===")

heatmap_colors_green <- colorRampPalette(
  c("#F7FCF5", "#C7E9C0", "#74C476", "#31A354", "#006D2C"))(100)
scale_01 <- function(m) {
  m <- as.matrix(m)
  out <- t(apply(m, 1, function(x) {
    rng <- range(x, na.rm = TRUE)
    if (!is.finite(diff(rng)) || diff(rng) == 0) return(rep(0, length(x)))
    (x - rng[1]) / (rng[2] - rng[1])
  })); dimnames(out) <- dimnames(m); out
}
get_gaps <- function(a) { g <- as.character(a$Category); which(g[-length(g)] != g[-1]) }
cluster_within_groups <- function(mat, a) {
  ro <- character(0)
  for (grp in levels(droplevels(a$Category))) {
    rows <- rownames(a)[a$Category == grp]
    if (length(rows) <= 1) { ro <- c(ro, rows); next }
    hc <- hclust(dist(mat[rows, , drop = FALSE]), method = "ward.D2"); ro <- c(ro, rows[hc$order])
  }
  ro
}
resolve_alias <- function(label, aliases, available) {
  cand <- c(label, if (nzchar(aliases)) trimws(strsplit(aliases, ",")[[1]]) else character(0))
  for (c0 in cand) if (c0 %in% available) return(c0)
  avn <- .norm(available)
  for (c0 in cand) { hit <- available[avn == .norm(c0)]; if (length(hit)) return(hit[1]) }
  NA_character_
}
prep_panel <- function(tb, cat_levels, available) {
  tb$feature <- mapply(resolve_alias, tb$label, tb$aliases, MoreArgs = list(available = available))
  miss <- tb$label[is.na(tb$feature)]
  if (length(miss)) message("  [panel] not found, skipped: ", paste(miss, collapse = ", "))
  tb <- tb[!is.na(tb$feature), ]
  tb$category <- factor(tb$category, levels = cat_levels)
  tb <- tb[order(tb$category), ]
  tb[!duplicated(tb$feature), ]
}

mono_cat_palette <- c(
  "Lineage" = "#90A4AE", "Classical" = "#42A5F5", "Intermediate" = "#7E57C2",
  "Non-classical" = "#AB47BC", "Inflammatory" = "#FF7043",
  "IFN-ISG" = "#FFCA28", "Antigen-presenting" = "#EC407A")
mono_group_palette <- c(
  "Classical monocyte" = "#42A5F5", "Inflammatory monocyte" = "#FF7043",
  "Antigen-presenting monocyte" = "#EC407A", "HLA-G+ monocyte" = "#AB47BC",
  "Non-classical/CD16+ monocyte" = "#26A69A", "REVIEW" = "#BDBDBD", "EXCLUDE" = "#9E9E9E")

rna_cat_levels <- c("Lineage","Classical","Intermediate","Non-classical",
                    "Inflammatory","IFN-ISG","Antigen-presenting")
rna_panel <- tibble::tribble(
  ~category,            ~label,     ~aliases,
  "Lineage","LYZ","",   "Lineage","S100A8","", "Lineage","S100A9","", "Lineage","LST1","",
  "Lineage","CTSS","",  "Lineage","FCN1","",   "Lineage","TYROBP","","Lineage","FCER1G","",
  "Classical","CD14","","Classical","VCAN","", "Classical","CCR2","","Classical","LRP1","",
  "Classical","SELL","",
  "Intermediate","FCGR3A","","Intermediate","MS4A7","","Intermediate","HLA-DRA","",
  "Intermediate","HLA-DRB1","","Intermediate","CD74","",
  "Non-classical","CX3CR1","","Non-classical","IFITM2","","Non-classical","LILRB1","",
  "Non-classical","RHOC","","Non-classical","C1QA","","Non-classical","C1QB","",
  "Inflammatory","IL1B","","Inflammatory","CXCL8","","Inflammatory","NLRP3","",
  "Inflammatory","CCL3","","Inflammatory","CCL4","","Inflammatory","PTGS2","",
  "Inflammatory","TNF","","Inflammatory","NFKBIA","",
  "IFN-ISG","ISG15","","IFN-ISG","IFIT1","","IFN-ISG","IFIT2","","IFN-ISG","IFIT3","",
  "IFN-ISG","IFI6","","IFN-ISG","MX1","","IFN-ISG","OAS1","","IFN-ISG","RSAD2","","IFN-ISG","IRF7","",
  "Antigen-presenting","HLA-DPA1","","Antigen-presenting","HLA-DPB1","",
  "Antigen-presenting","CIITA","","Antigen-presenting","CD86","")

adt_panel <- tibble::tribble(
  ~category,  ~label,    ~aliases,
  "Lineage","CD14","",    "Lineage","CD16","FCGR3A", "Lineage","CD11b","ITGAM",
  "Lineage","CD33","",    "Lineage","CD64","FCGR1A", "Lineage","HLA-DR","HLA-DRA,HLA.DRA,HLADR",
  "Lineage","CX3CR1","",  "Lineage","CD11c","ITGAX", "Lineage","CD86","")

# proposed-group column annotation (if section 5b ran)
col_annot <- NULL; col_annot_colors <- NULL
if ("Mono_group" %in% colnames(OPIS_MONO@meta.data)) {
  cg <- sapply(clusters, function(cl) as.character(OPIS_MONO$Mono_group[OPIS_MONO$Subcluster_ID == cl])[1])
  col_annot <- data.frame(Proposed = cg, row.names = clusters)
  col_annot_colors <- list(Proposed = mono_group_palette[intersect(names(mono_group_palette), unique(cg))])
}

build_green_heatmap <- function(obj, assay_name, panel_df, clust_levels, title, out_png,
                                cellwidth = 24, cellheight = 11) {
  DefaultAssay(obj) <- assay_name
  available <- rownames(obj[[assay_name]])
  panel_df  <- panel_df[panel_df$feature %in% available, , drop = FALSE]
  if (nrow(panel_df) < 2) { message("  <2 features for ", assay_name, " — skipped"); return(invisible()) }
  grp <- droplevels(factor(obj$Subcluster_ID)); clust_levels <- clust_levels[clust_levels %in% levels(grp)]
  if (assay_name == "ADT") {
    dat <- GetAssayData(obj, assay = "ADT", slot = "data")[panel_df$feature, , drop = FALSE]
    avg <- sapply(clust_levels, function(cl) rowMeans(dat[, which(grp == cl), drop = FALSE]))
  } else {
    avgL <- AverageExpression(obj, assays = assay_name, features = panel_df$feature,
                              group.by = "Subcluster_ID", layer = "data")[[assay_name]]
    colnames(avgL) <- sub("^g", "", colnames(avgL))
    avg <- avgL[panel_df$feature, clust_levels, drop = FALSE]
  }
  rownames(avg) <- panel_df$label; colnames(avg) <- clust_levels
  avg_scaled <- scale_01(avg)
  annot_df <- data.frame(Category = factor(panel_df$category, levels = unique(panel_df$category)),
                         row.names = panel_df$label)
  ro <- cluster_within_groups(avg_scaled, annot_df)
  avg_ord <- avg_scaled[ro, , drop = FALSE]; annot_ord <- annot_df[ro, , drop = FALSE]
  gaps <- get_gaps(annot_ord)
  clv <- levels(droplevels(annot_ord$Category)); cc <- mono_cat_palette[clv]; cc[is.na(cc)] <- "#BDBDBD"; names(cc) <- clv
  ann_colors <- c(list(Category = cc), if (!is.null(col_annot_colors)) col_annot_colors)
  cnts <- as.integer(table(grp)[clust_levels]); col_labels <- paste0(clust_levels, " (n=", cnts, ")")
  ca <- NULL
  if (!is.null(col_annot)) { ca <- col_annot[clust_levels, , drop = FALSE]; rownames(ca) <- col_labels }
  colnames(avg_ord) <- col_labels
  pheatmap(avg_ord, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none",
           color = heatmap_colors_green, border_color = "white",
           annotation_row = annot_ord, annotation_col = ca, annotation_colors = ann_colors,
           annotation_names_row = FALSE, gaps_row = gaps,
           cellwidth = cellwidth, cellheight = cellheight, fontsize = 10,
           fontsize_row = 8, fontsize_col = 9, angle_col = 45, main = title, filename = out_png)
  message("  saved: ", out_png)
}

build_green_heatmap(OPIS_MONO, "RNA",
                    prep_panel(rna_panel, rna_cat_levels, rownames(OPIS_MONO[["RNA"]])),
                    clusters, sprintf("Monocytes res %.1f | RNA markers (green)", chosen_res),
                    file.path(ann.dir, "marker_heatmap_RNA_GREEN.png"))
build_green_heatmap(OPIS_MONO, "ADT",
                    prep_panel(adt_panel, "Lineage", rownames(OPIS_MONO[["ADT"]])),
                    clusters, sprintf("Monocytes res %.1f | ADT markers (green)", chosen_res),
                    file.path(ann.dir, "marker_heatmap_ADT_GREEN.png"), cellheight = 16)
DefaultAssay(OPIS_MONO) <- "RNA"

# ============================ 6. Optional SingleR ============================
if (RUN_SINGLER) {
  if (!requireNamespace("SingleR", quietly = TRUE) || !requireNamespace("celldex", quietly = TRUE)) {
    message("  SingleR/celldex not installed — skipping.")
  } else {
    ref <- celldex::MonacoImmuneData()
    pred <- SingleR::SingleR(test = GetAssayData(OPIS_MONO, assay = "RNA", slot = "data"),
                             ref = ref, labels = ref$label.fine)
    OPIS_MONO$SingleR_label <- pred$labels
    ct <- as.matrix(table(OPIS_MONO$Subcluster_ID, OPIS_MONO$SingleR_label))
    write.csv(as.data.frame.matrix(ct), file.path(ann.tab, "SingleR_crosstab_counts.csv"))
    write.csv(round(as.data.frame.matrix(sweep(ct, 1, pmax(rowSums(ct), 1), "/")), 3),
              file.path(ann.tab, "SingleR_crosstab_prop.csv"))
  }
}

# ============================ 7. Save ========================================
qs_save(OPIS_MONO, file.path(load.path, sprintf("OPIS_MONO_res%.1f.qs2", chosen_res)))
message("\nDone @ res ", chosen_res, ". Outputs: ", out.base)
message("Saved object: ", file.path(load.path, sprintf("OPIS_MONO_res%.1f.qs2", chosen_res)))