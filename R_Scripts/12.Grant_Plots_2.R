############################################################

library(qs2)
library(Seurat)
library(dplyr)
library(stringr)
library(readxl)
library(ggplot2)
library(EnhancedVolcano)
library(SeuratExtend)
library(scCustomize)
library(ggpubr)
library(stringr)
library(enrichR)
# ---------------------------- #
# 0) Paths / Load OPIS object
# ---------------------------- #
setwd("/home/akshay-iyer/Documents/OPIS_ECCITEseq")

load.path <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/saved_R_data/"
opis_qs2  <- file.path(load.path, "OPIS_ALL_PostAnnotation_IFNModules.qs2")
OPIS_ALL  <- qs_read(opis_qs2)

# Enforce OUD order + colors
OPIS_ALL$OUD_status <- as.character(OPIS_ALL$OUD_status)
# If your labels are already OUD-/OUD+, keep them; otherwise comment this out
OPIS_ALL$OUD_status <- factor(OPIS_ALL$OUD_status, levels = c("OUD-", "OUD+"))
oud_cols <- c("OUD-" = "#81D9AC", "OUD+" = "#E17C7A")

# ---------------------------- #
# 1) Output dirs (from scratch)
# ---------------------------- #
out_base <- file.path(getwd(), "grant_plots")
out_ms   <- file.path(out_base, "Module_Scoring_CD4CD8", "ViolinPlots")
out_ccr5 <- file.path(out_base, "CCR5_CD4CD8", "ViolinPlots")
out_objs <- file.path(out_base, "Module_Scoring_CD4CD8", "Scored_Object")

dir.create(out_ms,   recursive = TRUE, showWarnings = FALSE)
dir.create(out_ccr5, recursive = TRUE, showWarnings = FALSE)
dir.create(out_objs, recursive = TRUE, showWarnings = FALSE)

# ---------------------------- #
# 2) CD4/CD8 cluster list (Manual_Annotation)
# ---------------------------- #
tcell_clusters <- c(
  "CD4+ Memory (Treg)",
  "CD4+ Naïve 1",
  "CD4+ Naïve 2",
  "CD4+ TCM",
  "CD8+ early NK-T-like cytotoxic",
  "CD8+ Early TEMRA",
  "CD8+ Late EM",
  "CD8+ Naïve ",
  "CD8+ NK-T-like cytotoxic",
  "CD8+ TEM (intermediate)",
  "CD8+ TEM (terminal)",
  "CD8+ TEMRA",
  "CD8+ Transitional memory"
)

cells_t <- colnames(OPIS_ALL)[OPIS_ALL$Manual_Annotation %in% tcell_clusters]
seu_t   <- subset(OPIS_ALL, cells = cells_t)

# Use annotated clusters for plotting order
seu_t$Manual_Annotation <- factor(seu_t$Manual_Annotation, levels = tcell_clusters)
Idents(seu_t) <- "Manual_Annotation"

# ---------------------------- #
# 3) Helpers
# ---------------------------- #
safe_id <- function(x) {
  y <- iconv(x, to = "ASCII//TRANSLIT")
  y <- gsub("[^A-Za-z0-9]+", "_", y, perl = TRUE)
  y <- gsub("_+", "_", y)
  y <- sub("^_", "", y); y <- sub("_$", "", y)
  y
}

safe_filename <- function(label) {
  x <- iconv(label, to = "ASCII//TRANSLIT")
  x <- gsub("[/\\\\]", "-", x)
  x <- gsub("[^A-Za-z0-9_. -]", "", x, perl = TRUE)
  x <- gsub("\\s+", "_", x)
  x
}

# Score ONE module into meta field "MS_<mod_id>" on a given object
score_module <- function(obj, genes, mod_id, assay_use = "RNA") {
  DefaultAssay(obj) <- assay_use
  present <- intersect(genes, rownames(obj[[assay_use]]))
  field   <- paste0("MS_", mod_id)
  
  if (length(present) == 0) {
    obj[[field]] <- NA_real_
    return(obj)
  }
  
  obj <- AddModuleScore(obj, features = list(present), name = paste0(field, "_"), assay = assay_use)
  tmp <- paste0(field, "_1")
  obj[[field]] <- obj[[tmp]]
  obj[[tmp]] <- NULL
  return(obj)
}

# SeuratExtend violin with OUD split + your colors
plot_vln2_oud <- function(obj, feature, out_png, title, split_by = "OUD_status") {
  # Requires Idents(obj) set to Manual_Annotation already
  p <- SeuratExtend::VlnPlot2(
    obj,
    features    = feature,
    cols        = oud_cols,
    split.by    = split_by,
    stat.method = "wilcox.test",
    show.mean   = TRUE,
    mean_colors = c("red", "blue")
  ) +
    ggtitle(title)
  
  ggsave(out_png, p, dpi = 500, width = 16, height = 7, bg = "white")
  invisible(p)
}

# ---------------------------- #
# 4) Modules (Savita + NEW activation module)
# ---------------------------- #
modules <- list(
  # NEW: activation module
  Activation_HLADRA_CD38 = c("HLA-DRA", "CD38"),
  
  # Savita modules
  Primary_IFN_I_module = c("ISG15","IFIT1","IFIT2","IFIT3","IFI6","IFI27","MX1","MX2","OAS1","OAS2","OAS3","BST2"),
  IFN_Effector_module  = c("ISG15","IFIT1","IFIT2","IFIT3","MX1","MX2","OAS1","OAS2","OAS3"),
  IFN_signalling_module = c("STAT1","STAT2","IRF7","IRF9","IRF1","IRF3"),
  IFN_negative_regulation_module = c("USP18","SOCS1","SOCS3"),
  CD8_NK_interface_genes = c("GZMB","PRF1","NKG7","GNLY","TYROBP"),
  Immune_stress_transcription_factors = c("JUN","FOS","FOSL2","ATF3","BATF"),
  IFN_linked_transcription_factors = c("STAT1","IRF1")
)

assay_use <- "RNA"

# ---------------------------- #
# 5) SCORING ONLY (CD4/CD8 cells only)
# ---------------------------- #
coverage <- data.frame()

for (mod_name in names(modules)) {
  genes  <- modules[[mod_name]]
  mod_id <- safe_id(mod_name)
  
  seu_t <- score_module(seu_t, genes, mod_id, assay_use = assay_use)
  field <- paste0("MS_", mod_id)
  
  present <- intersect(genes, rownames(seu_t[[assay_use]]))
  missing <- setdiff(genes, present)
  
  coverage <- rbind(
    coverage,
    data.frame(
      Module = mod_name,
      Meta_Field = field,
      Present = length(present),
      Missing = length(missing),
      Present_Genes = paste(present, collapse = ";"),
      Missing_Genes = paste(missing, collapse = ";"),
      stringsAsFactors = FALSE
    )
  )
}

write.csv(coverage, file.path(out_base, "Module_Scoring_CD4CD8", "ModuleScore_Coverage_CD4CD8.csv"), row.names = FALSE)

# Optional: store these module scores back onto OPIS_ALL (NA for non-CD4/CD8)
# Safely transfer module scores from seu_t back into OPIS_ALL meta.data
# (fills non-CD4/CD8 cells with NA)

stopifnot(all(colnames(seu_t) %in% colnames(OPIS_ALL)))

for (mod_name in names(modules)) {
  mod_id <- safe_id(mod_name)
  field  <- paste0("MS_", mod_id)
  
  # initialize full-length named vector
  v <- rep(NA_real_, ncol(OPIS_ALL))
  names(v) <- colnames(OPIS_ALL)
  
  # pull scores from subset (must be numeric vector)
  sc <- seu_t@meta.data[[field]]
  
  # assign by barcode
  v[colnames(seu_t)] <- as.numeric(sc)
  
  # write back to OPIS_ALL meta.data
  OPIS_ALL[[field]] <- v
}

# ---------------------------- #
# 6) PLOTTING ONLY (VIOLINS, CD4/CD8 only)
# ---------------------------- #

# A) Module-score violins (split by OUD)
for (mod_name in names(modules)) {
  mod_id <- safe_id(mod_name)
  field  <- paste0("MS_", mod_id)
  
  if (!field %in% colnames(seu_t@meta.data)) {
    message("[PLOT] Missing module field in CD4/CD8 subset (skipping): ", field)
    next
  }
  
  out_png <- file.path(out_ms, paste0("Vln_", safe_filename(mod_name), "_SplitBy_OUD.png"))
  plot_vln2_oud(
    seu_t,
    feature = field,
    out_png = out_png,
    title   = paste0(mod_name, " | CD4/CD8 only | Split by OUD_status"),
    split_by = "OUD_status"
  )
}

# B) CCR5 violins (split by OUD), CD4/CD8 only
if ("CCR5" %in% rownames(seu_t[["RNA"]])) {
  out_png <- file.path(out_ccr5, "Vln_CCR5_CD4CD8_SplitBy_OUD.png")
  plot_vln2_oud(
    seu_t,
    feature = "CCR5",
    out_png = out_png,
    title   = "CCR5 | CD4/CD8 only | Split by OUD_status",
    split_by = "OUD_status"
  )
} else {
  message("[CCR5] CCR5 not found in RNA assay rownames; skipping CCR5 violin.")
}

message("DONE. Outputs in:")
message("  Module violins: ", out_ms)
message("  CCR5 violins:   ", out_ccr5)
message("  Scored object:  ", out_objs)


############################ p 2 ################################################3
############################################################
# OPIS: ADT-defined Total CD8 / Total NK -> Module scoring
# (CODE ONLY UP TO MODULE SCORING; plotting/pathways later)
############################################################


# ---------------------------- #
# 1) Output dirs (fresh tree)
# ---------------------------- #
out_base <- file.path(getwd(), "grant_plots")
out_gate <- file.path(out_base, "ADT_Gating_TotalCD8_TotalNK")
out_ms   <- file.path(out_gate, "Module_Scoring")
dir.create(out_gate, recursive = TRUE, showWarnings = FALSE)
dir.create(out_ms,   recursive = TRUE, showWarnings = FALSE)

# ---------------------------- #
# 2) ADT marker names + cutoffs (YOUR PANEL)
# ---------------------------- #
# Total CD8: CD3+ CD8+
adt_cd3  <- "CD3D"
adt_cd8  <- "CD8A"

# Total NK: CD56+ CD16+ CD3negative
adt_cd56 <- "NCAM1"     # CD56
adt_cd16 <- "FCGR3A"    # CD16

# Optional: CD38+DR+ within Total CD8
adt_cd38 <- "CD38"
adt_dr   <- "HLA-DRA"   # DR proxy in your panel

# Cutoffs for positivity in ADT 'data' slot
# If ADT is CLR/DSB-like, 0 is a reasonable starting point.
cut_cd3  <- 0
cut_cd8  <- 0
cut_cd56 <- 0
cut_cd16 <- 0
cut_cd38 <- 0
cut_dr   <- 0

# ---------------------------- #
# 3) Helpers: safe_id + module scoring
# ---------------------------- #
safe_id <- function(x) {
  y <- iconv(x, to = "ASCII//TRANSLIT")
  y <- gsub("[^A-Za-z0-9]+", "_", y, perl = TRUE)
  y <- gsub("_+", "_", y)
  y <- sub("^_", "", y); y <- sub("_$", "", y)
  y
}

score_module <- function(obj, genes, mod_id, assay_use = "RNA") {
  DefaultAssay(obj) <- assay_use
  present <- intersect(genes, rownames(obj[[assay_use]]))
  field   <- paste0("MS_", mod_id)
  
  if (length(present) == 0) {
    obj[[field]] <- NA_real_
    return(obj)
  }
  
  obj <- AddModuleScore(obj, features = list(present), name = paste0(field, "_"), assay = assay_use)
  tmp <- paste0(field, "_1")
  obj[[field]] <- obj[[tmp]]
  obj[[tmp]] <- NULL
  return(obj)
}

# ---------------------------- #
# 4) ADT gating (Total CD8 / Total NK / CD38+DR+ CD8)
# ---------------------------- #
DefaultAssay(OPIS_ALL) <- "ADT"
adt_feats <- rownames(OPIS_ALL[["ADT"]])

need <- c(adt_cd3, adt_cd8, adt_cd56, adt_cd16, adt_cd38, adt_dr)
missing <- setdiff(need, adt_feats)
if (length(missing)) stop("Missing ADT features needed for gating: ", paste(missing, collapse = ", "))

adt_df <- FetchData(OPIS_ALL, vars = need, slot = "data")

# Total CD8 gate: CD3D+ & CD8A+
is_cd8_total <- (adt_df[[adt_cd3]] > cut_cd3) & (adt_df[[adt_cd8]] > cut_cd8)

# Total NK gate: NCAM1+ & FCGR3A+ & CD3D-
is_nk_total <- (adt_df[[adt_cd56]] > cut_cd56) &
  (adt_df[[adt_cd16]] > cut_cd16) &
  (adt_df[[adt_cd3]] <= cut_cd3)

# CD38+DR+ within Total CD8 (CD38+ & HLA-DRA+)
is_cd8_cd38dr <- is_cd8_total &
  (adt_df[[adt_cd38]] > cut_cd38) &
  (adt_df[[adt_dr]]   > cut_dr)

OPIS_ALL$Gate_TotalCD8_ADT          <- is_cd8_total
OPIS_ALL$Gate_TotalNK_ADT           <- is_nk_total
OPIS_ALL$Gate_CD8_CD38pos_DRpos_ADT <- is_cd8_cd38dr

message("ADT-gated counts:")
message("  Total CD8: ", sum(OPIS_ALL$Gate_TotalCD8_ADT, na.rm = TRUE))
message("  Total NK : ", sum(OPIS_ALL$Gate_TotalNK_ADT,  na.rm = TRUE))
message("  CD38+DR+ within CD8: ", sum(OPIS_ALL$Gate_CD8_CD38pos_DRpos_ADT, na.rm = TRUE))

# ---------------------------- #
# 5) Define RNA modules to score
# ---------------------------- #
modules <- list(
  Cytotoxicity_PRF1_GZMB = c("PRF1", "GZMB"),
  Cytotoxicity_Broad     = c("PRF1","GZMB","NKG7","GNLY","GZMH","CTSW"),
  FasL_FASLG             = c("FASLG"),
  
  Primary_IFN_I_module = c("ISG15","IFIT1","IFIT2","IFIT3","IFI6","IFI27","MX1","MX2","OAS1","OAS2","OAS3","BST2"),
  IFN_Effector_module  = c("ISG15","IFIT1","IFIT2","IFIT3","MX1","MX2","OAS1","OAS2","OAS3"),
  IFN_signalling_module = c("STAT1","STAT2","IRF7","IRF9","IRF1","IRF3"),
  IFN_negative_regulation_module = c("USP18","SOCS1","SOCS3")
)

# ---------------------------- #
# 6) Score modules on RNA for the gated populations
# ---------------------------- #
DefaultAssay(OPIS_ALL) <- "RNA"

cells_cd8 <- colnames(OPIS_ALL)[OPIS_ALL$Gate_TotalCD8_ADT %in% TRUE]
cells_nk  <- colnames(OPIS_ALL)[OPIS_ALL$Gate_TotalNK_ADT  %in% TRUE]

seu_cd8 <- subset(OPIS_ALL, cells = cells_cd8)
seu_nk  <- subset(OPIS_ALL, cells = cells_nk)

coverage <- data.frame()

score_many <- function(seu, pop_label) {
  for (nm in names(modules)) {
    mod_id <- safe_id(paste(pop_label, nm, sep = "__"))
    seu <- score_module(seu, modules[[nm]], mod_id, assay_use = "RNA")
    
    field <- paste0("MS_", mod_id)
    present <- intersect(modules[[nm]], rownames(seu[["RNA"]]))
    missing <- setdiff(modules[[nm]], present)
    
    coverage <<- rbind(
      coverage,
      data.frame(
        Population = pop_label,
        Module = nm,
        Meta_Field = field,
        Present = length(present),
        Missing = length(missing),
        Present_Genes = paste(present, collapse = ";"),
        Missing_Genes = paste(missing, collapse = ";"),
        stringsAsFactors = FALSE
      )
    )
  }
  seu
}

seu_cd8 <- score_many(seu_cd8, "TotalCD8")
seu_nk  <- score_many(seu_nk,  "TotalNK")

# ---------------------------- #
# 7) Write module scores back to OPIS_ALL (NA for non-gated)
# ---------------------------- #
write_back_scores <- function(parent, child) {
  stopifnot(all(colnames(child) %in% colnames(parent)))
  score_cols <- grep("^MS_", colnames(child@meta.data), value = TRUE)
  
  for (field in score_cols) {
    v <- rep(NA_real_, ncol(parent))
    names(v) <- colnames(parent)
    v[colnames(child)] <- as.numeric(child@meta.data[[field]])
    parent[[field]] <- v
  }
  parent
}

OPIS_ALL <- write_back_scores(OPIS_ALL, seu_cd8)
OPIS_ALL <- write_back_scores(OPIS_ALL, seu_nk)

# ---------------------------- #
# 8) Save coverage + scored object
# ---------------------------- #
write.csv(coverage, file.path(out_ms, "ModuleScore_Coverage_TotalCD8_TotalNK.csv"), row.names = FALSE)

# Output folder (within your existing out_gate)
out_gate_fp <- file.path(out_gate, "Gating_FeaturePlots")
dir.create(out_gate_fp, recursive = TRUE, showWarnings = FALSE)

# Ensure these are numeric 0/1 for plotting intensity
OPIS_ALL$Gate_TotalCD8_ADT_num          <- as.numeric(OPIS_ALL$Gate_TotalCD8_ADT %in% TRUE)
OPIS_ALL$Gate_TotalNK_ADT_num           <- as.numeric(OPIS_ALL$Gate_TotalNK_ADT %in% TRUE)
OPIS_ALL$Gate_CD8_CD38pos_DRpos_ADT_num <- as.numeric(OPIS_ALL$Gate_CD8_CD38pos_DRpos_ADT %in% TRUE)

pal_gate <- viridis::viridis(n = 10, option = "A")

# Total CD8
p_cd8 <- scCustomize::FeaturePlot_scCustom(
  OPIS_ALL,
  reduction  = "wnn.umap",
  features   = "Gate_TotalCD8_ADT_num",
  colors_use = pal_gate,
  order      = TRUE
) + ggtitle("ADT Gate: Total CD8 (CD3D+ CD8A+)")

ggsave(file.path(out_gate_fp, "Gate_TotalCD8_ADT_FeaturePlot.png"),
       p_cd8, dpi = 500, width = 10, height = 8, bg = "white")

# Total NK
p_nk <- scCustomize::FeaturePlot_scCustom(
  OPIS_ALL,
  reduction  = "wnn.umap",
  features   = "Gate_TotalNK_ADT_num",
  colors_use = pal_gate,
  order      = TRUE
) + ggtitle("ADT Gate: Total NK (NCAM1+ FCGR3A+ CD3D-)")

ggsave(file.path(out_gate_fp, "Gate_TotalNK_ADT_FeaturePlot.png"),
       p_nk, dpi = 500, width = 10, height = 8, bg = "white")

# CD38+DR+ within CD8
p_cd8_act <- scCustomize::FeaturePlot_scCustom(
  OPIS_ALL,
  reduction  = "wnn.umap",
  features   = "Gate_CD8_CD38pos_DRpos_ADT_num",
  colors_use = pal_gate,
  order      = TRUE
) + ggtitle("ADT Gate: CD8 CD38+ HLA-DRA+ (within Total CD8)")

ggsave(file.path(out_gate_fp, "Gate_CD8_CD38pos_DRpos_ADT_FeaturePlot.png"),
       p_cd8_act, dpi = 500, width = 10, height = 8, bg = "white")
############################################################
# Module-score VIOLIN plots (SeuratExtend::VlnPlot2)
# - Use ADT gates as the grouping on x-axis (NOT Manual_Annotation)
# - Split by OUD_status with colors:
#     OUD+ = #E17C7A
#     OUD- = #81D9AC
############################################################

# ---- OUD colors ----
oud_cols <- c("OUD+" = "#E17C7A", "OUD-" = "#81D9AC")

# ---- output ----
out_vln <- file.path(out_gate, "Module_Scoring", "ViolinPlots_ByGate")
dir.create(out_vln, recursive = TRUE, showWarnings = FALSE)

safe_filename <- function(label) {
  x <- iconv(label, to = "ASCII//TRANSLIT")
  x <- gsub("[/\\\\]", "-", x)
  x <- gsub("[^A-Za-z0-9_. -]", "", x, perl = TRUE)
  x <- gsub("\\s+", "_", x)
  x
}

# If you DON'T already have this in your script, include it:
safe_id <- function(x) {
  y <- iconv(x, to = "ASCII//TRANSLIT")
  y <- gsub("[^A-Za-z0-9]+", "_", y, perl = TRUE)
  y <- gsub("_+", "_", y)
  y <- sub("^_", "", y); y <- sub("_$", "", y)
  y
}

plot_vln2_oud <- function(obj, feature, out_png, title,
                          split_by = "OUD_status",
                          min_n_per_group = 10) {
  
  if (!split_by %in% colnames(obj@meta.data)) {
    message("[Vln] Missing split.by column: ", split_by, " | skipping: ", feature)
    return(invisible(NULL))
  }
  
  # enforce OUD- then OUD+ ordering if present
  lv <- unique(as.character(obj[[split_by]][, 1]))
  if (all(c("OUD-", "OUD+") %in% lv)) {
    obj[[split_by]][, 1] <- factor(as.character(obj[[split_by]][, 1]), levels = c("OUD-", "OUD+"))
  }
  
  # Count non-missing per split group (module scores live in meta.data)
  md <- obj@meta.data
  do_stats <- TRUE
  n_by <- integer(0)
  
  if (feature %in% colnames(md)) {
    dfc <- md[, c(split_by, feature), drop = FALSE]
    colnames(dfc) <- c(".split", ".val")
    dfc <- dfc[!is.na(dfc$.split) & !is.na(dfc$.val), , drop = FALSE]
    n_by <- table(as.character(dfc$.split))
    do_stats <- (length(n_by) >= 2) && all(n_by >= min_n_per_group)
  } else {
    # if feature is gene expression, we won't prefilter by feature
    n_by <- table(as.character(md[[split_by]]))
    do_stats <- (length(n_by) >= 2) && all(n_by >= min_n_per_group)
  }
  
  if (!do_stats) {
    message("[Vln] Skipping Wilcox (low N) for ", feature, " | counts: ",
            paste(names(n_by), as.integer(n_by), sep = "=", collapse = ", "))
  }
  
  # ---- IMPORTANT: VlnPlot2 can't take stat.method=NULL; use two different calls ----
  if (do_stats) {
    p <- SeuratExtend::VlnPlot2(
      obj,
      features     = feature,
      cols         = oud_cols,
      split.by     = split_by,
      stat.method  = "wilcox.test",
      show.mean    = TRUE,
      mean_colors  = c("red", "blue")
    )
  } else {
    # no stat.method argument at all
    p <- SeuratExtend::VlnPlot2(
      obj,
      features     = feature,
      cols         = oud_cols,
      split.by     = split_by,
      show.mean    = TRUE,
      mean_colors  = c("red", "blue")
    )
  }
  
  p <- p + ggtitle(title)
  
  ggsave(out_png, p, dpi = 500, width = 16, height = 7, bg = "white")
  invisible(p)
}

# ------------------------------------------------------------
# Make gated objects (if you haven't already)
# ------------------------------------------------------------
cells_cd8 <- colnames(OPIS_ALL)[OPIS_ALL$Gate_TotalCD8_ADT %in% TRUE]
cells_nk  <- colnames(OPIS_ALL)[OPIS_ALL$Gate_TotalNK_ADT  %in% TRUE]

seu_cd8 <- subset(OPIS_ALL, cells = cells_cd8)
seu_nk  <- subset(OPIS_ALL, cells = cells_nk)

# Set x-axis grouping to a constant gate label (still split by OUD)
seu_cd8$ADT_Gate_Label <- "TotalCD8"
seu_nk$ADT_Gate_Label  <- "TotalNK"
Idents(seu_cd8) <- "ADT_Gate_Label"
Idents(seu_nk)  <- "ADT_Gate_Label"

# ------------------------------------------------------------
# Find module score fields robustly
#   Strategy: for each module name, look for columns that contain its safe_id()
# ------------------------------------------------------------
get_ms_fields <- function(seu, modules_named_list) {
  ms_cols <- grep("^MS_", colnames(seu@meta.data), value = TRUE)
  if (!length(ms_cols)) return(character(0))
  
  out <- c()
  for (nm in names(modules_named_list)) {
    mid <- safe_id(nm)
    hit <- ms_cols[grepl(mid, ms_cols, fixed = TRUE)]
    if (length(hit)) out <- c(out, hit[1])  # take first deterministic match
  }
  unique(out)
}

ms_fields_cd8 <- get_ms_fields(seu_cd8, modules)
ms_fields_nk  <- get_ms_fields(seu_nk,  modules)

# ------------------------------------------------------------
# Identify module-score columns that exist AND are non-NA
# ------------------------------------------------------------
ms_cols_cd8 <- grep("^MS_", colnames(seu_cd8@meta.data), value = TRUE)
ms_cols_nk  <- grep("^MS_", colnames(seu_nk@meta.data),  value = TRUE)

# Keep only those with at least some non-NA values
ms_cols_cd8 <- ms_cols_cd8[sapply(ms_cols_cd8, function(x) sum(!is.na(seu_cd8@meta.data[[x]])) > 0)]
ms_cols_nk  <- ms_cols_nk[sapply(ms_cols_nk,  function(x) sum(!is.na(seu_nk@meta.data[[x]]))  > 0)]

message("[DEBUG] CD8 MS_ cols with data: ", length(ms_cols_cd8))
message("[DEBUG] NK  MS_ cols with data: ", length(ms_cols_nk))


# ------------------------------------------------------------
# Plot CD8
# ------------------------------------------------------------
out_cd8 <- file.path(out_vln, "TotalCD8")
dir.create(out_cd8, recursive = TRUE, showWarnings = FALSE)

for (field in ms_cols_cd8) {
  out_png <- file.path(out_cd8, paste0("Vln_TotalCD8_", safe_filename(field), "_SplitBy_OUD.png"))
  plot_vln2_oud(
    seu_cd8,
    feature = field,
    out_png = out_png,
    title   = paste0("ADT Gate: Total CD8 | ", field, " | Split by OUD_status")
  )
}
# ------------------------------------------------------------
# Plot NK (only columns with data)
# ------------------------------------------------------------
out_nk <- file.path(out_vln, "TotalNK")
dir.create(out_nk, recursive = TRUE, showWarnings = FALSE)

for (field in ms_cols_nk) {
  out_png <- file.path(out_nk, paste0("Vln_TotalNK_", safe_filename(field), "_SplitBy_OUD.png"))
  plot_vln2_oud(
    seu_nk,
    feature = field,
    out_png = out_png,
    title   = paste0("ADT Gate: Total NK | ", field, " | Split by OUD_status")
  )
}

############################################################
# OPIS: TOTALCD8 + TOTALNK (NEW ADT GATES) — 4 STEPS
#   1) DGE
#   2) EnrichR (write xlsx with DB tabs)
#   3) EnrichR curated plots: Pathways + TFs (colored by database)
#
# Assumptions:
#   - OPIS_ALL already loaded
#   - OPIS_ALL@meta.data has:
#       OUD_status  (values like "OUD+" / "OUD-")
#       Gate_TotalCD8_ADT  (TRUE/FALSE)
#       Gate_TotalNK_ADT   (TRUE/FALSE)
#
# Output root:
#   /home/akshay-iyer/Documents/OPIS_ECCITEseq/grant_plots/ADT_Gates_TotalCD8_TotalNK
############################################################

# ---------------------------- #
# Paths / output base
# ---------------------------- #
base_dir   <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq"
grant_base <- file.path(base_dir, "grant_plots", "ADT_Gates_TotalCD8_TotalNK")
dir.create(grant_base, recursive = TRUE, showWarnings = FALSE)

# curated workbook (your existing)
curated_xlsx <- file.path(base_dir, "saved_R_data", "Curated genes for figures.xlsx")
stopifnot(file.exists(curated_xlsx))

# Use these DBs (edit as needed)
enrich_dbs <- c(
  "TRRUST_Transcription_Factors_2019",
  "ChEA_2022",
  "TRANSFAC_and_JASPAR_PWMs",
  "KEGG_2021_Human",
  "WikiPathways_2024_Human",
  "GO_Biological_Process_2023",
  "MSigDB_Hallmark_2020",
  "Panther_2016",
  "Reactome_2022",
  "BioPlanet_2019"
)

# ---------------------------- #
# Helper: sanitize names
# ---------------------------- #
safe_id <- function(x) {
  y <- iconv(x, to = "ASCII//TRANSLIT")
  y <- gsub("[^A-Za-z0-9]+", "_", y, perl = TRUE)
  y <- gsub("_+", "_", y)
  y <- sub("^_", "", y); y <- sub("_$", "", y)
  y
}

# ---------------------------- #
# Helper: curated term+db reader (no header; col1 term; col2 db)
# ---------------------------- #
.read_curated_term_db <- function(xlsx, sheet) {
  dat <- readxl::read_excel(xlsx, sheet = sheet, col_names = FALSE)
  dat <- dat %>%
    transmute(
      term = as.character(.data[[1]]),
      db   = as.character(.data[[2]])
    ) %>%
    filter(!is.na(term), term != "", !is.na(db), db != "")
  dat
}

# ---------------------------- #
# Helper: DB tab matching (handles "TRRUST..._20" vs "...2019" etc.)
# ---------------------------- #
.match_enrichr_tab <- function(db_i, tabs, xlsx_path = NULL) {
  db_raw <- stringr::str_trim(as.character(db_i))
  if (is.na(db_raw) || db_raw == "") return(NA_character_)
  
  # FIRST: if we have the mapping tab, use it
  if (!is.null(xlsx_path) && file.exists(xlsx_path) && "DB_SHEET_MAP" %in% tabs) {
    mp <- tryCatch(readxl::read_excel(xlsx_path, sheet = "DB_SHEET_MAP"), error = function(e) NULL)
    if (!is.null(mp) && all(c("Database","Sheet") %in% names(mp))) {
      hit <- mp$Sheet[match(db_raw, mp$Database)]
      if (!is.na(hit) && hit %in% tabs) return(hit)
    }
  }
  
  # FALLBACK: old normalization matching
  tabs_trim <- stringr::str_trim(tabs)
  norm <- function(x) gsub("[^A-Za-z0-9]+", "", tolower(stringr::str_trim(x)))
  db_norm   <- norm(db_raw)
  tabs_norm <- norm(tabs_trim)
  
  idx <- which(tabs_norm == db_norm); if (length(idx)) return(tabs_trim[idx[1]])
  idx <- which(startsWith(tabs_norm, db_norm)); if (length(idx)) return(tabs_trim[idx[1]])
  idx <- which(startsWith(db_norm, tabs_norm)); if (length(idx)) return(tabs_trim[idx[1]])
  idx <- which(grepl(db_norm, tabs_norm, fixed = TRUE)); if (length(idx)) return(tabs_trim[idx[1]])
  NA_character_
}

# ---------------------------- #
# Helper: read an EnrichR tab (Excel sheet)
# ---------------------------- #
.read_enrichr_tab <- function(xlsx_path, sheet) {
  df <- readxl::read_excel(xlsx_path, sheet = sheet)
  if (!"Term" %in% names(df)) stop("No 'Term' column in: ", basename(xlsx_path), " | tab: ", sheet)
  if (!"Combined.Score" %in% names(df) && "Combined score" %in% names(df)) df$Combined.Score <- df[["Combined score"]]
  df
}

# ---------------------------- #
# Helper: plot curated terms (Pathways or TFs), colored by Database column from curated sheet
# ---------------------------- #
.plot_curated_enrich_terms <- function(curated_df, xlsx_path, pop_name, out_png, plot_title) {
  if (!file.exists(xlsx_path)) {
    message("[ENRICH] Missing enrichment xlsx: ", xlsx_path)
    return(invisible(NULL))
  }
  
  tabs <- readxl::excel_sheets(xlsx_path)
  
  hits <- lapply(seq_len(nrow(curated_df)), function(i) {
    term_i <- curated_df$term[i]
    db_i   <- curated_df$db[i]
    
    tab <- .match_enrichr_tab(db_i, tabs)
    if (is.na(tab)) return(NULL)
    
    df <- tryCatch(.read_enrichr_tab(xlsx_path, tab), error = function(e) NULL)
    if (is.null(df)) return(NULL)
    
    df$Term_l <- str_to_lower(str_trim(df$Term))
    term_l    <- str_to_lower(str_trim(term_i))
    
    row <- df[df$Term_l == term_l, , drop = FALSE]
    if (!nrow(row)) return(NULL)
    
    score <- suppressWarnings(as.numeric(row$Combined.Score[1]))
    metric <- "Combined.Score"
    
    if (!is.finite(score) || is.na(score)) {
      if ("Adjusted.P.value" %in% names(row)) {
        ap <- suppressWarnings(as.numeric(row$Adjusted.P.value[1]))
        if (is.finite(ap) && !is.na(ap) && ap > 0) {
          score <- -log10(ap)
          metric <- "-log10(Adjusted.P.value)"
        }
      }
    }
    if (!is.finite(score) || is.na(score)) return(NULL)
    
    data.frame(Term = term_i, Database = db_i, Value = score, Metric = metric, stringsAsFactors = FALSE)
  })
  
  hit_df <- bind_rows(hits)
  if (!nrow(hit_df)) {
    message("[ENRICH] No curated matches for ", pop_name, " | ", basename(xlsx_path))
    return(invisible(NULL))
  }
  
  metric_label <- if (any(hit_df$Metric == "Combined.Score")) "Combined.Score" else unique(hit_df$Metric)[1]
  
  hit_df <- hit_df %>%
    arrange(desc(Value)) %>%
    mutate(
      Term = factor(Term, levels = rev(unique(Term))),
      Database = factor(Database, levels = sort(unique(Database)))
    )
  
  p <- ggplot(hit_df, aes(x = Term, y = Value, fill = Database)) +
    geom_col(color = "black", linewidth = 0.2) +
    coord_flip() +
    labs(title = plot_title, x = NULL, y = metric_label, fill = "Database") +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text  = element_text(size = 10)
    )
  
  ggsave(out_png, p, dpi = 500, width = 11, height = max(5, 0.25 * nrow(hit_df) + 2), bg = "white")
  message("[ENRICH] Saved: ", out_png)
  invisible(hit_df)
}

# ---------------------------- #
# Helper: EnrichR runner -> Excel with DB tabs
# ---------------------------- #
run_enrichr_to_xlsx <- function(gene_symbols, dbs, out_xlsx) {
  gene_symbols <- unique(gene_symbols[!is.na(gene_symbols) & gene_symbols != ""])
  if (!length(gene_symbols)) {
    message("[ENRICHR] No genes supplied for: ", basename(out_xlsx))
    return(invisible(NULL))
  }
  
  # --- FIX: ensure enrichR site options exist ---
  if (is.null(getOption("enrichR.sites.base.address"))) {
    options(enrichR.sites.base.address = "https://maayanlab.cloud/Enrichr/")
  }
  if (is.null(getOption("enrichR.sites"))) {
    options(enrichR.sites = list(Enrichr = "https://maayanlab.cloud/Enrichr/"))
  }
  
  # Set site safely (works across enrichR versions)
  try(enrichR::setEnrichrSite("Enrichr"), silent = TRUE)
  if (is.null(getOption("enrichR.site"))) {
    options(enrichR.site = "https://maayanlab.cloud/Enrichr/")
  }
  
  enr <- tryCatch(enrichR::enrichr(gene_symbols, dbs), error = function(e) e)
  if (inherits(enr, "error")) {
    message("[ENRICHR] Failed: ", enr$message)
    return(invisible(NULL))
  }
  
  # ---- FIX: shorten sheet names to <=31 and ensure uniqueness ----
  short_sheet <- function(x, max_len = 31) {
    s <- gsub("[\\[\\]\\*\\?/\\\\:]", "_", x)     # Excel-illegal chars
    if (nchar(s) > max_len) s <- substr(s, 1, max_len)
    s
  }
  sheet_names <- vapply(dbs, short_sheet, character(1))
  # make unique (TRRUST... may collide if multiple long names)
  sheet_names <- make.unique(sheet_names, sep = "_")
  
  # write workbook + a mapping tab
  wb <- openxlsx::createWorkbook()
  
  map_df <- data.frame(Database = dbs, Sheet = sheet_names, stringsAsFactors = FALSE)
  openxlsx::addWorksheet(wb, "DB_SHEET_MAP")
  openxlsx::writeData(wb, "DB_SHEET_MAP", map_df)
  
  for (i in seq_along(dbs)) {
    db <- dbs[i]
    sh <- sheet_names[i]
    df <- enr[[db]]
    if (is.null(df) || !nrow(df)) df <- data.frame()
    openxlsx::addWorksheet(wb, sheetName = sh)
    openxlsx::writeData(wb, sheet = sh, x = df)
  }
  
  openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)
  message("[ENRICHR] Saved: ", out_xlsx)
  invisible(TRUE)
}

# ============================================================
# STEP 1) DGE (TotalNK + TotalCD8)
# ============================================================
out_step1 <- file.path(grant_base, "1_DGE")
dir.create(out_step1, recursive = TRUE, showWarnings = FALSE)

stopifnot(all(c("OUD_status","Gate_TotalCD8_ADT","Gate_TotalNK_ADT") %in% colnames(OPIS_ALL@meta.data)))

# cells for each gate
cells_totalcd8 <- colnames(OPIS_ALL)[OPIS_ALL$Gate_TotalCD8_ADT %in% TRUE]
cells_totalnk  <- colnames(OPIS_ALL)[OPIS_ALL$Gate_TotalNK_ADT  %in% TRUE]

run_dge_gate <- function(obj, cells_use, pop_name,
                         assay_use = "RNA",
                         test_use = "MAST",
                         latent_vars = c("nCount_RNA"),
                         logfc_threshold = 0,
                         min_pct = 0.1) {
  
  seu <- subset(obj, cells = cells_use)
  DefaultAssay(seu) <- assay_use
  
  seu$OUD_status <- as.character(seu$OUD_status)
  # enforce order if present
  if (all(c("OUD-","OUD+") %in% unique(seu$OUD_status))) {
    seu$OUD_status <- factor(seu$OUD_status, levels = c("OUD-","OUD+"))
  } else {
    seu$OUD_status <- factor(seu$OUD_status)
  }
  
  gtab <- table(seu$OUD_status, useNA = "ifany")
  message("[DGE] ", pop_name, " | counts: ", paste(names(gtab), gtab, collapse = ", "))
  
  if (!all(c("OUD-","OUD+") %in% names(gtab)) || any(gtab[c("OUD-","OUD+")] < 3)) {
    message("[DGE] Skipping ", pop_name, " (need >=3 cells per group).")
    return(invisible(NULL))
  }
  
  Idents(seu) <- "OUD_status"
  
  dge <- FindMarkers(
    seu,
    ident.1 = "OUD+",
    ident.2 = "OUD-",
    assay = assay_use,
    test.use = test_use,
    latent.vars = latent_vars,
    logfc.threshold = logfc_threshold,
    min.pct = min_pct
  )
  
  dge$Feature <- rownames(dge)
  dge
}

dge_totalnk  <- run_dge_gate(OPIS_ALL, cells_totalnk,  "TotalNK")
dge_totalcd8 <- run_dge_gate(OPIS_ALL, cells_totalcd8, "TotalCD8")

# save DGE csvs
if (!is.null(dge_totalnk)) {
  write.csv(dge_totalnk, file.path(out_step1, "TotalNK_OUD_Pos_vs_Neg_RNA.csv"), row.names = FALSE)
}
if (!is.null(dge_totalcd8)) {
  write.csv(dge_totalcd8, file.path(out_step1, "TotalCD8_OUD_Pos_vs_Neg_RNA.csv"), row.names = FALSE)
}

# ============================================================
# STEP 2) EnrichR (write Excel with DB tabs)
#   - uses significant genes (padj<0.05 & abs(log2FC)>=0.25)
#   - keeps SYMBOLs only; drops ENSG + MT-/RPS/RPL for stability
# ============================================================
out_step2 <- file.path(grant_base, "2_EnrichR_XLSX")
dir.create(out_step2, recursive = TRUE, showWarnings = FALSE)

prep_genes_for_enrichr <- function(dge_df,
                                   padj_cut = 0.05,
                                   abs_fc_cut = 0.25,
                                   drop_housekeeping = TRUE) {
  if (is.null(dge_df) || !nrow(dge_df)) return(character(0))
  
  fc_col <- if ("avg_log2FC" %in% colnames(dge_df)) "avg_log2FC" else "avg_logFC"
  
  df <- dge_df %>%
    filter(!is.na(p_val_adj), is.finite(p_val_adj), p_val_adj < padj_cut) %>%
    filter(!is.na(.data[[fc_col]]), is.finite(.data[[fc_col]]), abs(.data[[fc_col]]) >= abs_fc_cut) %>%
    mutate(Feature = as.character(Feature))
  
  genes <- unique(df$Feature)
  # remove ENSG
  genes <- genes[!str_detect(genes, "^ENSG")]
  # remove housekeeping
  if (drop_housekeeping) genes <- genes[!str_detect(genes, "^(MT-|RPL|RPS)")]
  unique(genes)
}

# TotalNK enrichr
if (!is.null(dge_totalnk)) {
  genes_nk <- prep_genes_for_enrichr(dge_totalnk)
  message("[ENRICHR] TotalNK genes: ", length(genes_nk))
  run_enrichr_to_xlsx(
    gene_symbols = genes_nk,
    dbs = enrich_dbs,
    out_xlsx = file.path(out_step2, "TotalNK_OUD_Pos_vs_Neg_RNA_Enrichment.xlsx")
  )
}

# TotalCD8 enrichr
if (!is.null(dge_totalcd8)) {
  genes_cd8 <- prep_genes_for_enrichr(dge_totalcd8)
  message("[ENRICHR] TotalCD8 genes: ", length(genes_cd8))
  run_enrichr_to_xlsx(
    gene_symbols = genes_cd8,
    dbs = enrich_dbs,
    out_xlsx = file.path(out_step2, "TotalCD8_OUD_Pos_vs_Neg_RNA_Enrichment.xlsx")
  )
}

# ============================================================
# STEP 3) EnrichR plots (NOT curated) — top terms per DB
#   - TF DBs = first 3 in enrich_dbs
#   - Pathway DBs = rest
#   - Plots colored by Database
# ============================================================
out_step3 <- file.path(grant_base, "3_EnrichR_TopHits_Plots")
out_pw    <- file.path(out_step3, "Pathways")
out_tf    <- file.path(out_step3, "TFs")
dir.create(out_pw, recursive = TRUE, showWarnings = FALSE)
dir.create(out_tf, recursive = TRUE, showWarnings = FALSE)

tf_dbs <- enrich_dbs[1:3]
pw_dbs <- enrich_dbs[4:length(enrich_dbs)]

# read mapping (since we shortened sheet names to <=31 chars)
.read_db_sheet_map <- function(xlsx_path) {
  tabs <- readxl::excel_sheets(xlsx_path)
  if (!"DB_SHEET_MAP" %in% tabs) stop("DB_SHEET_MAP tab not found in: ", basename(xlsx_path))
  mp <- readxl::read_excel(xlsx_path, sheet = "DB_SHEET_MAP")
  stopifnot(all(c("Database","Sheet") %in% names(mp)))
  mp
}

# read a DB tab by DB name using the map (handles shortened tabs)
.read_enrichr_db <- function(xlsx_path, db_name, mp) {
  sh <- mp$Sheet[match(db_name, mp$Database)]
  if (is.na(sh)) return(NULL)
  df <- readxl::read_excel(xlsx_path, sheet = sh)
  if (!nrow(df)) return(NULL)
  if (!"Term" %in% names(df)) return(NULL)
  if (!"Combined.Score" %in% names(df) && "Combined score" %in% names(df)) {
    df$Combined.Score <- df[["Combined score"]]
  }
  df$Database <- db_name
  df
}

# top-hits plot for a set of DBs (colored by Database)
.plot_top_hits <- function(xlsx_path, dbs, out_png, title, top_n_total = 15) {
  mp <- .read_db_sheet_map(xlsx_path)
  
  # pull all terms from each DB, compute a unified Value, then combine
  dfs <- lapply(dbs, function(db) {
    df <- tryCatch(.read_enrichr_db(xlsx_path, db, mp), error = function(e) NULL)
    if (is.null(df) || !nrow(df)) return(NULL)
    
    # metric preference: Combined.Score else -log10(Adjusted.P.value)
    if ("Combined.Score" %in% names(df)) {
      df$Value  <- suppressWarnings(as.numeric(df$Combined.Score))
      df$Metric <- "Combined.Score"
    } else if ("Adjusted.P.value" %in% names(df)) {
      ap <- suppressWarnings(as.numeric(df$Adjusted.P.value))
      df$Value  <- ifelse(is.finite(ap) & ap > 0, -log10(ap), NA_real_)
      df$Metric <- "-log10(Adjusted.P.value)"
    } else {
      return(NULL)
    }
    
    df %>%
      dplyr::select(Term, Database, Value, Metric) %>%
      dplyr::filter(is.finite(Value), !is.na(Value), !is.na(Term), Term != "")
  })
  
  all_df <- dplyr::bind_rows(dfs)
  if (!nrow(all_df)) {
    message("[PLOT] No enrichment hits to plot for: ", basename(xlsx_path))
    return(invisible(NULL))
  }
  
  # ---- KEY CHANGE: rank globally and keep TOP N TOTAL ----
  top_df <- all_df %>%
    dplyr::arrange(dplyr::desc(Value)) %>%
    dplyr::slice_head(n = top_n_total)
  
  metric_label <- if (any(top_df$Metric == "Combined.Score")) "Combined.Score" else unique(top_df$Metric)[1]
  
  # order terms for plotting
  top_df <- top_df %>%
    dplyr::arrange(Value) %>%  # so biggest ends up at top after coord_flip
    dplyr::mutate(Term = factor(Term, levels = unique(Term)))
  
  p <- ggplot(top_df, aes(x = Term, y = Value, fill = Database)) +
    geom_col(color = "black", linewidth = 0.2) +
    coord_flip() +
    labs(title = title, x = NULL, y = metric_label, fill = "Database") +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text  = element_text(size = 10)
    )
  
  ggsave(out_png, p, dpi = 500, width = 12, height = 7, bg = "white")
  message("[PLOT] Saved: ", out_png)
  invisible(top_df)
}

plot_enrichr_for_pop <- function(pop_name, enrich_xlsx_path) {
  # TF plot
  .plot_top_hits(
    xlsx_path = enrich_xlsx_path,
    dbs       = tf_dbs,
    out_png   = file.path(out_tf, paste0("TopTFs_", safe_id(pop_name), ".png")),
    title     = paste0(pop_name, " | Top TF enrichments"),
    top_n     = 15
  )
  
  # Pathway plot
  .plot_top_hits(
    xlsx_path = enrich_xlsx_path,
    dbs       = pw_dbs,
    out_png   = file.path(out_pw, paste0("TopPathways_", safe_id(pop_name), ".png")),
    title     = paste0(pop_name, " | Top pathway enrichments"),
    top_n     = 15
  )
}

# Run for TotalNK + TotalCD8
enrich_xlsx_nk  <- file.path(out_step2, "TotalNK_OUD_Pos_vs_Neg_RNA_Enrichment.xlsx")
enrich_xlsx_cd8 <- file.path(out_step2, "TotalCD8_OUD_Pos_vs_Neg_RNA_Enrichment.xlsx")

if (file.exists(enrich_xlsx_nk))  plot_enrichr_for_pop("TotalNK", enrich_xlsx_nk)
if (file.exists(enrich_xlsx_cd8)) plot_enrichr_for_pop("TotalCD8", enrich_xlsx_cd8)

message("[STEP3] DONE. Plots in: ", out_step3)
