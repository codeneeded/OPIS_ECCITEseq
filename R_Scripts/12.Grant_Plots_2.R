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
