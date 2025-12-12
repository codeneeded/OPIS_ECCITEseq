library(dplyr)
library(ggplot2)
library(EnhancedVolcano)

# --------------------------------------------
# Base paths (OPIS)
# --------------------------------------------
de_base <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/Differential_Expression"

dge_dir <- file.path(de_base, "DGE", "OUD_Pos_vs_Neg")
dpe_dir <- file.path(de_base, "DPE", "OUD_Pos_vs_Neg")

plot_base <- file.path(de_base, "Volcano_Plots")
plot_dir_rna <- file.path(plot_base, "RNA_OUD_Pos_vs_Neg")
plot_dir_adt <- file.path(plot_base, "ADT_OUD_Pos_vs_Neg")

dir.create(plot_dir_rna, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir_adt, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------
# Helpers
# --------------------------------------------
.pick_fc_col <- function(df) {
  if ("avg_log2FC" %in% names(df)) return("avg_log2FC")
  if ("avg_logFC"  %in% names(df)) return("avg_logFC")
  stop("No fold-change column found (expected avg_log2FC or avg_logFC).")
}

.ensure_feature_labels <- function(df) {
  # If read.csv wrote rownames into first column named "X", move it to Feature
  if ("X" %in% names(df) && !"Feature" %in% names(df)) {
    df$Feature <- df$X
  }
  # If Feature still missing, try rownames (won't exist if you didn't set row.names=1)
  if (!"Feature" %in% names(df)) {
    df$Feature <- rownames(df)
  }
  # Final fallback: create something
  if (all(is.na(df$Feature)) || all(df$Feature == "")) {
    df$Feature <- paste0("feature_", seq_len(nrow(df)))
  }
  df
}

.make_keyvals <- function(df, fc_col, fc_cut = 1.5, padj_cut = 0.01) {
  keyvals <- ifelse(
    abs(df[[fc_col]]) > fc_cut & df$p_val_adj < padj_cut, "#CD0BBC",
    ifelse(df$p_val_adj < padj_cut, "#28E2E5", "gray30")
  )
  keyvals[is.na(keyvals)] <- "gray30"
  names(keyvals)[keyvals == "gray30"]   <- "NS"
  names(keyvals)[keyvals == "#28E2E5"]  <- paste0("adj(p-value) < ", padj_cut)
  names(keyvals)[keyvals == "#CD0BBC"]  <- paste0("FC > ", fc_cut)
  keyvals
}

.plot_volcano_dir <- function(input_dir, out_dir, comparison_title,
                              fc_cut = 1.5, padj_cut = 0.01,
                              subtitle = "Differential Expression") {
  files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
  if (!length(files)) {
    message("No CSVs found in: ", input_dir)
    return(invisible(NULL))
  }
  
  for (f in files) {
    message("Volcano: ", basename(f))
    
    # IMPORTANT: many Seurat outputs were written with rownames; use row.names=1
    df <- tryCatch(
      read.csv(f, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE),
      error = function(e) {
        message("  Read failed: ", basename(f), " — trying without row.names=1")
        read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
      }
    )
    
    if (!all(c("p_val", "p_val_adj") %in% names(df))) {
      message("  Skipping (need p_val and p_val_adj): ", basename(f),
              " | cols: ", paste(names(df), collapse = ", "))
      next
    }
    
    df <- .ensure_feature_labels(df)
    fc_col <- .pick_fc_col(df)
    
    # drop missing
    df <- df %>%
      filter(!is.na(.data[[fc_col]]), !is.na(p_val), !is.na(p_val_adj)) %>%
      filter(is.finite(.data[[fc_col]]), is.finite(p_val), is.finite(p_val_adj))
    
    if (!nrow(df)) {
      message("  No rows after filtering: ", basename(f))
      next
    }
    
    keyvals <- .make_keyvals(df, fc_col = fc_col, fc_cut = fc_cut, padj_cut = padj_cut)
    
    # Cluster name = filename without .csv
    clust_name <- tools::file_path_sans_ext(basename(f))
    
    vp <- EnhancedVolcano(
      df,
      lab = df$Feature,
      x = fc_col,
      y = "p_val",
      xlab = bquote(~Log[2]~ 'fold change'),
      pCutoff  = padj_cut,
      FCcutoff = fc_cut,
      pointSize = 3.5,
      labSize   = 3.0,
      labCol    = "black",
      labFace   = "bold",
      colAlpha  = 4/5,
      legendPosition = "right",
      legendLabSize  = 14,
      legendIconSize = 4.0,
      drawConnectors = TRUE,
      widthConnectors = 0.8,
      colConnectors   = "black",
      title    = paste0(clust_name, ": ", comparison_title),
      subtitle = subtitle,
      colCustom = keyvals
    )
    
    ggsave(
      filename = file.path(out_dir, paste0("Volcano_", clust_name, ".png")),
      plot     = vp + guides(color = guide_legend(reverse = TRUE)),
      dpi      = 500,
      width    = 10,
      height   = 7,
      bg       = "white"
    )
  }
  
  invisible(NULL)
}

# --------------------------------------------
# Run volcano plots
# --------------------------------------------

.plot_volcano_dir(
  input_dir = dge_dir,
  out_dir   = plot_dir_rna,
  comparison_title = "OUD+ vs OUD-",
  fc_cut = 1.5,
  padj_cut = 0.01,
  subtitle = "Differential Gene Expression (RNA)"
)

.plot_volcano_dir(
  input_dir = dpe_dir,
  out_dir   = plot_dir_adt,
  comparison_title = "OUD+ vs OUD-",
  fc_cut = 1.5,
  padj_cut = 0.01,
  subtitle = "Differential Protein Expression (ADT)"
)
#############For Grant

# ----------------------------
# Paths
# ----------------------------
base_dir <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq"
dge_dir  <- file.path(base_dir, "Differential_Expression", "DGE", "OUD_Pos_vs_Neg")

grant_dir <- file.path(
  base_dir,
  "Differential_Expression",
  "Volcano_Plots",
  "For_Grant"
)
dir.create(grant_dir, recursive = TRUE, showWarnings = FALSE)

infile <- file.path(dge_dir, "Classical_CD14_mono_OUD_Pos_vs_Neg_RNA.csv")

# ----------------------------
# Genes to label (all up in OUD+)
# ----------------------------
highlight_genes <- c("ALDOA","PGK1","GAPDH","ACTG1","MT-ND5","GPX1","NDUFB9")

# ----------------------------
# Read DGE results
# ----------------------------
deg <- read.csv(infile, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

# FC col can vary by Seurat version
fc_col <- if ("avg_log2FC" %in% names(deg)) "avg_log2FC" else
  if ("avg_logFC" %in% names(deg)) "avg_logFC" else
    if ("Av_Log2FC" %in% names(deg)) "Av_Log2FC" else
      stop("No fold-change column found (expected avg_log2FC / avg_logFC / Av_Log2FC).")

deg$Gene <- rownames(deg)

# Keep finite rows
deg <- deg %>%
  filter(!is.na(.data[[fc_col]]), !is.na(p_val_adj)) %>%
  filter(is.finite(.data[[fc_col]]), is.finite(p_val_adj))

padj_cut <- 0.05

# ----------------------------
# Two-color scheme only:
#   - significant vs not significant
# ----------------------------
keyvals <- ifelse(deg$p_val_adj < padj_cut, "#28E2E5", "gray30")
keyvals[is.na(keyvals)] <- "gray30"
names(keyvals)[keyvals == "#28E2E5"] <- paste0("adj(p) < ", padj_cut)
names(keyvals)[keyvals == "gray30"]  <- "NS"

# ----------------------------
# Volcano plot
#   - No FC cutoff
#   - Label only the specified genes
# ----------------------------
vp <- EnhancedVolcano(
  deg,
  lab = deg$Gene,
  x = fc_col,
  y = "p_val_adj",
  xlab = bquote(~Log[2]~ "fold change (OUD+ vs OUD-)"),
  ylab = bquote(~-Log[10]~ "adj. p-value"),
  pCutoff  = padj_cut,
  FCcutoff = 0,                 # no fold-change cutoff
  pointSize = 3.2,
  labSize   = 4.0,
  labCol    = "black",
  labFace   = "bold",
  colAlpha  = 4/5,
  colCustom = keyvals,          # ONLY 2 colors used
  selectLab = highlight_genes,  # label these only
  drawConnectors = TRUE,
  widthConnectors = 0.8,
  colConnectors   = "black",
  boxedLabels = TRUE,
  legendPosition = "right",
  legendLabSize  = 13,
  legendIconSize = 4.0,
  title    = "Classical CD14+ Monocytes: OUD+ vs OUD-",
  subtitle = "Selected genes labeled (no logFC cutoff); significance by adj. p-value"
)

# Save PNG only
out_png <- file.path(grant_dir, "Volcano_Classical_CD14_mono_OUD_Pos_vs_Neg_RNA.png")

ggsave(
  filename = out_png,
  plot     = vp + guides(color = guide_legend(reverse = TRUE)),
  dpi      = 600,
  width    = 10.5,
  height   = 7.5,
  bg       = "white"
)

message("Saved: ", out_png)
