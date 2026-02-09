############################################################
# OPIS Grant Plots: Curated Volcano + Curated Pathways/TFs
# + Module Scoring 
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
opis_qs2  <- file.path(load.path, "OPIS_ALL_PostAnnotation.qs2")
OPIS_ALL  <- qs_read(opis_qs2)

# Input Excel
curated_xlsx <- file.path(load.path, "Curated genes for figures.xlsx")
stopifnot(file.exists(curated_xlsx))

# DGE directory (given)
dge_dir <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/Differential_Expression/DGE/OUD_Pos_vs_Neg"
stopifnot(dir.exists(dge_dir))

# Output base
out_base <- file.path(getwd(), "grant_plots")
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

out_volcano  <- file.path(out_base, "Volcano_CuratedGenes")
out_pathways <- file.path(out_base, "CuratedPathways_Plots")
out_tfs      <- file.path(out_base, "CuratedTFs_Plots")
out_modules  <- file.path(out_base, "Module_Scoring")

dir.create(out_volcano,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_pathways, recursive = TRUE, showWarnings = FALSE)
dir.create(out_tfs,      recursive = TRUE, showWarnings = FALSE)
dir.create(out_modules,  recursive = TRUE, showWarnings = FALSE)

# ---------------------------- #
# Helpers
# ---------------------------- #

.pick_fc_col <- function(df) {
  if ("avg_log2FC" %in% names(df)) return("avg_log2FC")
  if ("avg_logFC"  %in% names(df)) return("avg_logFC")
  # Some pipelines use "log2FC"
  if ("log2FC"     %in% names(df)) return("log2FC")
  stop("No fold-change column found (expected avg_log2FC, avg_logFC, or log2FC).")
}

.ensure_feature_labels <- function(df) {
  # If read.csv wrote rownames into first column named "X", move it to Feature
  if ("X" %in% names(df) && !"Feature" %in% names(df)) {
    df$Feature <- df$X
  }
  # If Feature still missing, try rownames
  if (!"Feature" %in% names(df)) {
    df$Feature <- rownames(df)
  }
  # Final fallback
  if (all(is.na(df$Feature)) || all(df$Feature == "")) {
    df$Feature <- paste0("feature_", seq_len(nrow(df)))
  }
  df
}

.read_dge <- function(csv_path) {
  df <- tryCatch(
    read.csv(csv_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE)
  )
  df <- .ensure_feature_labels(df)
  
  # Require at least p-values
  if (!all(c("p_val", "p_val_adj") %in% names(df))) {
    stop("Missing p_val and/or p_val_adj in: ", basename(csv_path))
  }
  
  fc_col <- .pick_fc_col(df)
  
  df <- df %>%
    filter(!is.na(.data[[fc_col]]), !is.na(p_val), !is.na(p_val_adj)) %>%
    filter(is.finite(.data[[fc_col]]), is.finite(p_val), is.finite(p_val_adj))
  
  list(df = df, fc_col = fc_col)
}

# Curated genes tab has NO header: one gene per row
.read_curated_genes_tab <- function(xlsx, sheet) {
  v <- readxl::read_excel(xlsx, sheet = sheet, col_names = FALSE) %>%
    dplyr::pull(1) %>%
    as.character() %>%
    str_trim() %>%
    .[. != "" & !is.na(.)]
  unique(v)
}

# Curated pathways/TF tabs: col1 = term, col2 = database, NO header
.read_curated_term_db <- function(xlsx, sheet) {
  df <- readxl::read_excel(xlsx, sheet = sheet, col_names = FALSE) %>%
    dplyr::select(term = 1, db = 2) %>%
    mutate(
      term = as.character(term) %>% str_trim(),
      db   = as.character(db) %>% str_trim()
    ) %>%
    filter(!is.na(term), term != "", !is.na(db), db != "")
  df
}

# Normalize cluster name into filename style you described:
# Example: "NK CD56dim" -> "NK_CD56dim"
.cluster_to_file_prefix <- function(cluster_name) {
  cluster_name %>%
    str_trim() %>%
    str_replace_all("\\s+", "_") %>%
    str_replace_all("[^A-Za-z0-9_]", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("^_|_$", "")
}

# Find enrichment result file(s) heuristically (so you don’t have to hardcode paths here).
# It searches within OPIS_ECCITEseq for CSV/TSV that include cluster + comparison + pathway/tf keywords.
.find_enrichment_file <- function(base_search_dir,
                                  cluster_prefix,
                                  kind = c("pathway","tf"),
                                  comparison_tag = "OUD_Pos_vs_Neg") {
  
  kind <- match.arg(kind)
  patt_kind <- if (kind == "pathway") "(pathway|enrich|enrichment|reactome|kegg|wikipathways|gobp|hallmark)"
  else "(tf|transcription|trrust|chea|jaspar|transfac)"
  
  # allow either underscores or spaces in underlying names
  cluster_regex <- str_replace_all(cluster_prefix, "_", "[ _-]?")
  
  # Candidate files
  files <- list.files(base_search_dir, recursive = TRUE, full.names = TRUE)
  files <- files[grepl("\\.(csv|tsv|txt)$", files, ignore.case = TRUE)]
  
  keep <- files[
    grepl(cluster_regex, files, ignore.case = TRUE) &
      grepl(comparison_tag, files, ignore.case = TRUE) &
      grepl(patt_kind, files, ignore.case = TRUE)
  ]
  
  # Prefer CSV first, then TSV
  keep <- keep[order(!grepl("\\.csv$", keep, ignore.case = TRUE))]
  
  if (!length(keep)) return(character(0))
  keep[1]  # return best single hit
}

# Read enrichment table (CSV/TSV) robustly
.read_enrichment_any <- function(path) {
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
}

# Choose plotting score column
.pick_score_cols <- function(df) {
  # Common Enrichr columns: "Adjusted.P.value", "Combined.Score"
  # or "adj_p", "combined_score", "p.adjust", etc.
  score_col <- NULL
  if ("Combined.Score" %in% names(df)) score_col <- "Combined.Score"
  if (is.null(score_col) && "combined_score" %in% names(df)) score_col <- "combined_score"
  if (is.null(score_col) && "Combined score" %in% names(df)) score_col <- "Combined score"
  
  adjp_col <- NULL
  if ("Adjusted.P.value" %in% names(df)) adjp_col <- "Adjusted.P.value"
  if (is.null(adjp_col) && "Adjusted P-value" %in% names(df)) adjp_col <- "Adjusted P-value"
  if (is.null(adjp_col) && "adj_p" %in% names(df)) adjp_col <- "adj_p"
  if (is.null(adjp_col) && "p_val_adj" %in% names(df)) adjp_col <- "p_val_adj"
  if (is.null(adjp_col) && "p.adjust" %in% names(df)) adjp_col <- "p.adjust"
  
  list(score_col = score_col, adjp_col = adjp_col)
}

# Match columns for term + db
.pick_term_db_cols <- function(df) {
  # Term column candidates
  term_candidates <- c("Term", "term", "Pathway", "pathway", "TF", "tf", "Transcription Factor", "transcription_factor")
  db_candidates   <- c("Gene_set", "gene_set", "Database", "database", "Library", "library", "db")
  
  term_col <- term_candidates[term_candidates %in% names(df)][1]
  db_col   <- db_candidates[db_candidates %in% names(df)][1]
  
  # If missing, try first column as term (common when exported weirdly)
  if (is.na(term_col) || is.null(term_col)) term_col <- names(df)[1]
  
  list(term_col = term_col, db_col = db_col)
}

# ---------------------------- #
# 1) Parse Excel sheets
# ---------------------------- #
sheets <- readxl::excel_sheets(curated_xlsx)

# Tabs naming rule you gave: "<cluster> curate(d) <something>"
# Be permissive about "curate" vs "curated" and casing.
is_dge_sheet <- grepl("curat(e|ed)\\s*DGE\\s*genes", sheets, ignore.case = TRUE)
is_pw_sheet  <- grepl("curat(e|ed)\\s*pathways",   sheets, ignore.case = TRUE)
is_tf_sheet  <- grepl("curat(e|ed)\\s*tfs?",       sheets, ignore.case = TRUE)

dge_sheets <- sheets[is_dge_sheet]
pw_sheets  <- sheets[is_pw_sheet]
tf_sheets  <- sheets[is_tf_sheet]

message("Found sheets:")
message("  DGE: ", paste(dge_sheets, collapse = " | "))
message("  PW : ", paste(pw_sheets,  collapse = " | "))
message("  TF : ", paste(tf_sheets,  collapse = " | "))

# Extract cluster name from sheet: everything before "curat"
.get_cluster_from_sheet <- function(sheet) {
  str_trim(str_split(sheet, regex("curat(e|ed)", ignore_case = TRUE))[[1]][1])
}

# ---------------------------- #
# 2) Volcano plots for curated genes
# ---------------------------- #
for (sh in dge_sheets) {
  cluster_name   <- .get_cluster_from_sheet(sh)
  cluster_prefix <- .cluster_to_file_prefix(cluster_name)
  
  dge_file <- file.path(dge_dir, paste0(cluster_prefix, "_OUD_Pos_vs_Neg_RNA.csv"))
  if (!file.exists(dge_file)) {
    message("[VOLCANO] Missing DGE file for ", cluster_name, " -> expected: ", dge_file)
    next
  }
  
  curated_genes <- .read_curated_genes_tab(curated_xlsx, sh)
  
  # Load DGE
  d <- tryCatch(.read_dge(dge_file), error = function(e) {
    message("[VOLCANO] Failed reading DGE for ", cluster_name, ": ", conditionMessage(e))
    return(NULL)
  })
  if (is.null(d)) next
  
  df     <- d$df
  fc_col <- d$fc_col
  
  # keep only curated genes that exist in DGE
  curated_present <- intersect(curated_genes, df$Feature)
  if (!length(curated_present)) {
    message("[VOLCANO] No curated genes found in DGE for ", cluster_name, " (sheet: ", sh, ")")
    next
  }
  
  # Label only curated genes; leave others unlabeled
  lab_vec <- ifelse(df$Feature %in% curated_present, df$Feature, "")
  
  # Build volcano
  vp <- EnhancedVolcano(
    df,
    lab = lab_vec,
    x   = fc_col,
    y   = "p_val_adj",              # use adjusted p for axis thresholding/labels
    xlab = bquote(~Log[2]~'fold change'),
    ylab = bquote(~-Log[10]~'adj p-value'),
    pCutoff  = 0.01,
    FCcutoff = 0.25,                # keep permissive since curated may include modest FC
    pointSize = 2.2,
    labSize   = 4.5,
    labCol    = "black",
    labFace   = "bold",
    colAlpha  = 0.85,
    legendPosition = "right",
    boxedLabels=T,
    legendLabSize  = 12,
    legendIconSize = 3.8,
    drawConnectors = TRUE,
    widthConnectors = 0.5,
    colConnectors   = "black",
    max.overlaps    = 100,
    title    = paste0(cluster_name, " | OUD Pos vs Neg (RNA)"),
    subtitle = "Curated genes labeled"
  )
  
  out_png <- file.path(out_volcano, paste0("Volcano_", cluster_prefix, "_CuratedGenes.png"))
  ggsave(out_png, vp, dpi = 500, width = 10, height = 7, bg = "white")
  message("[VOLCANO] Saved: ", out_png)
}

# ---------------------------- #
# 3) Curated pathway plots (DETERMINISTIC EnrichR Excel)
# ---------------------------- #

enrich_dir <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/Pathway_Analysis_EnrichR/RNA_OUD_Pos_vs_Neg/CSVs"
stopifnot(dir.exists(enrich_dir))

# helper: read one DB tab and keep required columns
.read_enrichr_tab <- function(xlsx_path, sheet) {
  df <- readxl::read_excel(xlsx_path, sheet = sheet)
  # Require at least Term
  if (!"Term" %in% names(df)) stop("No 'Term' column in: ", basename(xlsx_path), " | tab: ", sheet)
  
  # Normalize expected columns if present
  if (!"Combined.Score" %in% names(df) && "Combined score" %in% names(df)) {
    df$Combined.Score <- df[["Combined score"]]
  }
  df
}
.match_enrichr_tab <- function(db_i, tabs) {
  db_raw <- stringr::str_trim(as.character(db_i))
  if (is.na(db_raw) || db_raw == "") return(NA_character_)
  
  tabs_trim <- stringr::str_trim(tabs)
  
  # helper: normalize by removing ALL non-alphanumerics
  norm <- function(x) gsub("[^A-Za-z0-9]+", "", tolower(stringr::str_trim(x)))
  
  db_norm   <- norm(db_raw)
  tabs_norm <- norm(tabs_trim)
  
  # 1) exact normalized match
  idx <- which(tabs_norm == db_norm)
  if (length(idx)) return(tabs_trim[idx[1]])
  
  # 2) db is a prefix of tab (handles truncated years like ...20 vs ...2019)
  idx <- which(startsWith(tabs_norm, db_norm))
  if (length(idx)) return(tabs_trim[idx[1]])
  
  # 3) tab is a prefix of db (rare but safe)
  idx <- which(startsWith(db_norm, tabs_norm))
  if (length(idx)) return(tabs_trim[idx[1]])
  
  # 4) contains match either direction
  idx <- which(grepl(db_norm, tabs_norm, fixed = TRUE))
  if (length(idx)) return(tabs_trim[idx[1]])
  
  NA_character_
}

.plot_curated_enrich_terms <- function(curated_df, xlsx_path, cluster_name, out_png, plot_title) {
  # curated_df has columns term, db
  # Each db corresponds to a sheet/tab in xlsx_path
  
  if (!file.exists(xlsx_path)) {
    message("[ENRICH] Missing enrichment xlsx: ", xlsx_path)
    return(invisible(NULL))
  }
  
  tabs <- readxl::excel_sheets(xlsx_path)
  
  hits <- lapply(seq_len(nrow(curated_df)), function(i) {
    term_i <- curated_df$term[i]
    db_i   <- curated_df$db[i]
    
    tab <- .match_enrichr_tab(db_i, tabs)
    if (is.na(tab)) {
      message("[ENRICH] No tab match for database: ", db_i, " | available tabs: ", paste(tabs, collapse = " | "))
      return(NULL)
    }
    
    df  <- tryCatch(.read_enrichr_tab(xlsx_path, tab), error = function(e) NULL)
    if (is.null(df)) return(NULL)
    
    df$Term_l <- str_to_lower(str_trim(df$Term))
    term_l    <- str_to_lower(str_trim(term_i))
    
    row <- df[df$Term_l == term_l, , drop = FALSE]
    if (!nrow(row)) return(NULL)
    
    score  <- NA_real_
    metric <- NA_character_
    
    if ("Combined.Score" %in% names(row)) {
      score <- suppressWarnings(as.numeric(row$Combined.Score[1]))
      metric <- "Combined.Score"
    }
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
    
    data.frame(
      Term = term_i,
      Database = db_i,
      Value = score,
      Metric = metric,
      stringsAsFactors = FALSE
    )
  })
  
  hit_df <- bind_rows(hits)
  
  if (!nrow(hit_df)) {
    message("[ENRICH] No curated terms matched for ", cluster_name, " | ", basename(xlsx_path))
    return(invisible(NULL))
  }
  
  metric_label <- if (any(hit_df$Metric == "Combined.Score")) "Combined.Score" else unique(hit_df$Metric)[1]
  
  # Order by score (descending), keep factor order for plotting
  hit_df <- hit_df %>%
    arrange(desc(Value)) %>%
    mutate(
      Term = factor(Term, levels = rev(unique(Term))),
      Database = factor(Database, levels = sort(unique(Database)))
    )
  
  p <- ggplot(hit_df, aes(x = Term, y = Value, fill = Database)) +
    geom_col(color = "black", linewidth = 0.2) +
    coord_flip() +
    labs(
      title = plot_title,
      x = NULL,
      y = metric_label,
      fill = "Database"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.y = element_text(size = 13),   # <- THIS controls Term labels
      axis.text.x = element_text(size = 11),   # <- numeric axis (Value)
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text  = element_text(size = 10)
    )
  
  
  ggsave(
    out_png,
    p,
    dpi = 500,
    width = 11,
    height = max(5, 0.25 * nrow(hit_df) + 2),
    bg = "white"
  )
  message("[ENRICH] Saved: ", out_png)
  
  invisible(hit_df)
}


# Pathways sheets
for (sh in pw_sheets) {
  cluster_name   <- .get_cluster_from_sheet(sh)
  cluster_prefix <- .cluster_to_file_prefix(cluster_name)
  
  curated_pw <- .read_curated_term_db(curated_xlsx, sh)
  
  enrich_xlsx <- file.path(enrich_dir, paste0(cluster_prefix, "_OUD_Pos_vs_Neg_RNA_Enrichment.xlsx"))
  
  out_png <- file.path(out_pathways, paste0("CuratedPathways_", cluster_prefix, ".png"))
  .plot_curated_enrich_terms(
    curated_df = curated_pw,
    xlsx_path  = enrich_xlsx,
    cluster_name = cluster_name,
    out_png    = out_png,
    plot_title = paste0(cluster_name, " | Curated pathways")
  )
}
# ============================================================
# One-off replot for NK_CD56dim (exclude 2 pathways)
# ============================================================

cluster_name   <- "NK CD56dim"
cluster_prefix <- "NK_CD56dim"

# Pathways to remove (exact, case-insensitive)
exclude_terms <- c(
  "Natural Killer cell-mediated cytotoxicity",
  "Human immunodeficiency virus 1 infection",
  "Immune system signaling by interferons, interleukins, prolactin, and growth hormones",
  "Interactions Of Natural Killer Cells In Pancreatic Cancer WP5092"
)

# Read curated pathways for NK
curated_pw <- .read_curated_term_db(
  curated_xlsx,
  sheet = "NK CD56dim curated pathways"
)

# Drop excluded pathways
curated_pw <- curated_pw %>%
  mutate(term_norm = str_to_lower(str_trim(term))) %>%
  filter(!term_norm %in% str_to_lower(str_trim(exclude_terms))) %>%
  select(-term_norm)

# EnrichR Excel (already exists)
enrich_xlsx <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/Pathway_Analysis_EnrichR/RNA_OUD_Pos_vs_Neg/CSVs/NK_CD56dim_OUD_Pos_vs_Neg_RNA_Enrichment.xlsx"

# Output (overwrite or new file name if you want)
out_png <- file.path(
  out_pathways,
  "CuratedPathways_NK_CD56dim.png"
)

# Plot
.plot_curated_enrich_terms(
  curated_df   = curated_pw,
  xlsx_path    = enrich_xlsx,
  cluster_name = cluster_name,
  out_png      = out_png,
  plot_title   = "NK CD56dim | Curated pathways"
)
# ---------------------------- #
# 4) Curated TF plots (DETERMINISTIC EnrichR Excel)
# ---------------------------- #
for (sh in tf_sheets) {
  cluster_name   <- .get_cluster_from_sheet(sh)
  cluster_prefix <- .cluster_to_file_prefix(cluster_name)
  
  curated_tf <- .read_curated_term_db(curated_xlsx, sh)
  
  enrich_xlsx <- file.path(enrich_dir, paste0(cluster_prefix, "_OUD_Pos_vs_Neg_RNA_Enrichment.xlsx"))
  
  out_png <- file.path(out_tfs, paste0("CuratedTFs_", cluster_prefix, ".png"))
  .plot_curated_enrich_terms(
    curated_df = curated_tf,
    xlsx_path  = enrich_xlsx,
    cluster_name = cluster_name,
    out_png    = out_png,
    plot_title = paste0(cluster_name, " | Curated transcription factors")
  )
}

# ---------------------------- #
# 5) Module scoring 
# ---------------------------- #

##### Functions

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

pretty_module_name <- function(mod_name) {
  x <- mod_name
  x <- gsub("_", " ", x, fixed = TRUE)
  x
}

# Score ONE module into meta.data field "MS_<mod_id>"
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

######## Plotters
plot_module_feature_scCustom <- function(obj, module_field, folder,
                                         reduction = "wnn.umap",
                                         title_prefix = "RNA",
                                         pal_n = 10, pal_option = "A") {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  pal <- viridis::viridis(n = pal_n, option = pal_option)
  
  p <- scCustomize::FeaturePlot_scCustom(
    obj,
    reduction  = reduction,
    features   = module_field,
    colors_use = pal,
    order      = TRUE
  ) + ggtitle(paste(title_prefix, "|", module_field, "|", reduction))
  
  ggsave(file.path(folder, paste0("Feature_", module_field, ".png")),
         p, dpi = 500, width = 10, height = 8, bg = "white")
  invisible(p)
}

plot_module_feature_split <- function(obj, module_field, folder,
                                      split_by = "OUD_status",
                                      reduction = "wnn.umap",
                                      title_prefix = "RNA",
                                      pal_n = 10, pal_option = "A") {
  if (!split_by %in% colnames(obj@meta.data)) return(invisible(NULL))
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  pal <- viridis::viridis(n = pal_n, option = pal_option)
  
  p <- scCustomize::FeaturePlot_scCustom(
    obj,
    reduction = reduction,
    features  = module_field,
    split.by  = split_by,
    order     = TRUE,
    cols      = pal
  ) 
  
  ggsave(file.path(folder, paste0("Feature_", module_field, "_SplitBy_", split_by, ".png")),
         p, dpi = 500, width = 14, height = 7, bg = "white")
  invisible(p)
}

plot_module_vln_seuratExtend <- function(obj, module_field, folder,
                                         split_by = "OUD_status",
                                         title_prefix = "RNA") {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  
  if (split_by %in% colnames(obj@meta.data)) {
    p <- SeuratExtend::VlnPlot2(
      obj,
      features     = module_field,
      cols         = "default",
      split.by     = split_by,
      stat.method  = "wilcox.test",
      show.mean    = TRUE,
      mean_colors  = c("red", "blue")
    ) + ggtitle(paste(title_prefix, "|", module_field, "| Split by", split_by))
  } else {
    p <- SeuratExtend::VlnPlot2(
      obj,
      features     = module_field,
      cols         = "default",
      stat.method  = "wilcox.test",
      show.mean    = TRUE,
      mean_colors  = c("red", "blue")
    ) + ggtitle(paste(title_prefix, "|", module_field))
  }
  
  ggsave(file.path(folder, paste0("Vln_", module_field, ".png")),
         p, dpi = 500, width = 12, height = 7, bg = "white")
  invisible(p)
}

plot_module_box_by_cluster <- function(obj, module_field, folder,
                                       group_col   = "OUD_status",
                                       cluster_col = "Manual_Annotation",
                                       title = NULL,
                                       group_labels = NULL,
                                       group_colors = NULL) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  
  md <- obj@meta.data
  if (!all(c(group_col, cluster_col, module_field) %in% colnames(md))) {
    message("[BOX] Missing columns for ", module_field, " | Need: ",
            paste(c(group_col, cluster_col, module_field), collapse = ", "))
    return(invisible(NULL))
  }
  
  keep <- !is.na(md[[group_col]]) & !is.na(md[[cluster_col]]) & !is.na(md[[module_field]])
  md <- md[keep, , drop = FALSE]
  if (!nrow(md)) return(invisible(NULL))
  
  df <- data.frame(
    .group   = as.character(md[[group_col]]),
    .cluster = as.character(md[[cluster_col]]),
    .score   = as.numeric(md[[module_field]]),
    stringsAsFactors = FALSE
  )
  
  # -----------------------------
  # Robust ordering + default labels/colors
  # -----------------------------
  lvl <- unique(df$.group)
  
  # Prefer common 2-group orderings
  if (all(c("OUD-", "OUD+") %in% lvl)) {
    df$.group <- factor(df$.group, levels = c("OUD-", "OUD+"))
  } else if (all(c("Negative", "Positive") %in% lvl)) {
    df$.group <- factor(df$.group, levels = c("Negative", "Positive"))
  } else if (length(lvl) >= 2) {
    # deterministic fallback: alphabetical
    df$.group <- factor(df$.group, levels = sort(unique(df$.group)))
  } else {
    df$.group <- factor(df$.group)
  }
  
  # Defaults for OUD-/OUD+ if present
  default_labels_oud  <- c("OUD-" = "OUD-", "OUD+" = "OUD+")
  default_colors_oud  <- c("OUD-" = "#1180808B", "OUD+" = "#F54927")
  
  # Defaults for Negative/Positive if present
  default_labels_np <- c("Negative" = "Negative", "Positive" = "Positive")
  default_colors_np <- c("Negative" = "#1180808B", "Positive" = "#F54927")
  
  levs <- levels(df$.group)
  
  # If user didn't supply labels/colors, auto-fill for known schemes
  if (is.null(group_labels)) {
    if (all(names(default_labels_oud) %in% levs)) {
      group_labels <- default_labels_oud
    } else if (all(names(default_labels_np) %in% levs)) {
      group_labels <- default_labels_np
    }
  }
  
  if (is.null(group_colors)) {
    if (all(names(default_colors_oud) %in% levs)) {
      group_colors <- default_colors_oud
    } else if (all(names(default_colors_np) %in% levs)) {
      group_colors <- default_colors_np
    }
  }
  
  # If user supplied colors but names don't match levels, try to subset safely
  if (!is.null(group_colors)) {
    group_colors <- group_colors[names(group_colors) %in% levs]
    if (length(group_colors) == 0) group_colors <- NULL
  }
  if (!is.null(group_labels)) {
    group_labels <- group_labels[names(group_labels) %in% levs]
    if (length(group_labels) == 0) group_labels <- NULL
  }
  
  if (is.null(title)) title <- module_field
  pretty_title <- pretty_module_name(title)
  file_stub    <- safe_filename(module_field)
  
  p <- ggplot(df, aes(x = .group, y = .score, fill = .group)) +
    geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0.9, color = "black") +
    geom_jitter(width = 0.18, size = 0.4, alpha = 0.25) +
    ggpubr::stat_compare_means(
      method = "wilcox.test",
      label = "p.signif",
      hide.ns = TRUE,
      label.x.npc = "center",
      label.y.npc = 0.9,
      size = 7,
      color = "red"
    ) +
    facet_wrap(~ .cluster, scales = "free_y") +
    labs(title = pretty_title, x = group_col, y = "Module Score", fill = group_col) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 16) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(size = 14, face = "bold"),
      axis.text = element_text(color = "black"),
      legend.position = "top",
      legend.title = element_text(face = "bold")
    )
  
  # Apply labels (x-axis + legend) if provided
  if (!is.null(group_labels)) {
    p <- p + scale_x_discrete(labels = group_labels)
  }
  
  # Apply fills if provided; otherwise ggplot default palette will be used
  if (!is.null(group_colors)) {
    p <- p + scale_fill_manual(values = group_colors, breaks = names(group_colors))
  }
  
  ggsave(
    file.path(folder, paste0("Box_", file_stub, "_ByCluster.png")),
    p, width = 20, height = 18, dpi = 300, bg = "white"
  )
  
  invisible(p)
}

# ---------------------------- #
# 5) Run scoring + plots for each module
# ---------------------------- #
assay_use   <- "RNA"
reduction_use <- "wnn.umap"
group_col  <- "OUD_status"
cluster_col <- "seurat_clusters"  # change to your preferred cluster label if needed

has_oud <- group_col %in% colnames(OPIS_ALL@meta.data)

# optional: if your OUD labels are literally "Positive"/"Negative"
group_labels <- c("OUD-" = "OUD-", "OUD+" = "OUD+")
group_colors <- c("OUD-" = "#1180808B", "OUD+" = "#F54927")


# ---------------------------- #
# A) SCORING ONLY
# ---------------------------- #

coverage <- data.frame()

for (mod_name in names(modules)) {
  genes  <- modules[[mod_name]]
  mod_id <- safe_id(mod_name)
  
  # score
  OPIS_ALL <- score_module(OPIS_ALL, genes, mod_id, assay_use = assay_use)
  module_field <- paste0("MS_", mod_id)
  
  # coverage reporting
  present <- intersect(genes, rownames(OPIS_ALL[[assay_use]]))
  missing <- setdiff(genes, present)
  
  coverage <- rbind(
    coverage,
    data.frame(
      Module        = mod_name,
      Meta_Field    = module_field,
      Present       = length(present),
      Missing       = length(missing),
      Present_Genes = paste(present, collapse = ";"),
      Missing_Genes = paste(missing, collapse = ";"),
      stringsAsFactors = FALSE
    )
  )
}

# Save coverage + scored object (optional but recommended)
write.csv(coverage, file.path(out_modules, "ModuleScore_Coverage.csv"), row.names = FALSE)
qs2::qs_save(OPIS_ALL, file.path(out_modules, "OPIS_ALL_PostAnnotation_IFNModules.qs2"))


# ---------------------------- #
# B) PLOTTING ONLY
# ---------------------------- #
clustercells.high_col = "Manual_Annotation"
has_oud <- group_col %in% colnames(OPIS_ALL@meta.data)

for (mod_name in names(modules)) {
  mod_id <- safe_id(mod_name)
  module_field <- paste0("MS_", mod_id)
  
  # Skip if scoring missing (protects you if a module had 0 present genes)
  if (!module_field %in% colnames(OPIS_ALL@meta.data)) {
    message("[PLOT] Missing module field (skipping): ", module_field)
    next
  }
  
  # per-module folder
  mod_dir <- file.path(out_modules, mod_id)
  dir.create(mod_dir, recursive = TRUE, showWarnings = FALSE)
  
  # plots
  plot_module_feature_scCustom(
    OPIS_ALL, module_field, mod_dir,
    reduction = reduction_use, title_prefix = "RNA"
  )
  
  if (has_oud) {
    plot_module_feature_split(
      OPIS_ALL, module_field, mod_dir,
      split_by = group_col, reduction = reduction_use, title_prefix = "RNA"
    )
  }
  
  plot_module_vln_seuratExtend(
    OPIS_ALL, module_field, mod_dir,
    split_by = group_col, title_prefix = "RNA"
  )
  
  # Boxplots only make sense if you have a 2-group column; otherwise still draws but no stats
  if (has_oud) {
    plot_module_box_by_cluster(
      OPIS_ALL, module_field, mod_dir,
      group_col = group_col, cluster_col = cluster_col,
      title = mod_name,
      group_labels = group_labels,
      group_colors = group_colors
    )
  } else {
    plot_module_box_by_cluster(
      OPIS_ALL, module_field, mod_dir,
      group_col = group_col, cluster_col = cluster_col,
      title = mod_name
    )
  }
}  

##########3


Cluster_Highlight_Plot(seurat_object = OPIS_ALL, 
                       reduction  = "wnn.umap",
                       cluster_name = c("CD8+ TEMRA", "NK CD56dim"), 
                       highlight_color = c("#E0B0FF", '#9FE2BF'),
                       background_color = 'black',
                       alpha=4
                       )
ggsave(
  filename = file.path(out_base, "temra_dim_labled_umap.png"),
  width    = 8,
  height   = 7,
  dpi      = 400,
  bg='white'
)
