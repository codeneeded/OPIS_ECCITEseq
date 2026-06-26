# =============================================================================
# OPIS ECCITEseq — Monocytes: annotate + DGE (RNA) + DPE (ADT) + EnrichR
#
# 1. Apply the final res-0.4 annotation (Monocyte_annotation_README), drop the
#    T-doublet cluster (9), relabel the REVIEW cluster (6) to an analysable name.
# 2. DGE (RNA, MAST):  (a) population identity markers (one-vs-rest)
#                      (b) OUD+ vs OUD- within each population (+ all monocytes)
# 3. DPE (ADT, wilcox on DSB data): same two contrasts.
# 4. EnrichR pathway analysis on the RNA DGE (identity = all sig genes/pop;
#    OUD = up and down separately, since direction is meaningful).
#
# Significant table format matches the project standard: ONE file per analysis
# with both directions (p_val_adj < 0.05 & abs(avg_log2FC) > threshold).
#
# Reads : OPIS_MONO_res0.4.qs2     Writes: Monocytes/DGE_DPE/ , Monocytes/Pathway_EnrichR/
# =============================================================================

# ---- Toggles -----------------------------------------------------------------
RUN_IDENTITY_DGE <- TRUE     # RNA one-vs-rest population markers (MAST)
RUN_IDENTITY_DPE <- TRUE     # ADT one-vs-rest population markers (wilcox)
RUN_OUD_DGE      <- TRUE     # RNA OUD+ vs OUD- per population (MAST)
RUN_OUD_DPE      <- TRUE     # ADT OUD+ vs OUD- per population (wilcox)
RUN_ABUNDANCE    <- TRUE     # ClusterDistrBar cluster-distribution / abundance plots
RUN_ENRICHR      <- TRUE     # EnrichR on the RNA DGE (needs internet + enrichR + openxlsx)

# ---- Config ------------------------------------------------------------------
group_var      <- "Mono_label"   # analysis grouping; switch to "Mono_fine" for sub-states
DGE_TEST       <- "MAST"         # RNA test (project standard)
DGE_LATENT     <- NULL           # e.g. "orig.ident" to adjust for donor (see caveat at end)
dge_lfc        <- 0.25           # |avg_log2FC| significance cut, RNA
dpe_lfc        <- 0.25           # |avg_log2FC| significance cut, ADT (lower if too strict)
oud_min_cells  <- 20             # min cells per OUD group within a population to test

suppressPackageStartupMessages({
  library(Seurat); library(SeuratExtend); library(qs2)
  library(dplyr); library(tidyr); library(tibble); library(ggplot2)
})

# ---- Paths -------------------------------------------------------------------
load.path <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/saved_R_data"
mono.base <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/subclustering/pre_annotation/Monocytes"
dge.root  <- file.path(mono.base, "DGE_DPE")
enr.root  <- file.path(mono.base, "Pathway_EnrichR")
d_idg <- file.path(dge.root, "Identity_DGE_RNA"); d_idp <- file.path(dge.root, "Identity_DPE_ADT")
d_oug <- file.path(dge.root, "OUD_DGE_RNA");      d_oup <- file.path(dge.root, "OUD_DPE_ADT")
abund.dir <- file.path(mono.base, "abundance")
for (d in c(d_idg, d_idp, d_oug, d_oup, abund.dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. LOAD + APPLY ANNOTATION
# =============================================================================
message("Loading OPIS_MONO_res0.4.qs2 ...")
OPIS_MONO <- qs_read(file.path(load.path, "OPIS_MONO_res0.4.qs2"))
if (inherits(OPIS_MONO[["RNA"]], "Assay5"))
  OPIS_MONO[["RNA"]] <- tryCatch(JoinLayers(OPIS_MONO[["RNA"]]), error = function(e) OPIS_MONO[["RNA"]])

mono.fine <- c(
  "0"="Classical monocyte (HSP/stress)", "1"="Classical monocyte (VCAN+/platelet-assoc)",
  "2"="Inflammatory monocyte", "3"="Classical monocyte (hypoxic/glycolytic)",
  "4"="MHC-II-high antigen-presenting monocyte", "5"="CD16+ monocyte (ISG-high)",
  "6"="Classical monocyte (CCR2/IRF8-high) [REVIEW]", "7"="HLA-G+ monocyte",
  "8"="Non-classical monocyte (CD16++)", "9"="EXCLUDE - T-cell doublet/contamination")
mono.group <- c(
  "0"="Classical monocyte","1"="Classical monocyte","3"="Classical monocyte",
  "2"="Inflammatory monocyte","4"="Antigen-presenting monocyte","7"="HLA-G+ monocyte",
  "5"="Non-classical/CD16+ monocyte","8"="Non-classical/CD16+ monocyte",
  "6"="REVIEW","9"="EXCLUDE")

sc <- as.character(OPIS_MONO$Subcluster_ID)
OPIS_MONO$Mono_fine  <- factor(unname(mono.fine [sc]))
OPIS_MONO$Mono_group <- factor(unname(mono.group[sc]))

# Drop the T-cell doublet cluster (9 / EXCLUDE); keep cluster 6 but give the
# REVIEW group an analysable label.
OPIS_MONO <- subset(OPIS_MONO, subset = Subcluster_ID != "9")
lab <- as.character(OPIS_MONO$Mono_group)
lab[lab == "REVIEW"] <- "Classical mono (CCR2/IRF8-high)"
OPIS_MONO$Mono_label <- factor(lab)
OPIS_MONO$Mono_fine  <- droplevels(OPIS_MONO$Mono_fine)
OPIS_MONO$Mono_group <- droplevels(OPIS_MONO$Mono_group)
Idents(OPIS_MONO) <- group_var
message("Populations (", group_var, "): ",
        paste(levels(factor(OPIS_MONO[[group_var]][,1])), collapse = " | "))

# OUD level resolution
oud_levels <- levels(factor(OPIS_MONO$OUD_status))
oud_pos <- oud_levels[grepl("\\+|pos", oud_levels, ignore.case = TRUE)][1]
oud_neg <- setdiff(oud_levels, oud_pos)[1]
message("OUD contrast: ", oud_pos, " vs ", oud_neg)

# Composition: population x OUD (context for the OUD DGE)
comp <- as.data.frame.matrix(table(OPIS_MONO[[group_var]][,1], OPIS_MONO$OUD_status))
write.csv(comp %>% rownames_to_column("population"),
          file.path(dge.root, "population_x_OUD_counts.csv"), row.names = FALSE)
pct <- as.data.frame.matrix(round(prop.table(table(OPIS_MONO[[group_var]][,1],
                                                   OPIS_MONO$OUD_status), margin = 2) * 100, 2))
write.csv(pct %>% rownames_to_column("population"),
          file.path(dge.root, "population_x_OUD_percentWithinOUD.csv"), row.names = FALSE)

# =============================================================================
# 1b. ABUNDANCE / cluster distribution (SeuratExtend::ClusterDistrBar)
# =============================================================================
if (RUN_ABUNDANCE) {
  message("\n=== abundance / cluster distribution ===")
  cl  <- OPIS_MONO@meta.data[[group_var]]
  sv  <- function(p, file, w = 10, h = 6)
    ggsave(file.path(abund.dir, file), p, width = w, height = h, dpi = 300, bg = "white")
  
  # By donor (orig.ident): how each population is made up across donors
  sv(ClusterDistrBar(origin = OPIS_MONO$orig.ident, cluster = cl) +
       ggtitle("Population composition by donor (%)"),
     "abundance_byDonor_percent.png", w = 11)
  sv(ClusterDistrBar(origin = OPIS_MONO$orig.ident, cluster = cl, percent = FALSE) +
       ggtitle("Population counts by donor"),
     "abundance_byDonor_counts.png", w = 11)
  
  # By OUD status: is a population OUD-skewed?
  sv(ClusterDistrBar(origin = OPIS_MONO$OUD_status, cluster = cl) +
       ggtitle("OUD proportion within each population (%)"),
     "abundance_byOUD_percent.png")
  sv(ClusterDistrBar(origin = OPIS_MONO$OUD_status, cluster = cl, percent = FALSE) +
       ggtitle("Counts by OUD within each population"),
     "abundance_byOUD_counts.png")
  
  # Reverse: monocyte composition within each OUD group (compositional shift)
  sv(ClusterDistrBar(origin = cl, cluster = OPIS_MONO$OUD_status) +
       ggtitle("Monocyte composition within each OUD group (%)"),
     "composition_withinOUD_percent.png", w = 8)
  message("  abundance plots -> ", abund.dir)
}

# ---- shared significant-table writer ----------------------------------------
write_sig <- function(full, key_col, lfc, dir, prefix) {
  write.csv(full, file.path(dir, paste0(prefix, "_full.csv")), row.names = FALSE)
  sig <- full %>% filter(p_val_adj < 0.05, abs(avg_log2FC) > lfc) %>%
    mutate(direction = ifelse(avg_log2FC > 0, "up", "down")) %>%
    arrange(.data[[key_col]], desc(avg_log2FC))
  write.csv(sig, file.path(dir, paste0(prefix, "_significant.csv")), row.names = FALSE)
  message("  ", prefix, ": ", nrow(sig), " significant (up=",
          sum(sig$avg_log2FC > 0), ", down=", sum(sig$avg_log2FC < 0), ")")
  sig
}

# =============================================================================
# 2. IDENTITY DGE (RNA, one-vs-rest, MAST)
# =============================================================================
sig_id_rna <- NULL
if (RUN_IDENTITY_DGE) {
  message("\n=== Identity DGE (RNA, ", DGE_TEST, ") ===")
  if (DGE_TEST == "MAST" && !requireNamespace("MAST", quietly = TRUE))
    stop("MAST not installed: BiocManager::install('MAST') or set DGE_TEST <- 'wilcox'.")
  DefaultAssay(OPIS_MONO) <- "RNA"; Idents(OPIS_MONO) <- group_var
  full <- FindAllMarkers(OPIS_MONO, assay = "RNA", test.use = DGE_TEST, latent.vars = DGE_LATENT,
                         only.pos = FALSE, logfc.threshold = 0, min.pct = 0.1, verbose = TRUE)
  sig_id_rna <- write_sig(full, "cluster", dge_lfc, d_idg, "Identity_DGE_RNA")
}

# =============================================================================
# 3. IDENTITY DPE (ADT, one-vs-rest, wilcox on DSB data)
# =============================================================================
if (RUN_IDENTITY_DPE) {
  message("\n=== Identity DPE (ADT, wilcox) ===")
  DefaultAssay(OPIS_MONO) <- "ADT"; Idents(OPIS_MONO) <- group_var
  full <- FindAllMarkers(OPIS_MONO, assay = "ADT", test.use = "wilcox",
                         only.pos = FALSE, logfc.threshold = 0, min.pct = 0.1, verbose = FALSE)
  write_sig(full, "cluster", dpe_lfc, d_idp, "Identity_DPE_ADT")
  DefaultAssay(OPIS_MONO) <- "RNA"
}

# =============================================================================
# 4/5. OUD CONTRASTS (OUD+ vs OUD-), per population + all monocytes
# =============================================================================
run_oud <- function(assay, test, lfc, dir, prefix) {
  DefaultAssay(OPIS_MONO) <- assay
  pops <- levels(factor(OPIS_MONO[[group_var]][,1]))
  collect <- list()
  # all monocytes
  groups <- c("ALL_monocytes", pops)
  for (pop in groups) {
    sub <- if (pop == "ALL_monocytes") OPIS_MONO else
      OPIS_MONO[, OPIS_MONO[[group_var]][,1] == pop]
    tb <- table(sub$OUD_status)
    if (length(tb) < 2 || min(tb[c(oud_pos, oud_neg)], na.rm = TRUE) < oud_min_cells) {
      message("  skip ", pop, " (cells/group too few: ",
              paste(names(tb), tb, sep="=", collapse=", "), ")"); next
    }
    Idents(sub) <- "OUD_status"
    m <- tryCatch(FindMarkers(sub, ident.1 = oud_pos, ident.2 = oud_neg, assay = assay,
                              test.use = test, latent.vars = if (assay=="RNA") DGE_LATENT else NULL,
                              logfc.threshold = 0, min.pct = 0.1, verbose = FALSE),
                  error = function(e) { message("  ", pop, " failed: ", conditionMessage(e)); NULL })
    if (is.null(m) || nrow(m) == 0) next
    m$gene <- rownames(m); m$cluster <- pop          # 'cluster' col so write_sig works
    collect[[pop]] <- m
  }
  DefaultAssay(OPIS_MONO) <- "RNA"
  if (length(collect) == 0) { message("  no OUD contrasts produced for ", prefix); return(NULL) }
  full <- bind_rows(collect)
  write_sig(full, "cluster", lfc, dir, prefix)   # positive avg_log2FC = up in OUD+
}

sig_oud_rna <- NULL
if (RUN_OUD_DGE) {
  message("\n=== OUD DGE (RNA, ", DGE_TEST, "): ", oud_pos, " vs ", oud_neg, " ===")
  sig_oud_rna <- run_oud("RNA", DGE_TEST, dge_lfc, d_oug, "OUD_DGE_RNA")
}
if (RUN_OUD_DPE) {
  message("\n=== OUD DPE (ADT, wilcox): ", oud_pos, " vs ", oud_neg, " ===")
  run_oud("ADT", "wilcox", dpe_lfc, d_oup, "OUD_DPE_ADT")
}

# =============================================================================
# 6. ENRICHR pathway analysis on RNA DGE
# =============================================================================
databases <- c(
  "TRRUST_Transcription_Factors_2019", "ChEA_2022", "TRANSFAC_and_JASPAR_PWMs",
  "KEGG_2021_Human", "WikiPathways_2024_Human", "GO_Biological_Process_2023",
  "MSigDB_Hallmark_2020", "Panther_2016", "Reactome_2022", "BioPlanet_2019")
tf_databases      <- c("TRRUST_Transcription_Factors_2019", "ChEA_2022", "TRANSFAC_and_JASPAR_PWMs")
pathway_databases <- setdiff(databases, tf_databases)

init_enrichr <- function() {
  if (!requireNamespace("enrichR", quietly = TRUE) || !requireNamespace("openxlsx", quietly = TRUE)) {
    message("  enrichR/openxlsx not installed — skipping EnrichR."); return(FALSE) }
  suppressPackageStartupMessages({ library(enrichR); library(openxlsx) })  # attach sets base URL
  ok <- tryCatch({ enrichR::setEnrichrSite("Enrichr"); isTRUE(getOption("enrichR.live")) },
                 error = function(e) { message("  Enrichr error: ", conditionMessage(e)); FALSE })
  if (!ok) message("  Enrichr unreachable — skipping."); ok
}
run_enrichr_on_genes <- function(genes, label, base_output) {
  genes <- unique(genes[!is.na(genes) & nzchar(genes)])
  if (length(genes) < 5) { message("  <5 genes for ", label, " — skipped"); return(invisible()) }
  message("  - Enrichr: ", length(genes), " genes | ", label)
  enr <- enrichR::enrichr(genes, databases)
  dir.create(file.path(base_output, "CSVs"),  recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(base_output, "Plots"), recursive = TRUE, showWarnings = FALSE)
  wb <- openxlsx::createWorkbook()
  for (db in names(enr)) { openxlsx::addWorksheet(wb, substr(db,1,31)); openxlsx::writeData(wb, substr(db,1,31), enr[[db]]) }
  openxlsx::saveWorkbook(wb, file.path(base_output, "CSVs", paste0(label, "_Enrichment.xlsx")), overwrite = TRUE)
  top_tf <- list(); top_path <- list()
  for (db in names(enr)) {
    d <- enr[[db]]; if ("Combined Score" %in% colnames(d)) d <- dplyr::rename(d, Combined.Score = `Combined Score`)
    if (!"Combined.Score" %in% colnames(d)) next
    s <- d %>% filter(Adjusted.P.value < 0.05); if (nrow(s) == 0) next
    tt <- s %>% arrange(desc(Combined.Score)) %>% slice_head(n = 10) %>% mutate(Database = db)
    if (db %in% tf_databases) top_tf[[db]] <- tt else top_path[[db]] <- tt
  }
  mk_bar <- function(dfl, ttl, png) {
    df <- bind_rows(dfl); if (!"Combined.Score" %in% colnames(df) || nrow(df) == 0) return(invisible())
    df <- df %>% arrange(desc(Combined.Score)) %>% slice_head(n = 20)
    p <- ggplot(df, aes(reorder(Term, Combined.Score), Combined.Score, fill = Database)) +
      geom_col() + scale_y_log10() + coord_flip() + labs(title = ttl, x = NULL, y = "log10(Combined Score)") +
      theme_minimal() + theme(axis.text.y = element_text(size = 10))
    ggsave(png, p, width = 12, height = 10, dpi = 300, bg = "white")
  }
  mk_bar(top_tf,   paste("TFs -", label),      file.path(base_output, "Plots", paste0(label, "_TFs.png")))
  mk_bar(top_path, paste("Pathways -", label), file.path(base_output, "Plots", paste0(label, "_Pathways.png")))
}

if (RUN_ENRICHR) {
  message("\n=== EnrichR ===")
  if (!init_enrichr()) {
    message("  EnrichR unavailable — skipping pathway analysis.")
  } else {
    # Identity markers: all significant genes per population (both directions together)
    if (!is.null(sig_id_rna)) {
      for (pop in unique(sig_id_rna$cluster)) {
        g <- sig_id_rna$gene[sig_id_rna$cluster == pop]
        run_enrichr_on_genes(g, paste0("IDENTITY_", gsub("[^A-Za-z0-9]+", "_", pop)),
                             file.path(enr.root, "Identity"))
      }
    } else if (file.exists(file.path(d_idg, "Identity_DGE_RNA_significant.csv"))) {
      s <- read.csv(file.path(d_idg, "Identity_DGE_RNA_significant.csv"), check.names = FALSE)
      for (pop in unique(s$cluster))
        run_enrichr_on_genes(s$gene[s$cluster == pop],
                             paste0("IDENTITY_", gsub("[^A-Za-z0-9]+","_",pop)), file.path(enr.root, "Identity"))
    }
    # OUD: up and down separately (direction is meaningful: up = higher in OUD+)
    oud_sig <- if (!is.null(sig_oud_rna)) sig_oud_rna else
      tryCatch(read.csv(file.path(d_oug, "OUD_DGE_RNA_significant.csv"), check.names = FALSE),
               error = function(e) NULL)
    if (!is.null(oud_sig)) {
      for (pop in unique(oud_sig$cluster)) {
        tag <- gsub("[^A-Za-z0-9]+", "_", pop)
        run_enrichr_on_genes(oud_sig$gene[oud_sig$cluster == pop & oud_sig$avg_log2FC > 0],
                             paste0("OUD_", tag, "_UP"),   file.path(enr.root, "OUD"))
        run_enrichr_on_genes(oud_sig$gene[oud_sig$cluster == pop & oud_sig$avg_log2FC < 0],
                             paste0("OUD_", tag, "_DOWN"), file.path(enr.root, "OUD"))
      }
    }
  }
}

# =============================================================================
# 7. SAVE annotated object
# =============================================================================
qs_save(OPIS_MONO, file.path(load.path, "OPIS_MONO_Annotated.qs2"))
message("\nDone. DGE/DPE -> ", dge.root, " ; EnrichR -> ", enr.root)
message("Saved: ", file.path(load.path, "OPIS_MONO_Annotated.qs2"))

# -----------------------------------------------------------------------------
# CAVEAT (OUD contrasts): FindMarkers tests single cells, so cross-donor OUD
# effects are pseudoreplicated (donors, not cells, are the unit of replication).
# Treat the per-population OUD lists as exploratory; for manuscript claims confirm
# with a pseudobulk test (aggregate counts per donor x population, then
# DESeq2/edgeR/limma-voom on donor-level pseudobulk). Set DGE_LATENT <- "orig.ident"
# to at least adjust for donor in MAST as a partial mitigation.
# -----------------------------------------------------------------------------