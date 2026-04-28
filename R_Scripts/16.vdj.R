###############################################################################
# OPIS TCR Pipeline
# ----------------------------------------------------------------------------
# Adapted from the TARA/EARTH infant longitudinal TCR script.
# OPIS is NOT longitudinal and has two subset Seurat objects (CD4 and CD8)
# instead of one per-cohort WNN object.
#
# Things you MAY need to edit (search for "EDIT:" in the script):
#   1. The contig-file path inside each sample folder (auto-detected, but
#      verify after the first run).
#   2. The cell-barcode rename regex in the "Merge with Seurat" section —
#      this depends on how your OPIS objects were named.
#   3. Any "Condition" / "Group" metadata column you want to stratify on.
###############################################################################

# ---- Load Libraries ---------------------------------------------------------
library(Seurat)
library(SingleCellExperiment)
library(scater)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(ggplot2)
library(gridExtra)
library(SeuratWrappers)
library(ggrepel)
library(patchwork)
library(scCustomize)
library(circlize)
library(ComplexHeatmap)
library(readxl)
library(scRepertoire)
library(igraph)
library(Cairo)
library(RColorBrewer)
library(Polychrome)
library(qs2)
library(Trex)

# ---- Paths ------------------------------------------------------------------
tcr.path  <- "/media/akshay-iyer/Elements/OPIS _Comparisons_Sequencing_Data/cellranger_out/TCR"
load.path <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/saved_R_data"
out.root  <- "~/Documents/OPIS_ECCITEseq/VDJ/TCR"

# Make sure output subdirectories exist
dir.create(out.root, recursive = TRUE, showWarnings = FALSE)
for (sub in c("Clonal_Visualizations", "CD3_Composition",
              "Clonal_Diversity", "Network_Analysis",
              "Seurat_Plots/Hyperexpansion", "Seurat_Plots/Clones_per_Cluster",
              "Seurat_Plots/Clonal_Overlay", "Trex")) {
  dir.create(file.path(out.root, sub), recursive = TRUE, showWarnings = FALSE)
}

###############################################################################
# 1. Read contig files -------------------------------------------------------
###############################################################################

# Each subfolder = one sample
sample_names <- list.dirs(tcr.path, full.names = FALSE, recursive = FALSE)
sample_names <- sample_names[sample_names != ""]
message("Samples detected: ", paste(sample_names, collapse = ", "))

# EDIT: auto-find filtered_contig_annotations.csv within each sample folder.
# If a sample has multiple matches (shouldn't, for TCR only), first match wins.
find_contig <- function(sample_dir) {
  hits <- list.files(sample_dir,
                     pattern = "^filtered_contig_annotations\\.csv$",
                     recursive = TRUE, full.names = TRUE)
  if (length(hits) == 0) stop("No filtered_contig_annotations.csv under ", sample_dir)
  hits[1]
}

contig_list_raw <- lapply(sample_names, function(s) {
  f <- find_contig(file.path(tcr.path, s))
  message("Reading ", s, "  <-  ", f)
  read.csv(f)
})
names(contig_list_raw) <- sample_names

# ---- Combine into scRepertoire object --------------------------------------
combined.TCR <- combineTCR(contig_list_raw, samples = sample_names)

# OPTIONAL: attach Condition/Group metadata for downstream stratification.
# EDIT: replace the placeholder below with a vector the same length
# as sample_names that encodes your groups (e.g. Responder / NonResponder,
# Treated / Untreated, etc.). For now everything is "OPIS".
condition_vec <- rep("OPIS", length(sample_names))
combined.TCR <- addVariable(combined.TCR,
                            variable.name = "Condition",
                            variables      = condition_vec)

###############################################################################
# 2. Basic Clonal Visualisations --------------------------------------------
###############################################################################
setwd(file.path(out.root, "Clonal_Visualizations"))

clonalQuant(combined.TCR, cloneCall = "strict", chain = "both", scale = FALSE)
ggsave("OPIS_Unique_Clones_Strict_TRAB_raw.png",    width = 18, height = 10)

clonalQuant(combined.TCR, cloneCall = "strict", chain = "both", scale = TRUE)
ggsave("OPIS_Unique_Clones_Strict_TRAB_scaled.png", width = 18, height = 10)

# If you add Condition/Group, you can split by it:
clonalQuant(combined.TCR, cloneCall = "strict", chain = "both",
            group.by = "Condition", scale = FALSE)
ggsave("OPIS_Unique_Clones_Strict_TRAB_byCondition.png", width = 12, height = 9)

clonalAbundance(combined.TCR, cloneCall = "strict", scale = FALSE)
ggsave("OPIS_Clonal_Abundance_raw.png",    width = 14, height = 10)
clonalAbundance(combined.TCR, cloneCall = "strict", scale = TRUE)
ggsave("OPIS_Clonal_Abundance_scaled.png", width = 14, height = 10)

clonalLength(combined.TCR, cloneCall = "aa", chain = "both")
ggsave("OPIS_Clonal_Length_raw.png",    width = 14, height = 10)
clonalLength(combined.TCR, cloneCall = "aa", chain = "both", scale = TRUE)
ggsave("OPIS_Clonal_Length_scaled.png", width = 14, height = 10)

clonalHomeostasis(combined.TCR, cloneCall = "strict",
                  cloneSize = c(Rare = 0.001, Small = 0.01, Medium = 0.1,
                                Large = 0.3, Hyperexpanded = 1))
ggsave("OPIS_Clonal_Homeostasis.png", width = 18, height = 10)

###############################################################################
# 3. CDR3 Composition --------------------------------------------------------
###############################################################################
setwd(file.path(out.root, "CD3_Composition"))

percentAA(combined.TCR, chain = "TRA", aa.length = 20)
ggsave("OPIS_Percent_AA_TRA.png", width = 22, height = 18)
percentAA(combined.TCR, chain = "TRB", aa.length = 20)
ggsave("OPIS_Percent_AA_TRB.png", width = 22, height = 18)

positionalEntropy(combined.TCR, chain = "both", aa.length = 20)
ggsave("OPIS_Positional_Entropy_TRAB.png", width = 18, height = 12)

for (g in c("TRAV", "TRBV", "TRAJ", "TRBJ")) {
  vizGenes(combined.TCR, x.axis = g, y.axis = NULL,
           plot = "heatmap", scale = TRUE)
  ggsave(paste0("OPIS_Heatmap_", g, "_gene.png"), width = 12, height = 9)
}

percentKmer(combined.TCR, cloneCall = "aa", chain = "TRA",
            motif.length = 3, top.motifs = 25)
ggsave("OPIS_Heatmap_TRA_kmer.png", width = 12, height = 9)
percentKmer(combined.TCR, cloneCall = "aa", chain = "TRB",
            motif.length = 3, top.motifs = 25)
ggsave("OPIS_Heatmap_TRB_kmer.png", width = 12, height = 9)

###############################################################################
# 4. Clonal Diversity & Overlap ---------------------------------------------
###############################################################################
setwd(file.path(out.root, "Clonal_Diversity"))

clonalDiversity(combined.TCR, cloneCall = "strict")
ggsave("OPIS_Clonal_Diversity_Strict.png", width = 12, height = 9)

clonalOverlap(combined.TCR, cloneCall = "strict", method = "morisita")
ggsave("OPIS_Clonal_Overlap_morisita.png", width = 14, height = 11)
clonalOverlap(combined.TCR, cloneCall = "strict", method = "raw")
ggsave("OPIS_Clonal_Overlap_raw.png",      width = 14, height = 11)

###############################################################################
# 5. Clonal Clustering Networks ---------------------------------------------
###############################################################################
setwd(file.path(out.root, "Network_Analysis"))

custom_colors <- c(
  "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
  "#FFFF33","#A65628","#F781BF","#999999","#66C2A5",
  "#FC8D62","#8DA0CB","#A6D854","#FFD92F","#E5C494"
)

set_colors_and_legend <- function(igraph.object) {
  color.legend     <- unique(igraph::V(igraph.object)$group)
  ordered_indices  <- order(nchar(color.legend), decreasing = TRUE)
  color.legend     <- color.legend[ordered_indices]
  col_samples      <- custom_colors[seq_along(color.legend)][ordered_indices]
  list(color.legend = color.legend, col_samples = col_samples)
}

plot_igraph <- function(igraph.object, col_samples, color.legend, title_text) {
  plot(igraph.object,
       vertex.size     = sqrt(igraph::V(igraph.object)$size),
       vertex.label    = NA,
       edge.arrow.size = 0,
       edge.curved     = 0.3,
       vertex.color    = col_samples)
  legend("topleft", legend = color.legend, pch = 16,
         col = unique(col_samples), bty = "n")
  title(title_text, cex.main = 1.5, font.main = 2)
}

# Build a single whole-cohort network per chain (you can also subset to
# specific sample groups if needed, as the template does).
igraph.TRA <- clonalCluster(combined.TCR, chain = "TRA", sequence = "aa",
                            group.by = "sample", threshold = 0.90,
                            exportGraph = TRUE)
igraph.TRB <- clonalCluster(combined.TCR, chain = "TRB", sequence = "aa",
                            group.by = "sample", threshold = 0.85,
                            exportGraph = TRUE)

TRA_colors <- set_colors_and_legend(igraph.TRA)
TRB_colors <- set_colors_and_legend(igraph.TRB)

CairoPNG("OPIS_TRA_network.png", width = 7200, height = 5500, res = 500)
plot_igraph(igraph.TRA, TRA_colors$col_samples, TRA_colors$color.legend,
            "OPIS TRA Sequences (Threshold 0.90)")
dev.off()

CairoPNG("OPIS_TRB_network.png", width = 7200, height = 5500, res = 500)
plot_igraph(igraph.TRB, TRB_colors$col_samples, TRB_colors$color.legend,
            "OPIS TRB Sequences (Threshold 0.85)")
dev.off()

###############################################################################
# 6. Load Seurat objects and attach TCR info --------------------------------
###############################################################################
OPIS_CD8 <- qs2::qs_read(file.path(load.path, "OPIS_CD8_WithPseudotime.qs2"))
OPIS_CD4 <- qs2::qs_read(file.path(load.path, "OPIS_CD4_WithPseudotime.qs2"))

# ---- EDIT: barcode renaming to match scRepertoire format -------------------
# scRepertoire prefixes cell barcodes with "<sample>_" (what you passed to
# combineTCR). Your Seurat objects must use the same convention.
#
# Inspect first: head(Cells(OPIS_CD8)); head(OPIS_CD8$orig.ident)
#
# The template used:
#   modified_barcodes <- gsub(".*_.*_.*_(.*)", "\\1", rownames(obj[[]]))
#   modified_barcodes <- paste0(obj$orig.ident, "_", modified_barcodes)
#
# That regex strips 3 underscore-separated prefixes. OPIS may have a different
# prefix count. Pick / adapt one of these options:

rename_for_scRep <- function(obj) {
  bc <- Cells(obj)
  # Option A (most common): barcodes already look like "AAACCTG...-1"; just prepend sample
  # bare_bc <- bc
  # Option B: strip everything up to & including the last underscore
  bare_bc <- sub(".*_", "", bc)
  new_bc  <- paste0(obj$orig.ident, "_", bare_bc)
  RenameCells(obj, new.names = new_bc)
}

OPIS_CD8 <- rename_for_scRep(OPIS_CD8)
OPIS_CD4 <- rename_for_scRep(OPIS_CD4)

# Sanity check: how many Seurat cells will find a TCR match?
rep_cells <- unlist(lapply(combined.TCR, function(x) x$barcode))
message("CD8 overlap with TCR: ",
        sum(Cells(OPIS_CD8) %in% rep_cells), " / ", ncol(OPIS_CD8))
message("CD4 overlap with TCR: ",
        sum(Cells(OPIS_CD4) %in% rep_cells), " / ", ncol(OPIS_CD4))
# If these are ~0, the rename_for_scRep() regex is wrong — adjust and rerun.

# ---- Attach clonal data -----------------------------------------------------
clone_size_bins <- c(Single = 1, Small = 5, Medium = 20,
                     Large = 100, Hyperexpanded = 500)

OPIS_CD8 <- combineExpression(combined.TCR, OPIS_CD8,
                              cloneCall  = "strict",
                              chain      = "both",
                              group.by   = "sample",
                              cloneSize  = clone_size_bins,
                              proportion = FALSE)

OPIS_CD4 <- combineExpression(combined.TCR, OPIS_CD4,
                              cloneCall  = "strict",
                              chain      = "both",
                              group.by   = "sample",
                              cloneSize  = clone_size_bins,
                              proportion = FALSE)

###############################################################################
# 7. Seurat Hyperexpansion Plots --------------------------------------------
###############################################################################
setwd(file.path(out.root, "Seurat_Plots/Hyperexpansion"))

colorblind_vector <- hcl.colors(n = 7, palette = "plasma", fixup = TRUE)

# EDIT: reduction name — template uses "wnn.umap"; your OPIS objects may use
# "umap" or "umap.harmony". Check with Reductions(OPIS_CD8).
red_name <- "wnn.umap"

make_hyperexp_plot <- function(obj, tag) {
  DimPlot_scCustom(obj, group.by = "cloneSize", reduction = red_name) +
    scale_color_manual(values = rev(colorblind_vector[c(1,3,4,5,7)]))
  ggsave(paste0("OPIS_", tag, "_Hyperexpansion.png"), width = 8, height = 7)
}
make_hyperexp_plot(OPIS_CD8, "CD8")
make_hyperexp_plot(OPIS_CD4, "CD4")

# Split by Condition if you populated it on the Seurat object
if ("Condition" %in% colnames(OPIS_CD8[[]])) {
  DimPlot_scCustom(OPIS_CD8, group.by = "cloneSize",
                   reduction = red_name, split.by = "Condition",
                   split_seurat = TRUE) +
    scale_color_manual(values = rev(colorblind_vector[c(1,3,4,5,7)]))
  ggsave("OPIS_CD8_Hyperexpansion_byCondition.png", width = 14, height = 7)
}
if ("Condition" %in% colnames(OPIS_CD4[[]])) {
  DimPlot_scCustom(OPIS_CD4, group.by = "cloneSize",
                   reduction = red_name, split.by = "Condition",
                   split_seurat = TRUE) +
    scale_color_manual(values = rev(colorblind_vector[c(1,3,4,5,7)]))
  ggsave("OPIS_CD4_Hyperexpansion_byCondition.png", width = 14, height = 7)
}

###############################################################################
# 8. Clones per Cluster -----------------------------------------------------
###############################################################################
setwd(file.path(out.root, "Seurat_Plots/Clones_per_Cluster"))

# EDIT: cluster column names — template used "snn.louvianmlr_1" and
# "predicted.celltype.l2". Adjust to whatever you have on OPIS_CD4/CD8.
cluster_cols <- c("seurat_clusters", "predicted.celltype.l2")

run_occupy <- function(obj, tag) {
  for (cc in cluster_cols) {
    if (!(cc %in% colnames(obj[[]]))) {
      message("Skipping ", cc, " (not in ", tag, ")")
      next
    }
    clonalOccupy(obj, x.axis = cc, label = FALSE)
    ggsave(paste0("OPIS_", tag, "_Clonal_Occupancy_", cc, ".png"),
           width = 20, height = 11)
    
    clonalOccupy(obj, x.axis = cc, proportion = TRUE, label = FALSE)
    ggsave(paste0("OPIS_", tag, "_Clonal_Occupancy_proportion_", cc, ".png"),
           width = 20, height = 11)
    
    tab <- clonalOccupy(obj, x.axis = cc, exportTable = TRUE)
    write.csv(tab,
              paste0("OPIS_", tag, "_Clones_per_Cluster_", cc, ".csv"),
              row.names = FALSE)
  }
}
run_occupy(OPIS_CD8, "CD8")
run_occupy(OPIS_CD4, "CD4")

###############################################################################
# 9. Clonal Overlay ---------------------------------------------------------
###############################################################################
setwd(file.path(out.root, "Seurat_Plots/Clonal_Overlay"))

clonalOverlay(OPIS_CD8, reduction = red_name,
              cutpoint = 1, bins = 10, facet.by = "orig.ident") +
  guides(color = "none")
ggsave("OPIS_CD8_Clonal_Overlay_By_Sample.png", width = 20, height = 14)

clonalOverlay(OPIS_CD4, reduction = red_name,
              cutpoint = 1, bins = 10, facet.by = "orig.ident") +
  guides(color = "none")
ggsave("OPIS_CD4_Clonal_Overlay_By_Sample.png", width = 20, height = 14)

if ("Condition" %in% colnames(OPIS_CD8[[]])) {
  clonalOverlay(OPIS_CD8, reduction = red_name, cutpoint = 10, bins = 25,
                facet.by = "Condition") + guides(color = "none")
  ggsave("OPIS_CD8_Clonal_Overlay_By_Condition.png", width = 16, height = 9)
}
if ("Condition" %in% colnames(OPIS_CD4[[]])) {
  clonalOverlay(OPIS_CD4, reduction = red_name, cutpoint = 10, bins = 25,
                facet.by = "Condition") + guides(color = "none")
  ggsave("OPIS_CD4_Clonal_Overlay_By_Condition.png", width = 16, height = 9)
}

###############################################################################
# 10. Trex — epitope annotation -------------------------------------------
###############################################################################
setwd(file.path(out.root, "Trex"))

annotate_and_tidy <- function(obj, tag) {
  ed0 <- annotateDB(obj, chains = "TRB")
  ed1 <- annotateDB(obj, chains = "TRB", edit.distance = 1)
  ed2 <- annotateDB(obj, chains = "TRB", edit.distance = 2)
  
  to_df <- function(o) {
    o@meta.data %>%
      dplyr::mutate(cells = rownames(.)) %>%
      dplyr::select(
        cells,
        PID       = orig.ident,
        any_of(c("Condition", "seurat_clusters", "predicted.celltype.l2")),
        CTstrict,
        clonalFrequency,
        TRB_Epitope.target,
        TRB_Epitope.sequence,
        TRB_Epitope.species,
        TRB_Tissue,
        TRB_Cell.type,
        TRB_Database
      ) %>%
      dplyr::filter(!is.na(clonalFrequency) & clonalFrequency > 1)
  }
  
  df0 <- to_df(ed0) %>% rename_with(~ paste0(., "_ED0"), starts_with("TRB_"))
  df1 <- to_df(ed1) %>% rename_with(~ paste0(., "_ED1"), starts_with("TRB_"))
  df2 <- to_df(ed2) %>% rename_with(~ paste0(., "_ED2"), starts_with("TRB_"))
  
  key_cols <- intersect(c("cells","PID","Condition","seurat_clusters",
                          "predicted.celltype.l2","CTstrict","clonalFrequency"),
                        colnames(df0))
  
  combined <- df0 %>%
    full_join(df1, by = key_cols) %>%
    full_join(df2, by = key_cols)
  
  write_csv(combined, paste0("OPIS_", tag, "_Trex_TRB_Epitope_Database.csv"))
  
  # Long format: epitope target
  long_target <- combined %>%
    pivot_longer(cols = starts_with("TRB_Epitope.target"),
                 names_to = "Edit_Distance",
                 names_pattern = "TRB_Epitope.target_(.*)",
                 values_to = "Epitope_Target") %>%
    mutate(Epitope_Target = ifelse(is.na(Epitope_Target),
                                   "Unknown", Epitope_Target)) %>%
    distinct(CTstrict, Edit_Distance, Epitope_Target, clonalFrequency)
  
  # Plot: all targets
  ggplot(long_target,
         aes(x = Epitope_Target, y = clonalFrequency, fill = CTstrict)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~ Edit_Distance, ncol = 1) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          legend.position = "none",
          strip.text = element_text(face = "bold")) +
    ylab("Clonal Frequency") + xlab("Epitope Target") +
    ggtitle(paste0("OPIS ", tag, ": Clonal Frequency per Epitope Target"))
  ggsave(paste0("OPIS_", tag, "_TRB_ClonalFreq_vs_Target.png"),
         width = 18, height = 13, dpi = 300, bg = "white")
  
  # Plot: excluding Unknown
  ggplot(long_target %>% filter(Epitope_Target != "Unknown"),
         aes(x = Epitope_Target, y = clonalFrequency, fill = CTstrict)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~ Edit_Distance, ncol = 1) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          legend.position = "none",
          strip.text = element_text(face = "bold")) +
    ylab("Clonal Frequency") + xlab("Epitope Target") +
    ggtitle(paste0("OPIS ", tag,
                   ": Clonal Frequency per Epitope Target (Excluding Unknown)"))
  ggsave(paste0("OPIS_", tag, "_TRB_ClonalFreq_vs_Target_noUnknown.png"),
         width = 18, height = 13, dpi = 300, bg = "white")
  
  # Per-PID summary table
  summary_by_PID <- long_target %>%
    left_join(combined %>% select(CTstrict, PID) %>% distinct(),
              by = "CTstrict") %>%
    group_by(Edit_Distance, PID, Epitope_Target) %>%
    summarise(Total_Clonal_Frequency = sum(clonalFrequency, na.rm = TRUE),
              .groups = "drop") %>%
    arrange(Edit_Distance, PID, desc(Total_Clonal_Frequency))
  write_csv(summary_by_PID,
            paste0("OPIS_", tag, "_Trex_Epitope_Specificity_By_Sample.csv"))
  
  # Species pie charts per PID (no Age facet — OPIS isn't longitudinal)
  long_target_full <- combined %>%
    pivot_longer(cols = starts_with("TRB_Epitope.target"),
                 names_to = "Edit_Distance",
                 names_pattern = "TRB_Epitope.target_(.*)",
                 values_to = "Epitope_Target")
  long_species_full <- combined %>%
    pivot_longer(cols = starts_with("TRB_Epitope.species"),
                 names_to = "Edit_Distance",
                 names_pattern = "TRB_Epitope.species_(.*)",
                 values_to = "Epitope_Species")
  long_species <- long_target_full %>%
    select(cells, PID, CTstrict, clonalFrequency,
           Edit_Distance, Epitope_Target) %>%
    left_join(long_species_full %>%
                select(cells, CTstrict, Edit_Distance, Epitope_Species),
              by = c("cells","CTstrict","Edit_Distance")) %>%
    mutate(Epitope_Target  = ifelse(is.na(Epitope_Target),  "Unknown", Epitope_Target),
           Epitope_Species = ifelse(is.na(Epitope_Species), "Unknown", Epitope_Species)) %>%
    distinct()
  
  pie_data <- long_species %>%
    group_by(PID, Edit_Distance, Epitope_Species) %>%
    summarise(Total_Clonal_Frequency = sum(clonalFrequency, na.rm = TRUE),
              .groups = "drop") %>%
    group_by(PID, Edit_Distance) %>%
    mutate(Percent = Total_Clonal_Frequency / sum(Total_Clonal_Frequency)) %>%
    ungroup() %>%
    mutate(Epitope_Species_short = sapply(strsplit(Epitope_Species, ";"),
                                          function(x) if (length(x) > 2) paste0(x[1], ";", x[2], "...")
                                          else paste(x, collapse = ";")))
  
  all_labels <- unique(pie_data$Epitope_Species_short)
  palette <- Polychrome::createPalette(
    N = length(all_labels),
    seedcolors = c("#E41A1C","#377EB8","#4DAF4A","#FF7F00"))
  names(palette) <- all_labels
  
  pie_dir <- file.path(out.root, "Trex", paste0(tag, "_PID_Pie_Chart"))
  dir.create(pie_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (pid in unique(pie_data$PID)) {
    plot_data <- pie_data %>% filter(PID == pid) %>% droplevels()
    if (nrow(plot_data) == 0) next
    n_cols <- length(unique(plot_data$Edit_Distance))
    
    p <- ggplot(plot_data,
                aes(x = "", y = Percent, fill = Epitope_Species_short)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar(theta = "y") +
      facet_grid(cols = vars(Edit_Distance)) +
      scale_fill_manual(values = palette) +
      theme_void(base_size = 14) +
      theme(legend.position = "bottom",
            legend.title = element_blank(),
            legend.text = element_text(size = 12),
            legend.key.height = unit(0.6, "lines"),
            strip.text.x = element_text(size = 14, face = "bold",
                                        margin = margin(t = 6)),
            plot.title = element_text(size = 16, face = "bold",
                                      hjust = 0.5, margin = margin(b = 12)),
            plot.margin = margin(t = 12, r = 10, b = 24, l = 10)) +
      guides(fill = guide_legend(nrow = 4, byrow = TRUE)) +
      ggtitle(paste0("OPIS ", tag, " — Epitope Species Specificity: ", pid))
    
    ggsave(file.path(pie_dir, paste0(pid, "_pie.png")),
           plot = p, bg = "white",
           width = max(5, n_cols * 5), height = 6, dpi = 300)
  }
  
  return(list(ed0 = ed0, ed1 = ed1, ed2 = ed2, table = combined))
}

OPIS_CD8_Trex <- annotate_and_tidy(OPIS_CD8, "CD8")
OPIS_CD4_Trex <- annotate_and_tidy(OPIS_CD4, "CD4")

###############################################################################
# 11. Save ----------------------------------------------------------------
###############################################################################
qs2::qs_save(list(
  OPIS_CD8     = OPIS_CD8,
  OPIS_CD8_ED0 = OPIS_CD8_Trex$ed0,
  OPIS_CD8_ED1 = OPIS_CD8_Trex$ed1,
  OPIS_CD8_ED2 = OPIS_CD8_Trex$ed2
), file = file.path(load.path, "OPIS_CD8_TCR_Combined.qs2"))

qs2::qs_save(list(
  OPIS_CD4     = OPIS_CD4,
  OPIS_CD4_ED0 = OPIS_CD4_Trex$ed0,
  OPIS_CD4_ED1 = OPIS_CD4_Trex$ed1,
  OPIS_CD4_ED2 = OPIS_CD4_Trex$ed2
), file = file.path(load.path, "OPIS_CD4_TCR_Combined.qs2"))

message("OPIS TCR pipeline complete.")

###############################################################################
# OPIS TCR Pipeline — Extensions
# ============================================================================
# Append after section 11 of the main TCR script (or source() after that
# script has finished running and the OPIS_CD4 / OPIS_CD8 / OPIS_CD?_Trex
# objects are in the environment).
#
# Adds:
#   12.  Config & helpers
#   13.  Per-participant log clone-size distribution
#   14.  Clonal expansion: OUD- vs OUD+ (stacked bar + UMAP contour)
#   15.  DGE: OUD+ clones vs OUD- clones (CD4 & CD8)
#   16.  DGE: clones vs non-clones in CD8 TEMRA & CD8 Innate-like
#   17.  HIV / CMV / Unknown specificity classification (helper)
#   18.  UMAP overlay: HIV-specific CD8 cells, split by OUD
#   19.  Cluster frequency of HIV-specific cells per CD8 cluster x OUD
#   20.  % cells in expanded clones per specificity, split by OUD
#   21.  DGE: HIV-specific vs CMV-specific cells (CD8)
#   22.  Waffle charts (replacement for the broken pie chart panel)
###############################################################################

# Required if you're sourcing this stand-alone (already loaded by main script):
suppressPackageStartupMessages({
  library(Seurat)
  library(scRepertoire)
  library(Trex)
  library(tidyverse)
  library(cowplot)
  library(scCustomize)
  library(Polychrome)
  library(scales)
  library(qs2)
})

###############################################################################
# OPTIONAL: STANDALONE LOAD --------------------------------------------------
# Uncomment this block if you're running this file independently (not pasted
# at the end of the main TCR script). Reads the qs2 files written at the end
# of the main script and rebuilds the variables this script expects.
###############################################################################

# load.path <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/saved_R_data"
# out.root  <- "~/Documents/OPIS_ECCITEseq/VDJ/TCR"
# red_name  <- "wnn.umap"   # change if your reduction is named differently
#
# cd8_pkg <- qs2::qs_read(file.path(load.path, "OPIS_CD8_TCR_Combined.qs2"))
# cd4_pkg <- qs2::qs_read(file.path(load.path, "OPIS_CD4_TCR_Combined.qs2"))
#
# OPIS_CD8 <- cd8_pkg$OPIS_CD8
# OPIS_CD4 <- cd4_pkg$OPIS_CD4
#
# OPIS_CD8_Trex <- list(ed0 = cd8_pkg$OPIS_CD8_ED0,
#                       ed1 = cd8_pkg$OPIS_CD8_ED1,
#                       ed2 = cd8_pkg$OPIS_CD8_ED2)
# OPIS_CD4_Trex <- list(ed0 = cd4_pkg$OPIS_CD4_ED0,
#                       ed1 = cd4_pkg$OPIS_CD4_ED1,
#                       ed2 = cd4_pkg$OPIS_CD4_ED2)
#
# # Rebuild the long-format `$table` that section 22 (waffle charts) expects.
# # Mirrors the to_df()/rename_with() block from the main script.
# build_combined_trex_table <- function(ed0, ed1, ed2) {
#   to_df <- function(o) {
#     o@meta.data %>%
#       dplyr::mutate(cells = rownames(.)) %>%
#       dplyr::select(cells,
#                     PID = orig.ident,
#                     dplyr::any_of(c("OUD_status","celltype_annotation",
#                                     "seurat_clusters")),
#                     CTstrict, clonalFrequency,
#                     dplyr::starts_with("TRB_")) %>%
#       dplyr::filter(!is.na(clonalFrequency) & clonalFrequency > 1)
#   }
#   df0 <- to_df(ed0) %>% dplyr::rename_with(~ paste0(., "_ED0"), dplyr::starts_with("TRB_"))
#   df1 <- to_df(ed1) %>% dplyr::rename_with(~ paste0(., "_ED1"), dplyr::starts_with("TRB_"))
#   df2 <- to_df(ed2) %>% dplyr::rename_with(~ paste0(., "_ED2"), dplyr::starts_with("TRB_"))
#   key_cols <- intersect(c("cells","PID","OUD_status","celltype_annotation",
#                           "seurat_clusters","CTstrict","clonalFrequency"),
#                         colnames(df0))
#   df0 %>% dplyr::full_join(df1, by = key_cols) %>%
#           dplyr::full_join(df2, by = key_cols)
# }
# OPIS_CD8_Trex$table <- build_combined_trex_table(
#   OPIS_CD8_Trex$ed0, OPIS_CD8_Trex$ed1, OPIS_CD8_Trex$ed2)
# OPIS_CD4_Trex$table <- build_combined_trex_table(
#   OPIS_CD4_Trex$ed0, OPIS_CD4_Trex$ed1, OPIS_CD4_Trex$ed2)

###############################################################################
# 12. Config & helpers ------------------------------------------------------
###############################################################################

# OUD metadata
oud_col      <- "OUD_status"
oud_pos      <- "OUD+"
oud_neg      <- "OUD-"
oud_palette  <- c(`OUD-` = "#3B7FB8", `OUD+` = "#D95F02")

# Cell-type annotation column (same one used in module-score script)
celltype_col <- "celltype_annotation"

# Specific clusters of interest for clones-vs-nonclones DGE
temra_label  <- "CD8+ TEMRA"
innate_label <- "Innate-like T"

# Trex object to use for specificity calls (ED0 = exact, ED1 = ≤1 mismatch, ED2 = ≤2)
# ED1 is a reasonable balance; change to ed0 for stricter calls
trex_obj_cd8 <- OPIS_CD8_Trex$ed1
trex_obj_cd4 <- OPIS_CD4_Trex$ed1

# Reduction (already set as red_name in main script; reuse if present)
if (!exists("red_name")) red_name <- "wnn.umap"

# Output dirs
oud_dir    <- file.path(out.root, "OUD_Comparisons");  dir.create(oud_dir,    recursive = TRUE, showWarnings = FALSE)
dge_dir    <- file.path(out.root, "DGE");              dir.create(dge_dir,    recursive = TRUE, showWarnings = FALSE)
spec_dir   <- file.path(out.root, "Specificity");      dir.create(spec_dir,   recursive = TRUE, showWarnings = FALSE)
waffle_dir <- file.path(out.root, "Trex", "Waffle");   dir.create(waffle_dir, recursive = TRUE, showWarnings = FALSE)

# --- Helper: classify cells as expanded / non-expanded -----------------------
# combineExpression() turns cloneSize into labels like "Single (0 < X <= 1)".
# Anything beyond the smallest bin counts as expanded.
mark_expanded <- function(obj) {
  cs <- as.character(obj$cloneSize)
  obj$has_TCR     <- !is.na(cs)
  obj$is_expanded <- obj$has_TCR & !grepl("^Single", cs)
  obj
}
OPIS_CD8     <- mark_expanded(OPIS_CD8)
OPIS_CD4     <- mark_expanded(OPIS_CD4)
trex_obj_cd8 <- mark_expanded(trex_obj_cd8)
trex_obj_cd4 <- mark_expanded(trex_obj_cd4)

# --- Helper: safe cluster-cell extraction ------------------------------------
cells_in_cluster <- function(obj, cluster_label) {
  rownames(obj@meta.data)[
    !is.na(obj@meta.data[[celltype_col]]) &
      obj@meta.data[[celltype_col]] == cluster_label
  ]
}

###############################################################################
# 13. Per-participant log clone-size distribution ---------------------------
###############################################################################
# For each PID, compute one row per unique clone with its size, normalize to
# clones-per-1000-cells in that PID, then log-transform.
# Plot one violin per PID coloured by OUD status.

clone_size_distribution <- function(obj, tag) {
  
  df <- obj@meta.data %>%
    dplyr::filter(!is.na(CTstrict)) %>%
    dplyr::group_by(orig.ident, CTstrict) %>%
    dplyr::summarise(n_cells = dplyr::n(), .groups = "drop_last") %>%
    dplyr::mutate(per_1k     = n_cells / sum(n_cells) * 1000,
                  log_n      = log10(n_cells),
                  log_per_1k = log10(per_1k)) %>%
    dplyr::ungroup()
  
  meta_pid <- obj@meta.data %>%
    dplyr::distinct(orig.ident, .data[[oud_col]])
  df <- df %>% dplyr::left_join(meta_pid, by = "orig.ident")
  
  write.csv(df,
            file.path(oud_dir, paste0("CloneSizeDistribution_", tag, ".csv")),
            row.names = FALSE)
  
  # Order PIDs by OUD group then by name
  pid_order <- meta_pid %>%
    dplyr::arrange(.data[[oud_col]], orig.ident) %>%
    dplyr::pull(orig.ident)
  df$orig.ident <- factor(df$orig.ident, levels = pid_order)
  
  p <- ggplot(df, aes(x = orig.ident, y = log_per_1k,
                      fill = .data[[oud_col]])) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.75,
                color = "grey20", linewidth = 0.3) +
    geom_jitter(width = 0.12, size = 0.4, alpha = 0.45, color = "grey25") +
    scale_fill_manual(values = oud_palette, name = "OUD") +
    labs(x = NULL,
         y = expression(log[10]("clones per 1000 cells")),
         title = paste0(tag,
                        ": clone-size distribution per participant")) +
    theme_cowplot(font_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title  = element_text(face = "bold", hjust = 0.5))
  
  ggsave(file.path(oud_dir, paste0("CloneSizeDistribution_", tag, ".png")),
         p, width = max(8, length(pid_order) * 0.45),
         height = 6, dpi = 300, bg = "white")
}
clone_size_distribution(OPIS_CD8, "CD8")
clone_size_distribution(OPIS_CD4, "CD4")

###############################################################################
# 14. Clonal expansion: OUD- vs OUD+ ----------------------------------------
###############################################################################
# (a) Stacked bar of cloneSize categories per OUD group

expansion_stacked <- function(obj, tag) {
  
  df <- obj@meta.data %>%
    dplyr::filter(!is.na(cloneSize)) %>%
    dplyr::count(.data[[oud_col]], cloneSize) %>%
    dplyr::group_by(.data[[oud_col]]) %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    dplyr::ungroup()
  
  desired <- c("Hyperexpanded","Large","Medium","Small","Single")
  bin_lvl <- levels(factor(df$cloneSize))
  ordered <- bin_lvl[order(match(
    sapply(strsplit(bin_lvl, " "), `[`, 1), desired))]
  df$cloneSize <- factor(df$cloneSize, levels = ordered)
  
  bin_pal <- rev(RColorBrewer::brewer.pal(
    n = max(3, length(ordered)), name = "Spectral"))[seq_along(ordered)]
  names(bin_pal) <- ordered
  
  p <- ggplot(df, aes(x = .data[[oud_col]], y = prop, fill = cloneSize)) +
    geom_col(width = 0.7, color = "grey15", linewidth = 0.3) +
    scale_y_continuous(labels = scales::percent_format(),
                       expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(values = bin_pal, name = "Clone size") +
    labs(x = NULL, y = "Proportion of cells",
         title = paste0(tag,
                        ": clonal expansion composition by OUD status")) +
    theme_cowplot(font_size = 14) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  ggsave(file.path(oud_dir,
                   paste0("ClonalExpansion_StackedBar_", tag, ".png")),
         p, width = 7, height = 6, dpi = 300, bg = "white")
}
expansion_stacked(OPIS_CD8, "CD8")
expansion_stacked(OPIS_CD4, "CD4")

# (b) UMAP contour overlay faceted by OUD
# cutpoint = 2 -> clonalOverlay keeps only cells with clonalFrequency >= 2,
# so singletons are excluded from the density estimate.
contour_by_oud <- function(obj, tag) {
  if (!(red_name %in% Reductions(obj))) {
    message("Reduction '", red_name, "' not found in ", tag,
            "; skipping contour."); return(invisible())
  }
  p <- clonalOverlay(obj, reduction = red_name,
                     cutpoint = 2, bins = 25, facet.by = oud_col) +
    guides(color = "none") +
    ggtitle(paste0(tag,
                   ": clonal density contours by OUD (singletons excluded)")) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  ggsave(file.path(oud_dir, paste0("ClonalOverlay_", tag, "_byOUD.png")),
         p, width = 14, height = 7, dpi = 300, bg = "white")
}
contour_by_oud(OPIS_CD8, "CD8")
contour_by_oud(OPIS_CD4, "CD4")

###############################################################################
# 15. DGE: OUD+ clones vs OUD- clones ---------------------------------------
###############################################################################
dge_oud_clones <- function(obj, tag) {
  
  cells_clonal <- rownames(obj@meta.data)[obj@meta.data$is_expanded]
  if (length(cells_clonal) < 10) {
    message("Too few clonal cells for ", tag, " OUD DGE."); return(NULL)
  }
  obj_c <- subset(obj, cells = cells_clonal)
  n_pos <- sum(obj_c@meta.data[[oud_col]] == oud_pos, na.rm = TRUE)
  n_neg <- sum(obj_c@meta.data[[oud_col]] == oud_neg, na.rm = TRUE)
  if (n_pos < 5 || n_neg < 5) {
    message("Too few cells per OUD group in ", tag,
            " (OUD+=", n_pos, ", OUD-=", n_neg, ")"); return(NULL)
  }
  
  Idents(obj_c) <- obj_c@meta.data[[oud_col]]
  m <- FindMarkers(obj_c, ident.1 = oud_pos, ident.2 = oud_neg,
                   logfc.threshold = 0.1, min.pct = 0.1)
  m$gene <- rownames(m)
  m <- m %>% dplyr::arrange(p_val_adj, dplyr::desc(abs(avg_log2FC)))
  
  write.csv(m, file.path(dge_dir,
                         paste0(tag, "_OUDpos_vs_OUDneg_clones.csv")), row.names = FALSE)
  m
}
dge_cd8_oud <- dge_oud_clones(OPIS_CD8, "CD8")
dge_cd4_oud <- dge_oud_clones(OPIS_CD4, "CD4")

###############################################################################
# 16. DGE: clones vs non-clones in CD8 TEMRA & CD8 Innate-like ---------------
###############################################################################
dge_clones_vs_nonclones <- function(obj, cluster_label, file_tag) {
  
  cells <- cells_in_cluster(obj, cluster_label)
  if (length(cells) < 20) {
    message("Cluster '", cluster_label, "' has only ", length(cells),
            " cells; skipping."); return(NULL)
  }
  obj_sub <- subset(obj, cells = cells)
  obj_sub <- subset(obj_sub,
                    cells = rownames(obj_sub@meta.data)[obj_sub$has_TCR])
  
  n_clo <- sum(obj_sub$is_expanded, na.rm = TRUE)
  n_non <- sum(!obj_sub$is_expanded, na.rm = TRUE)
  if (n_clo < 5 || n_non < 5) {
    message("Too few cells in ", cluster_label,
            " (clonal=", n_clo, ", non-clonal=", n_non, ")"); return(NULL)
  }
  
  obj_sub$clone_status <- ifelse(obj_sub$is_expanded, "Clonal", "NonClonal")
  Idents(obj_sub) <- obj_sub$clone_status
  m <- FindMarkers(obj_sub, ident.1 = "Clonal", ident.2 = "NonClonal",
                   logfc.threshold = 0.1, min.pct = 0.1)
  m$gene <- rownames(m)
  m <- m %>% dplyr::arrange(p_val_adj, dplyr::desc(abs(avg_log2FC)))
  
  write.csv(m, file.path(dge_dir,
                         paste0("CD8_", file_tag, "_Clonal_vs_NonClonal.csv")),
            row.names = FALSE)
  m
}
dge_temra  <- dge_clones_vs_nonclones(OPIS_CD8, temra_label,  "TEMRA")
dge_innate <- dge_clones_vs_nonclones(OPIS_CD8, innate_label, "InnateLike")

###############################################################################
# 17. HIV / CMV specificity classification ----------------------------------
###############################################################################
# Sets a 3-level factor `specificity` on a Trex object:
#   HIV   if TRB_Epitope.species mentions HIV
#   CMV   if TRB_Epitope.species mentions CMV / Cytomegalovirus
#   Unknown otherwise
# Cells with no TCR also fall in Unknown — keep that in mind for denominators.

classify_specificity <- function(obj) {
  meta <- obj@meta.data
  spec <- rep("Unknown", nrow(meta))
  
  if ("TRB_Epitope.species" %in% colnames(meta)) {
    sp <- as.character(meta[["TRB_Epitope.species"]])
    is_hiv <- grepl("HIV|Human immunodeficiency virus", sp, ignore.case = TRUE)
    is_cmv <- grepl("CMV|Cytomegalovirus",              sp, ignore.case = TRUE)
    spec[is_hiv]            <- "HIV"
    spec[is_cmv & !is_hiv]  <- "CMV"
  } else {
    warning("TRB_Epitope.species not found; all cells will be 'Unknown'.")
  }
  obj$specificity <- factor(spec, levels = c("HIV","CMV","Unknown"))
  obj
}
trex_obj_cd8 <- classify_specificity(trex_obj_cd8)
trex_obj_cd4 <- classify_specificity(trex_obj_cd4)

cat("\nCD8 specificity counts:\n");
print(table(trex_obj_cd8$specificity, trex_obj_cd8@meta.data[[oud_col]],
            useNA = "ifany"))
cat("\nCD4 specificity counts:\n");
print(table(trex_obj_cd4$specificity, trex_obj_cd4@meta.data[[oud_col]],
            useNA = "ifany"))

###############################################################################
# 18. UMAP overlay: HIV-specific CD8, split by OUD --------------------------
###############################################################################
hiv_cells_cd8 <- rownames(trex_obj_cd8@meta.data)[
  trex_obj_cd8$specificity == "HIV"]
cmv_cells_cd8 <- rownames(trex_obj_cd8@meta.data)[
  trex_obj_cd8$specificity == "CMV"]

if (length(hiv_cells_cd8) > 0 && red_name %in% Reductions(trex_obj_cd8)) {
  p <- DimPlot(trex_obj_cd8, reduction = red_name,
               cells.highlight = list(`HIV-specific` = hiv_cells_cd8),
               cols.highlight  = "#C4463A", cols = "grey85",
               sizes.highlight = 1.6, pt.size = 0.4,
               split.by = oud_col) +
    ggtitle(paste0("CD8: HIV-specific cells (n = ",
                   length(hiv_cells_cd8), ") by OUD")) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  ggsave(file.path(spec_dir, "HIVspecific_UMAP_CD8_byOUD.png"),
         p, width = 14, height = 7, dpi = 300, bg = "white")
} else {
  message("No HIV-specific CD8 cells found — skipping UMAP overlay.")
}

if (length(cmv_cells_cd8) > 0 && red_name %in% Reductions(trex_obj_cd8)) {
  p <- DimPlot(trex_obj_cd8, reduction = red_name,
               cells.highlight = list(`CMV-specific` = cmv_cells_cd8),
               cols.highlight  = "#3B7FB8", cols = "grey85",
               sizes.highlight = 1.6, pt.size = 0.4,
               split.by = oud_col) +
    ggtitle(paste0("CD8: CMV-specific cells (n = ",
                   length(cmv_cells_cd8), ") by OUD")) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  ggsave(file.path(spec_dir, "CMVspecific_UMAP_CD8_byOUD.png"),
         p, width = 14, height = 7, dpi = 300, bg = "white")
}

###############################################################################
# 19. Cluster frequency of HIV-specific cells per CD8 cluster x OUD ---------
###############################################################################
freq_df <- trex_obj_cd8@meta.data %>%
  dplyr::filter(!is.na(.data[[celltype_col]]),
                !is.na(.data[[oud_col]])) %>%
  dplyr::group_by(.data[[celltype_col]], .data[[oud_col]]) %>%
  dplyr::summarise(total_cells = dplyr::n(),
                   n_HIV   = sum(specificity == "HIV"),
                   n_CMV   = sum(specificity == "CMV"),
                   pct_HIV = 100 * n_HIV / total_cells,
                   pct_CMV = 100 * n_CMV / total_cells,
                   .groups = "drop")
write.csv(freq_df, file.path(spec_dir,
                             "CD8_HIV_CMV_freq_byCluster_OUD.csv"), row.names = FALSE)

p_hiv <- ggplot(freq_df,
                aes(x = .data[[celltype_col]], y = pct_HIV,
                    fill = .data[[oud_col]])) +
  geom_col(position = position_dodge(width = 0.8),
           width = 0.7, color = "grey15", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", pct_HIV)),
            position = position_dodge(width = 0.8),
            vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = oud_palette, name = "OUD") +
  labs(x = NULL, y = "% HIV-specific cells per cluster",
       title = "CD8: HIV-specific frequency per cluster, by OUD") +
  theme_cowplot(font_size = 14) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"),
        plot.title  = element_text(face = "bold", hjust = 0.5))
ggsave(file.path(spec_dir, "CD8_HIV_freq_byCluster_OUD.png"),
       p_hiv, width = 12, height = 6, dpi = 300, bg = "white")

# Optional: same plot for CMV — useful for context
p_cmv <- ggplot(freq_df,
                aes(x = .data[[celltype_col]], y = pct_CMV,
                    fill = .data[[oud_col]])) +
  geom_col(position = position_dodge(width = 0.8),
           width = 0.7, color = "grey15", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", pct_CMV)),
            position = position_dodge(width = 0.8),
            vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = oud_palette, name = "OUD") +
  labs(x = NULL, y = "% CMV-specific cells per cluster",
       title = "CD8: CMV-specific frequency per cluster, by OUD") +
  theme_cowplot(font_size = 14) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"),
        plot.title  = element_text(face = "bold", hjust = 0.5))
ggsave(file.path(spec_dir, "CD8_CMV_freq_byCluster_OUD.png"),
       p_cmv, width = 12, height = 6, dpi = 300, bg = "white")

###############################################################################
# 20. % cells in expanded clones per specificity, split by OUD --------------
###############################################################################
exp_df <- trex_obj_cd8@meta.data %>%
  dplyr::filter(has_TCR, !is.na(.data[[oud_col]])) %>%
  dplyr::group_by(.data[[oud_col]], specificity) %>%
  dplyr::summarise(total_cells  = dplyr::n(),
                   expanded     = sum(is_expanded),
                   pct_expanded = 100 * expanded / total_cells,
                   .groups = "drop")
write.csv(exp_df, file.path(spec_dir,
                            "CD8_PctExpanded_BySpecificity_OUD.csv"), row.names = FALSE)

p <- ggplot(exp_df, aes(x = specificity, y = pct_expanded,
                        fill = .data[[oud_col]])) +
  geom_col(position = position_dodge(width = 0.8),
           width = 0.7, color = "grey15", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", pct_expanded, total_cells)),
            position = position_dodge(width = 0.8),
            vjust = -0.2, size = 3.5, lineheight = 0.95) +
  scale_fill_manual(values = oud_palette, name = "OUD") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = "Antigen specificity", y = "% cells in expanded clones",
       title = "CD8: clonal expansion per specificity, by OUD") +
  theme_cowplot(font_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
ggsave(file.path(spec_dir, "CD8_PctExpanded_BySpecificity_OUD.png"),
       p, width = 8, height = 6, dpi = 300, bg = "white")

###############################################################################
# 21. DGE: HIV-specific vs CMV-specific cells (CD8) -------------------------
###############################################################################
if (length(hiv_cells_cd8) >= 5 && length(cmv_cells_cd8) >= 5) {
  Idents(trex_obj_cd8) <- trex_obj_cd8$specificity
  dge_hiv_cmv <- FindMarkers(trex_obj_cd8,
                             ident.1 = "HIV", ident.2 = "CMV",
                             logfc.threshold = 0.1, min.pct = 0.1)
  dge_hiv_cmv$gene <- rownames(dge_hiv_cmv)
  dge_hiv_cmv <- dge_hiv_cmv %>%
    dplyr::arrange(p_val_adj, dplyr::desc(abs(avg_log2FC)))
  write.csv(dge_hiv_cmv,
            file.path(dge_dir, "CD8_HIV_vs_CMV_DGE.csv"), row.names = FALSE)
} else {
  message("Skipping HIV vs CMV DGE: HIV=", length(hiv_cells_cd8),
          ", CMV=", length(cmv_cells_cd8))
}

###############################################################################
# 22. Waffle charts (replacement for the broken pie panel) ------------------
###############################################################################
# Each tile = 1% of cells (weighted by clonalFrequency) in that PID/edit-distance
# combination, colored by epitope species. Faceted by edit distance per PID.
# Uses base ggplot2 (no `waffle` package needed).

build_waffle_data <- function(combined_df) {
  
  long_target <- combined_df %>%
    tidyr::pivot_longer(cols = tidyr::starts_with("TRB_Epitope.species"),
                        names_to     = "Edit_Distance",
                        names_pattern = "TRB_Epitope.species_(.*)",
                        values_to    = "Epitope_Species") %>%
    dplyr::mutate(
      Epitope_Species = ifelse(is.na(Epitope_Species),
                               "Unknown", Epitope_Species),
      # Compress >2-element semicolon lists to keep legend manageable
      Epitope_Species_short = vapply(strsplit(Epitope_Species, ";"),
                                     function(x) if (length(x) > 2) paste0(x[1], ";", x[2], "...")
                                     else paste(x, collapse = ";"), character(1)))
  
  long_target %>%
    dplyr::group_by(PID, Edit_Distance, Epitope_Species_short) %>%
    dplyr::summarise(Total = sum(clonalFrequency, na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::group_by(PID, Edit_Distance) %>%
    dplyr::mutate(pct = 100 * Total / sum(Total)) %>%
    dplyr::ungroup()
}

# Build a 100-tile (10x10) grid for one PID, one edit distance
make_tile_grid <- function(species_vec, pct_vec) {
  # Largest-remainder rounding so tiles sum to exactly 100
  raw    <- pct_vec
  floors <- floor(raw)
  rem    <- raw - floors
  short  <- 100 - sum(floors)
  if (short > 0) {
    bumps <- order(rem, decreasing = TRUE)[seq_len(short)]
    floors[bumps] <- floors[bumps] + 1L
  }
  ord <- order(pct_vec, decreasing = TRUE)
  rep(species_vec[ord], floors[ord])
}

waffle_plot_one_pid <- function(pct_data, palette, pid, tag) {
  
  pid_data <- pct_data %>% dplyr::filter(PID == pid)
  if (nrow(pid_data) == 0) return(invisible(NULL))
  
  tile_df <- pid_data %>%
    dplyr::group_by(Edit_Distance) %>%
    dplyr::summarise(
      grid = list(make_tile_grid(Epitope_Species_short, pct)),
      .groups = "drop") %>%
    tidyr::unnest(cols = grid) %>%
    dplyr::group_by(Edit_Distance) %>%
    dplyr::mutate(idx = dplyr::row_number(),
                  x   = ((idx - 1) %% 10) + 1,
                  y   = floor((idx - 1) / 10) + 1) %>%
    dplyr::ungroup()
  
  ed_levels <- sort(unique(tile_df$Edit_Distance))
  tile_df$Edit_Distance <- factor(tile_df$Edit_Distance, levels = ed_levels)
  
  ggplot(tile_df, aes(x = x, y = y, fill = grid)) +
    geom_tile(color = "white", linewidth = 0.7) +
    facet_wrap(~ Edit_Distance, ncol = length(ed_levels)) +
    scale_fill_manual(values = palette, na.value = "grey92",
                      name = "Epitope species") +
    coord_equal() +
    scale_y_reverse() +
    labs(title = paste0("OPIS ", tag, " — ", pid),
         subtitle = "Each tile = 1% of cells (weighted by clonal frequency)",
         x = NULL, y = NULL) +
    theme_void(base_size = 13) +
    theme(strip.text       = element_text(size = 12, face = "bold",
                                          margin = margin(b = 6)),
          legend.position  = "bottom",
          legend.text      = element_text(size = 10),
          legend.key.size  = unit(0.6, "lines"),
          plot.title       = element_text(face = "bold", hjust = 0.5),
          plot.subtitle    = element_text(hjust = 0.5,
                                          margin = margin(b = 8)),
          plot.background  = element_rect(fill = "white", color = NA),
          panel.spacing    = unit(1.2, "lines")) +
    guides(fill = guide_legend(nrow = 4, byrow = TRUE,
                               override.aes = list(color = NA)))
}

run_waffles <- function(combined_table, tag) {
  
  pct_data <- build_waffle_data(combined_table)
  if (nrow(pct_data) == 0) {
    message("No waffle data for ", tag, "."); return(invisible())
  }
  
  all_labels <- sort(unique(pct_data$Epitope_Species_short))
  palette <- as.character(Polychrome::createPalette(
    N = max(length(all_labels), 3),
    seedcolors = c("#E41A1C","#377EB8","#4DAF4A","#FF7F00")))[
      seq_along(all_labels)]
  names(palette) <- all_labels
  # Force "Unknown" to a neutral grey
  if ("Unknown" %in% names(palette)) palette["Unknown"] <- "grey75"
  
  out_dir <- file.path(waffle_dir, tag)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  pids <- unique(pct_data$PID)
  for (pid in pids) {
    p <- waffle_plot_one_pid(pct_data, palette, pid, tag)
    if (is.null(p)) next
    n_facets <- length(unique(pct_data$Edit_Distance[pct_data$PID == pid]))
    ggsave(file.path(out_dir, paste0(pid, "_waffle.png")),
           p, width = max(6, n_facets * 4.2),
           height = 7, dpi = 300, bg = "white")
  }
  message("Wrote ", length(pids), " waffle charts -> ", out_dir)
}

run_waffles(OPIS_CD8_Trex$table, "CD8")
run_waffles(OPIS_CD4_Trex$table, "CD4")

message("\nOPIS TCR extension analyses complete.")
message("  OUD comparisons -> ", oud_dir)
message("  DGE             -> ", dge_dir)
message("  Specificity     -> ", spec_dir)
message("  Waffle charts   -> ", waffle_dir)
