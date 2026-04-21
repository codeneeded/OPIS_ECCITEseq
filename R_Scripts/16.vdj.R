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