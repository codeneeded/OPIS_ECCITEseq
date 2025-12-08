###Citeseq Pipeline
#Cite seurat and ds packages
#Load Required Libraries
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(hdf5r)
library(dsb)
library(data.table)
library(ggplot2)
library(hdf5r)
library(qs2) ### For more efficient save/loading

###############################################Path + GLobal Variables ##################################################
## Set path to load data
setwd("/media/akshay-iyer/Elements/OPIS _Comparisons_Sequencing_Data/cellranger_out/GEX_FB")  # Set working directory

# Define input and output paths
in.path <- "/media/akshay-iyer/Elements/OPIS _Comparisons_Sequencing_Data/cellranger_out/GEX_FB/"
out.path <- "/home/akshay-iyer/Documents/OPIS_ECCITEseq/saved_R_data/"  # Update as needed
r.str <- "/outs/raw_feature_bc_matrix.h5"
f.str <- "/outs/filtered_feature_bc_matrix.h5"

# Function to get all folder names within a specified path
get_folder_names <- function(in.path) {
  # List all directories within the specified path without recursion
  folder_names <- list.dirs(in.path, full.names = FALSE, recursive = FALSE)
  
  # Filter out the root path itself if present
  folder_names <- folder_names[folder_names != ""]
  
  return(folder_names)
}

# List folder names

f_names <- get_folder_names(in.path)

# Print the folder names
print(f_names)

####################################### Load Data + Basic Pre-processing ##############################################################
# Standard protein names (desired output)

panel_to_gene <- c(
  "Hu.CD101"           = "CD101",
  "Hu.CD103"           = "ITGAE",
  "Hu.CD105_43A3"      = "ENG",
  "Hu.CD107a"          = "LAMP1",
  "Hu.CD112"           = "NECTIN2",
  "Hu.CD119"           = "IFNGR1",
  "Hu.CD11a"           = "ITGAL",
  "Hu.CD11b"           = "ITGAM",
  "Hu.CD11c"           = "ITGAX",
  "Hu.CD122"           = "IL2RB",
  "Hu.CD123"           = "IL3RA",
  "Hu.CD124"           = "IL4R",
  "Hu.CD127"           = "IL7R",
  "Hu.CD13"            = "ANPEP",
  "Hu.CD134"           = "TNFRSF4",
  "Hu.CD137"           = "TNFRSF9",
  "Hu.CD141"           = "THBD",
  "Hu.CD146"           = "MCAM",
  "Hu.CD14_M5E2"       = "CD14",
  "Hu.CD152"           = "CTLA4",
  "Hu.CD154"           = "CD40LG",
  "Hu.CD155"           = "PVR",
  "Hu.CD158"           = "KIR2DL1",
  "Hu.CD158b"          = "KIR2DL3",
  "Hu.CD158e1"         = "KIR3DL1",
  "Hu.CD16"            = "FCGR3A",
  "Hu.CD161"           = "KLRB1",
  "Hu.CD163"           = "CD163",
  "Hu.CD169"           = "SIGLEC1",
  "Hu.CD18"            = "ITGB2",
  "Hu.CD183"           = "CXCR3",
  "Hu.CD185"           = "CXCR5",
  "Hu.CD19"            = "CD19",
  "Hu.CD194"           = "CCR4",
  "Hu.CD195"           = "CCR5",
  "Hu.CD196"           = "CCR6",
  "Hu.CD1c"            = "CD1C",
  "Hu.CD1d"            = "CD1D",
  "Hu.CD2"             = "CD2",
  "Hu.CD20_2H7"        = "MS4A1",
  "Hu.CD21"            = "CR2",
  "Hu.CD22"            = "CD22",
  "Hu.CD223"           = "LAG3",
  "Hu.CD224"           = "GGT1",
  "Hu.CD226_11A8"      = "CD226",
  "Hu.CD23"            = "FCER2",
  "Hu.CD24"            = "CD24",
  "Hu.CD244"           = "CD244",
  "Hu.CD25"            = "IL2RA",
  "Hu.CD26"            = "DPP4",
  "Hu.CD267"           = "TNFRSF13B",
  "Hu.CD268"           = "TNFRSF13C",
  "Hu.CD27"            = "CD27",
  "Hu.CD270"           = "TNFRSF14",
  "Hu.CD272"           = "BTLA",
  "Hu.CD274"           = "CD274",
  "Hu.CD279"           = "PDCD1",
  "Hu.CD28"            = "CD28",
  "Hu.CD29"            = "ITGB1",
  "Hu.CD303"           = "CLEC4C",
  "Hu.CD31"            = "PECAM1",
  "Hu.CD314"           = "KLRK1",
  "Hu.CD319"           = "SLAMF7",
  "Hu.CD32"            = "FCGR2A",
  "Hu.CD328"           = "SIGLEC7",
  "Hu.CD33"            = "CD33",
  "Hu.CD335"           = "NCR1",
  "Hu.CD35"            = "CR1",
  "Hu.CD352"           = "SLAMF6",
  "Hu.CD36"            = "CD36",
  "Hu.CD38_HIT2"       = "CD38",
  "Hu.CD39"            = "ENTPD1",
  "Hu.CD3_UCHT1"       = "CD3D",
  "Hu.CD40"            = "CD40",
  "Hu.CD41"            = "ITGA2B",
  "Hu.CD42b"           = "GP1BB",
  "Hu.CD45RA"          = "CD45RA",
  "Hu.CD45RO"          = "CD45RO",
  "Hu.CD45_HI30"       = "CD45",
  "Hu.CD47"            = "CD47",
  "Hu.CD48"            = "CD48",
  "Hu.CD49a"           = "ITGA1",
  "Hu.CD49b"           = "ITGA2",
  "Hu.CD49d"           = "ITGA4",
  "Hu.CD4_RPA.T4"      = "CD4",
  "Hu.CD5"             = "CD5",
  "Hu.CD52"            = "CD52",
  "Hu.CD54"            = "ICAM1",
  "Hu.CD56"            = "NCAM1",
  "Hu.CD57"            = "B3GAT1",
  "Hu.CD58"            = "CD58",
  "Hu.CD62L"           = "SELL",
  "Hu.CD62P"           = "SELP",
  "Hu.CD64"            = "FCGR1A",
  "Hu.CD69"            = "CD69",
  "Hu.CD7"             = "CD7",
  "Hu.CD71"            = "TFRC",
  "Hu.CD73"            = "NT5E",
  "Hu.CD79b"           = "CD79B",
  "Hu.CD8"             = "CD8A",
  "Hu.CD81"            = "CD81",
  "Hu.CD82"            = "CD82",
  "Hu.CD83"            = "CD83",
  "Hu.CD85j"           = "LILRB1",
  "Hu.CD86"            = "CD86",
  "Hu.CD88"            = "C5AR1",
  "Hu.CD94"            = "KLRD1",
  "Hu.CD95"            = "FAS",
  "Hu.CD99"            = "CD99",
  "Hu.CLEC12A"         = "CLEC12A",
  "Hu.CX3CR1"          = "CX3CR1",
  "Hu.FceRIa"          = "FCER1A",
  "Hu.GPR56"           = "ADGRG1",
  "Hu.HLA.ABC"         = "HLA-A",
  "Hu.HLA.DR"          = "HLA-DRA",
  "Hu.HLA.E"           = "HLA-E",
  "Hu.Ig.LightChain.k" = "IGKC",
  "Hu.Ig.LightChain.l" = "IGLC1", 
  "Hu.IgD"             = "IGHD",
  "Hu.IgM"             = "IGHM",
  "Hu.KLRG1"           = "KLRG1",
  "Hu.LOX.1"           = "OLR1",
  "Hu.TCR.AB"          = "TCR-AB",
  "Hu.TCR.Va7.2"       = "TCR-vA7.2",
  "Hu.TCR.Vd2"         = "TCR-vD2",
  "Hu.TIGIT"           = "TIGIT",
  "HuMs.CD44"          = "CD44",
  "HuMs.CD49f"         = "ITGA6",
  "HuMs.integrin.b7"   = "ITGB7",
  "HuMsRt.CD278"       = "ICOS",
  "Isotype_HTK888"     = "Hamster-IgG",
  "Isotype_MOPC.173"   = "Mouse-IgG2a",
  "Isotype_MOPC.21"    = "Mouse-IgG1",
  "Isotype_MPC.11"     = "Mouse-IgG2b",
  "Isotype_RTK2071"    = "Rat-IgG1",
  "Isotype_RTK2758"    = "Rat-IgG2a",
  "Isotype_RTK4530"    = "Rat-IgG2b"
)



# Loop through the files and standardize protein names
for (i in f_names) {
  # Read data from 10x outputs
  raw <- Read10X_h5(paste0(in.path, i, r.str))
  cells <- Read10X_h5(paste0(in.path, i, f.str))
  
  # Replace antibody capture names with standardized names based on position
  raw$`Antibody Capture`@Dimnames[[1]]   <- unname(panel_to_gene[ raw$`Antibody Capture`@Dimnames[[1]] ])
  cells$`Antibody Capture`@Dimnames[[1]] <- unname(panel_to_gene[ cells$`Antibody Capture`@Dimnames[[1]] ])
  # define a vector of cell-containing barcodes and remove them from unfiltered data 
  stained_cells <- colnames(cells$`Gene Expression`)
  background <- setdiff(colnames(raw$`Gene Expression`), stained_cells)
  
  # Assign final processed data
  prot <- raw$`Antibody Capture`
  rna <- raw$`Gene Expression`
  
  # Create metadata
  rna.size <- log10(Matrix::colSums(rna))
  prot.size <- log10(Matrix::colSums(prot))
  nCount_RNA <- Matrix::colSums(rna) 
  nCount_ADT <- Matrix::colSums(prot) 
  nFeature_RNA <- Matrix::colSums(rna > 0) 
  nFeature_ADT <- Matrix::colSums(prot > 0) 
  mtgene <- grep(pattern = "^MT-", rownames(rna), value = TRUE)
  mt.prop <- Matrix::colSums(rna[mtgene, ]) / Matrix::colSums(rna)
  
  # Combine metadata into a data frame
  md <- as.data.frame(cbind(nCount_RNA, nFeature_RNA, nCount_ADT, nFeature_ADT, rna.size, prot.size, mt.prop))
  
  # add indicator for barcodes Cell Ranger called as cells
  md$drop.class <- ifelse(rownames(md) %in% stained_cells, 'cell', 'background')
  
  # remove barcodes with no evidence of capture in the experiment
  md <- md[md$rna.size > 0 & md$prot.size > 0, ]
  
  # Save metadata and matrices
  assign(paste0(i, ".md"), md)
  assign(paste0(i, ".prot"), prot)
  assign(paste0(i, ".rna"), rna)
}


###  Droplet Settings
# Define the base path
base_path <- "~/Documents/OPIS_ECCITEseq"

# Define the directories to create
qc_dir <- file.path(base_path, "QC")
droplet_settings_dir <- file.path(qc_dir, "Droplet_Settings")

# Create directories if they do not exist
if (!dir.exists(qc_dir)) {
  dir.create(qc_dir, recursive = TRUE)
  cat("Created directory:", qc_dir, "\n")
}

if (!dir.exists(droplet_settings_dir)) {
  dir.create(droplet_settings_dir, recursive = TRUE)
  cat("Created directory:", droplet_settings_dir, "\n")
}

for (i in f_names) {
  # Set metadata object and protein object
  md <- eval(parse(text = paste0(i, ".md")))
  prot <- eval(parse(text = paste0(i, ".prot")))
  rna <- eval(parse(text = paste0(i, ".rna")))
  
  # Output plot for detected genes vs protein library size
  png(file.path(droplet_settings_dir, paste0(i, "_genevsprotlibsize.png")),
      width = 800, height = 600)
  
  p <- ggplot(md, aes(x = log10(nFeature_RNA), y = prot.size)) + 
    theme_bw() + 
    geom_bin2d(bins = 300) + 
    scale_fill_viridis_c(option = "C") + 
    facet_wrap(~drop.class)
  
  print(p)
  dev.off()
}

############################################### Droplet Thresholds for DSB Norm###################################################

# Define the function to process each entry
process_entry <- function(entry_name) {
  # Dynamically create the variable name for the data frame
  md_var_name <- paste0(entry_name, ".md")
  
  # Get the data frame from the variable name
  md <- get(md_var_name)
  
  # Filter the rows based on the given conditions
  filtered_md <- md[md$prot.size > 1.5 & md$prot.size < 4 & md$rna.size < 2.5, ]
  
  # Extract the row names of the filtered rows
  background_drops <- rownames(filtered_md)
  
  # Dynamically assign the result to the new variable
  assign(paste0(entry_name, ".background_drops"), background_drops, envir = .GlobalEnv)
}


# Loop through each entry and process it
for (entry_name in f_names) {
  process_entry(entry_name)
}


###################################Initial QC + DSB Normalization###################################################


for (i in f_names) {
  #Set metadata object and protein object to prevent constant eval parse calling  
  md <- eval(parse(text= paste0(i,'.md')))
  prot <- eval(parse(text = paste0(i,'.prot')))
  rna <- eval(parse(text = paste0(i,'.rna')))
  background_drops <- eval(parse(text = paste0(i,'.background_drops')))
  
  background.adt.mtx = as.matrix(prot[ , background_drops])
  
  cellmd = md[md$drop.class == 'cell', ]  #Define Cell Metadata on only cells not background  
  
  # filter drops with + / - 3 median absolute deviations from the median library size
  rna.mult = (3*mad(cellmd$rna.size))
  prot.mult = (3*mad(cellmd$prot.size))
  rna.lower = median(cellmd$rna.size) - rna.mult
  rna.upper = median(cellmd$rna.size) + rna.mult
  prot.lower = median(cellmd$prot.size) - prot.mult
  prot.upper = median(cellmd$prot.size) + prot.mult
  
  # filter rows based on droplet qualty control metrics
  qc_cells = rownames(
    cellmd[cellmd$prot.size > prot.lower & 
             cellmd$prot.size < prot.upper & 
             cellmd$rna.size > rna.lower & 
             cellmd$rna.size < rna.upper & 
             cellmd$mt.prop < 0.25, ]
  )
  
  # Output thresholds for quality control metrics as in any standard scRNAseq analysis
  png(paste0(i,'_qc_thresholds.png'),width = 800, height = 600)
  plot_aes = list(theme_bw(), geom_point(shape = 21 , stroke = 0, size = 0.7), scale_fill_viridis_c(option = "C"))
  p1 = ggplot(cellmd, aes(x = rna.size )) + geom_histogram(bins = 50) + theme_bw() + xlab("log10 RNA library size")
  p2 = ggplot(cellmd, aes(x = mt.prop)) + geom_histogram(bins = 50) + theme_bw() + xlab("mitochondrial read proportion")
  p3 = ggplot(cellmd, aes(x = log10(nFeature_RNA), y = rna.size, fill = mt.prop )) + plot_aes
  p4 = ggplot(cellmd, aes(x = nFeature_RNA, y = prot.size, fill = mt.prop )) + plot_aes
  print (p1+p2+p3+p4)
  dev.off()
  cell.adt.raw = as.matrix(prot[ , qc_cells])
  cell.rna.raw = rna[ ,qc_cells]
  cellmd = cellmd[qc_cells, ]
  
  #Proteins without Staining
  pm = sort(apply(cell.adt.raw, 1, max))
  pm2 = apply(background.adt.mtx, 1, max)
  head(pm2)
  
  #Assign it to your final output
  assign(paste0(i,'.cell.adt.raw'),cell.adt.raw)
  assign(paste0(i,'.cell.rna.raw'),cell.rna.raw)
  assign(paste0(i,'.background.adt.mtx'),background.adt.mtx)
  assign(paste0(i,'.cellmd'),cellmd)
  assign(paste0(i,'.pm'),pm)
}


# Check if you need to remove proteins without staining
#https://www.rdocumentation.org/packages/dsb/versions/0.3.0
#prot.expres.total <- rbindlist(adt.list)

#In this case we do not

### DSB Normalisation
#Set isotype control

isotype.controls <- c('Mouse-IgG1', 'Mouse-IgG2a','Mouse-IgG2b', 'Rat-IgG2b'
                      ,'Rat-IgG1', 'Rat-IgG2a', 'Hamster-IgG')




#normalize protein data for the cell containing droplets with the dsb method. 
# container to collect logs
drop_log <- list()

for (i in f_names) {
  message("Running DSB for: ", i)
  
  cell.adt.raw       <- eval(parse(text = paste0(i, ".cell.adt.raw")))
  background.adt.mtx <- eval(parse(text = paste0(i, ".background.adt.mtx")))
  
  ## --- 1) Compute max and 99th percentile ---
  max_cell <- apply(cell.adt.raw,       1, max)
  max_bg   <- apply(background.adt.mtx, 1, max)
  
  q99_cell <- apply(cell.adt.raw,       1, function(x) quantile(x, 0.99, na.rm = TRUE))
  q99_bg   <- apply(background.adt.mtx, 1, function(x) quantile(x, 0.99, na.rm = TRUE))
  
  ## --- 2) Identify proteins to drop ---
  zero_drop   <- (max_cell == 0 & max_bg == 0)
  lowq99_drop <- (q99_cell < 1 & q99_bg < 1)
  
  # combine criteria
  non_staining <- zero_drop | lowq99_drop
  
  # NEVER drop isotypes
  non_staining[rownames(cell.adt.raw) %in% isotype.controls] <- FALSE
  
  # extract names
  dropped <- rownames(cell.adt.raw)[non_staining]
  
  if (length(dropped) > 0) {
    message("Dropping proteins in ", i, ": ", paste(dropped, collapse = ", "))
    
    ## --- 3) Log dropped proteins with metadata ---
    drop_log[[i]] <- data.frame(
      sample      = i,
      protein     = dropped,
      reason      = ifelse(zero_drop[dropped], "zero", "low_q99"),
      max_cell    = max_cell[dropped],
      max_bg      = max_bg[dropped],
      q99_cell    = q99_cell[dropped],
      q99_bg      = q99_bg[dropped],
      stringsAsFactors = FALSE
    )
    
    ## --- 4) Actually drop them ---
    cell.adt.raw       <- cell.adt.raw[!non_staining, , drop = FALSE]
    background.adt.mtx <- background.adt.mtx[!non_staining, , drop = FALSE]
  }
  
  ## --- 5) Run DSB ---
  cells.dsb.norm <- DSBNormalizeProtein(
    cell_protein_matrix       = cell.adt.raw,
    empty_drop_matrix         = background.adt.mtx,
    denoise.counts            = TRUE,
    use.isotype.control       = TRUE,
    isotype.control.name.vec  = isotype.controls
  )
  
  cells.dsb.norm <- Matrix(as.matrix(cells.dsb.norm), sparse = TRUE)
  assign(paste0(i, ".cells.dsb.norm"), cells.dsb.norm)
}

setwd('/home/akshay-iyer/Documents/OPIS_ECCITEseq/QC/Droplet_Settings')

## --- 6) Save full drop log to CSV ---
if (length(drop_log) > 0) {
  drop_log_df <- do.call(rbind, drop_log)
  write.csv(drop_log_df, file = "DSB_dropped_proteins_log.csv", row.names = FALSE)
  message("Saved: DSB_dropped_proteins_log.csv")
} else {
  message("No proteins were dropped across samples.")
}



########################################Create + Merge Seurat Object (norm dsb+RNA)########################################

# Create Seurat Object

for (i in f_names) {
  cellmd <- eval(parse(text= paste0(i,'.cellmd')))
  cell.adt.raw <- eval(parse(text= paste0(i,'.cell.adt.raw')))
  cells.dsb.norm <- eval(parse(text= paste0(i,'.cells.dsb.norm')))
  cell.rna.raw <- eval(parse(text= paste0(i,'.cell.rna.raw')))
  
  # integrating with Seurat
  stopifnot(isTRUE(all.equal(rownames(cellmd), colnames(cell.adt.raw))))
  stopifnot(isTRUE(all.equal(rownames(cellmd), colnames(cell.rna.raw))))
  
  # create Seurat object note: min.cells is a gene filter, not a cell filter
  s = Seurat::CreateSeuratObject(counts = cell.rna.raw, 
                                 meta.data = cellmd,
                                 assay = "RNA", 
                                 min.cells = 20)
  
  # add dsb normalized matrix "dsb_norm_prot" to the "CITE" assay data slot
  s[["ADT"]] <- CreateAssayObject(data = cells.dsb.norm)
  s$orig.ident <- i
  #Assign
  assign(paste0(i,'.seurat'),s)
}


#Merge Seurat Objects

# Dynamically construct the merge() function
merged_seurat <- merge(
  x = eval(parse(text = paste0(f_names[1], ".seurat"))), 
  y = lapply(f_names[-1], function(name) eval(parse(text = paste0(name, ".seurat")))), 
  add.cell.id = f_names
)

###################################################### Edit Seurat Metadata ###############################################


### Add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

### Compute percent mitochondrial genes per cell
merged_seurat <- PercentageFeatureSet(merged_seurat, pattern = "^MT-", col.name = "percent_mito")

### Compute percentage of ribosomal genes per cell
merged_seurat <- PercentageFeatureSet(merged_seurat, pattern = "^RP[SL]", col.name = "percent_ribo")

# Percentage hemoglobin genes - includes all genes starting with HB except HBP.
merged_seurat <- PercentageFeatureSet(merged_seurat, "^HB[^(P)]", col.name = "percent_hb")

#Same for platelets
merged_seurat <- PercentageFeatureSet(merged_seurat, "PECAM1|PF4", col.name = "percent_plat")
levels(as.factor(merged_seurat$orig.ident))

# Define your OUD+ sample IDs
oud_positive <- c("S001", "S002", "S003", "S007", "S008")

# Add a new metadata column
merged_seurat$OUD_status <- ifelse(
  merged_seurat$orig.ident %in% oud_positive,
  "OUD+",
  "OUD-"
)

# Make it a factor (optional but useful)
merged_seurat$OUD_status <- factor(merged_seurat$OUD_status, levels = c("OUD-", "OUD+"))

# Create .RData object to load at any time
qs_save(merged_seurat, file=paste0(out.path,"Seuratv5_CITEseq_dsbnorm_merged_Seurat.qs2"))

########################################################################################################################