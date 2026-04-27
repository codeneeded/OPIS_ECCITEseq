# OPIS ECCITE-seq — Single-Cell Multi-Omic Immune Profiling in HIV and Opioid Use Disorder

> ECCITE-seq (paired scRNA-seq + surface protein + TCR profiling) analysis of immune cell composition, transcriptional states, and T cell clonotype dynamics in people with HIV and opioid use disorder. Code repository for a manuscript in preparation.

---

## Overview

This repository contains the R analysis pipelines for a single-cell multi-omic study applying **ECCITE-seq** (Enhanced CRISPR-compatible CITE-seq: simultaneous scRNA-seq + antibody-derived tag [ADT] surface protein profiling + TCR V(D)J sequencing) to PBMCs from participants enrolled in the **OPIS (Opioid Immunity Study)** at the University of Miami.

The OPIS study (NCT04304768) enrolled HIV-infected and HIV-uninfected adults with and without opioid use disorder (OUD) to characterize how chronic opioid exposure shapes immune function — particularly in the context of HIV-driven immune dysregulation. This repository extends the published serological and functional findings from the OPIS cohort to the single-cell transcriptomic and phenotypic level, providing mechanistic resolution into how OUD and HIV status interact to remodel the immune landscape.

**Manuscript in preparation.** No publication data are available yet.

---

## Study Cohort — OPIS (NCT04304768)

| Group | Description |
|---|---|
| **HIV+ OP+** | HIV-infected, opioid-dependent |
| **HIV+ OP−** | HIV-infected, non-opioid user |
| **HIV− OP+** | HIV-uninfected, opioid-dependent |
| **HIV− OP−** | HIV-uninfected, non-opioid user (controls) |

- **Sponsor:** University of Miami — Pahwa Laboratory
- **Intervention:** Participants received seasonal quadrivalent influenza vaccination (Fluzone Quadrivalent)
- **Primary question:** How does chronic opioid use modulate immune cell composition, transcriptional programs, and T cell clonotype dynamics — and how does HIV status interact with OUD to shape these effects?

---

## Repository Structure

```
OPIS_ECCITEseq/
│
├── R_Scripts/                         # Main Seurat v5 analysis pipelines
├── QC/                                # Quality control outputs
│                                      # (filtering thresholds, library metrics,
│                                      #  doublet removal)
├── Integration/                       # Multi-sample batch correction and
│                                      # multimodal WNN integration
├── Annotation/                        # Cell type annotation outputs
│                                      # (UMAPs, marker heatmaps, cluster labels)
├── subclustering/                     # Fine-grained subclustering of major
│                                      # immune lineages (T, B, NK, monocyte)
├── Differential_Expression/           # Differential gene expression results
│                                      # across HIV and OUD status groups
├── Pathway_Analysis_EnrichR/          # Pathway enrichment analyses (EnrichR)
├── VDJ/TCR/                          # TCR V(D)J clonotype analysis
│                                      # (scRepertoire)
├── Manuscript_Figures/               # Publication-ready figures
├── grant_plots/                       # Figures generated for grant applications
│
├── Isotype.xlsx                       # Isotype control data for DSB ADT
│                                      # normalization
├── OPIS_ECCITEseq.Rproj              # RStudio project file
└── README.md
```

---

## Analysis Pipeline

### 1. Quality Control (`QC/`)
- Per-sample filtering on minimum gene count, maximum mitochondrial RNA fraction
- Library-level QC metrics (median genes/cell, total UMI, ADT detection rate)
- Doublet detection and removal

### 2. Normalization & Integration (`Integration/`)
- **RNA:** Log-normalization and scaling
- **ADT (surface protein):** DSB (denoised and scaled by background) normalization using isotype controls from `Isotype.xlsx`
- Batch correction across participants and sequencing runs
- Multimodal integration using **Seurat v5 WNN** (weighted nearest neighbor), combining RNA and ADT modalities into a joint embedding

### 3. Cell Type Annotation (`Annotation/`)
- Broad PBMC lineage annotation using canonical RNA markers and ADT surface protein co-expression
- Major populations: CD4+ T cells, CD8+ T cells, B cells, NK cells, monocytes, dendritic cells, and other lineages
- Cluster resolution optimized and validated against reference-based predictions

### 4. Subclustering (`subclustering/`)
- Fine-grained resolution of major immune lineages to resolve functional subpopulations:
  - **T cells:** Naïve, central memory, effector memory, TEMRA, exhausted, regulatory subsets
  - **B cells:** Naïve, memory, plasmablasts, transitional subsets
  - **NK cells:** CD56bright, CD56dim, adaptive/memory-like subsets
  - **Monocytes:** Classical, intermediate, non-classical subsets

### 5. Differential Expression (`Differential_Expression/`)
- Pairwise DGE comparisons across the four OPIS groups (HIV × OUD status)
- Primary contrasts: HIV+ OP+ vs HIV+ OP−; HIV− OP+ vs HIV− OP−; HIV+ vs HIV−; OP+ vs OP−
- Results exported per cell type and per comparison

### 6. Pathway Analysis (`Pathway_Analysis_EnrichR/`)
- Gene set enrichment of DGE results using **EnrichR**
- Databases queried: GO Biological Process, KEGG, Reactome, WikiPathways
- Pathway-level visualization across immune cell populations and group comparisons

### 7. TCR Repertoire Analysis (`VDJ/TCR/`)
- V(D)J clonotype assembly and integration with scRNA-seq data using **scRepertoire**
- Clonal expansion analysis across HIV and OUD groups
- Clonotype sharing between groups and across cell populations
- Diversity metrics (Shannon entropy, clonal dominance)

---

## Scientific Questions

1. **Immune cell composition** — How do OUD and HIV status, independently and interactively, alter the frequencies of major and minor PBMC subpopulations?
2. **Transcriptional programs** — Which gene expression programs (activation, exhaustion, inflammation, interferon signaling) are dysregulated by OUD and/or HIV at the single-cell level?
3. **Surface phenotype** — How do ADT-resolved surface protein profiles distinguish immune states across the four groups beyond what RNA alone captures?
4. **T cell clonotypes** — Does OUD alter T cell clonal expansion, diversity, or clonotype composition in the presence or absence of HIV?
5. **Pathway-level effects** — Which biological pathways are enriched in DGE signatures, and do OUD and HIV converge on shared or distinct mechanisms?
6. **Mechanistic basis of vaccine responses** — Can single-cell immune profiling provide a mechanistic explanation for the previously observed effect of OUD on influenza vaccine antibody responses?

---

## Related Publication

The following peer-reviewed publication from the OPIS cohort provides the serological and functional context for this single-cell study:

- **Dang CM, Nelson CM, Pahwa RN, Tookes HE, Feaster DJ, Singh P, Rodriguez AE, Forrest DW, Nakamura N, Ghanta PP, Jayaweera DT, Iyer A, Pallikkuth S, Pahwa SG.** (2025). Chronic opioid use is associated with higher antibody response to influenza vaccination in people living with HIV. *Frontiers in Immunology*, 16, 1686103. [DOI: 10.3389/fimmu.2025.1686103](https://doi.org/10.3389/fimmu.2025.1686103)
  > Serological analysis of HAI antibody titres across the four OPIS groups, demonstrating a paradoxical association between OUD and elevated influenza vaccine antibody responses. The ECCITE-seq analyses in this repository are designed to provide single-cell mechanistic resolution of this finding.

---

## Dependencies

All scripts are written in **R**. Key packages:

| Package | Purpose |
|---|---|
| `Seurat` (v5) | WNN multimodal clustering, DGE, visualization |
| `dsb` | ADT normalization with isotype controls |
| `scRepertoire` | TCR clonotype analysis |
| `enrichR` | Pathway enrichment analysis |
| `dplyr` / `ggplot2` / `patchwork` | Data wrangling and visualization |
| `pheatmap` / `ggrepel` | Heatmaps and labeled plots |

---

## Status

**Manuscript in preparation.** This repository is under active development. Code, figures, and data will be updated as the project progresses. Raw sequencing data will be deposited at NCBI GEO upon publication.

> ⚠️ Raw sequencing data and participant-level metadata are not included in this repository due to privacy constraints. Contact the study team for data access inquiries.

---

## Citation

A full citation will be added upon publication. In the interim, please acknowledge this repository and the OPIS study (NCT04304768) and the Pahwa Laboratory, University of Miami Miller School of Medicine.

---

## Contact

For questions, please open an [issue](https://github.com/codeneeded/OPIS_ECCITEseq/issues) in this repository.
