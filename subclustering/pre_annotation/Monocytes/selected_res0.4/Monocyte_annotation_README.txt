# OPIS ECCITEseq — Monocyte subclustering annotation (res 0.4)

Object: `OPIS_MONO_res0.4.qs2`  ·  10 subclusters (0–9), 7,280 cells
Parent clusters: 1 (Classical CD14+ mono), 15 (Classical CD14+ IL1B-high), 17 (CD16+ NCM)
Evidence used: curated RNA + ADT `FindAllMarkers`, signature module scores, cluster
correlation + dendrogram, per-cluster QC. WNN (RNA + DSB ADT) embedding; ADT averaged DSB-safe.

---

## Important caveats (read first)

1. **"Antigen_Presenting" is not a discriminating signature here.** It is built from
   HLA-DRA/DRB1/DPA1/DPB1/CD74/CIITA/CD86 — markers all monocytes express highly — so it
   tops the automated ranking for several clusters (0, 1, 7, 9) and the auto-generated
   `suggested_names.csv` is therefore misleading. All labels below are made from
   *cluster-specific* markers, not the AP signature.
2. **CD16 = ADT `FCGR3A`** (the panel uses gene-style protein names: FCGR3A, ITGAM=CD11b,
   ITGAX=CD11c, FCGR1A=CD64, HLA-DRA=HLA-DR). The curated ADT dot plot only rendered the 4
   names that matched literally; the calls below use the full ADT `FindAllMarkers` table.
3. **IFN/ISG is a weak axis** (raw module scores ≤0.19); the row-scaled signature heatmap
   exaggerates it. Cluster 5 is the only cluster with a real ISG program.

---

## Per-cluster annotation

### Classical monocytes (CD14+ CD16−) — the bulk: clusters 0, 1, 3
All three are CD14-protein-high (ADT CD14 strongly enriched, e.g. cluster 0 logFC 4.1),
classical-signature positive, CD16(FCGR3A)-negative. They differ only by an overlaid
activation/metabolic program, and are mutually 0.986–0.987 correlated with 115–137 DE genes
— i.e. one cell type in different states, not distinct identities.

- **0 (1435) — Classical monocyte, heat-shock/stress-high.** Top markers are almost
  entirely chaperones: HSPB1, HSPA1A/B, HSP90AA1/AB1, HSPD1, DNAJB1, BAG3, AHSA1, STIP1.
  A stress/dissociation program on a CD14+ backbone.
- **1 (1328) — Classical monocyte, VCAN+ / platelet-associated.** Classical core (VCAN,
  THBS1, PMP22) but carries a strong platelet signal (RNA PPBP/PF4/ITGA2B/GP9/TUBB1 top;
  ADT SELP/CD62P). Most likely monocyte–platelet aggregates (biologically relevant in
  HIV/OUD) and/or platelet ambient RNA — flag, but the cells are monocytes (keep).
- **3 (880) — Classical monocyte, hypoxic/glycolytic.** HILPDA, DDIT4, ADM, NDRG1, SLC2A3,
  HK2, PER1 = HIF/hypoxia + glycolysis program on a CD14+ backbone.

### 2 (1026) — Inflammatory monocyte
Highest Inflammatory score; defined by secreted chemokines/cytokines: CCL2, CCL7, CXCL1,
CXCL3, IL1B, SERPINB2, CLEC5A, C3AR1. ADT CD36-high. Very distinct from classical
(900–1500 DE genes vs classical clusters). Keep as its own population.

### 4 (726) — MHC-II-high antigen-presenting monocyte (cDC2-like lean)
Defining feature is high MHC-II: HLA-DQA1/DQB1, HLA-DPA1/DPB1, HLA-DMB, CD74 (RNA + ADT
HLA-DRA logFC 2.5), plus RETN, EPAS1. Carries a weak cDC2 lean (CLEC10A, CD1C). With
monocyte-only parents this is best read as an **MHC-II-high antigen-presenting monocyte**
rather than a bona fide cDC2 (canonical CLEC9A/XCR1/FCER1A are not convincingly present).

### 7 (313) — HLA-G+ (immunoregulatory) monocyte
Cleanly defined by HLA-G (top marker) with AHRR, RUNX3, ACP5, CYP1B1; ADT CD33 (logFC 7.1),
CD58, CD64(FCGR1A), CD226. CD16(FCGR3A)-negative, so classical-lineage but a distinct
HLA-G+ immunomodulatory state. Distinct enough to keep separate (does not merge with any
cluster at r ≥ 0.90).

### CD16+ monocytes (non-classical / intermediate): clusters 5, 8
The only clusters with a **negative** Classical score and high Intermediate + NonClassical
scores; ADT CD16(FCGR3A) is decisive. 5↔8 are 0.971 correlated but 190 DE genes apart, so
they are kept as two states of the CD16+ compartment.

- **5 (707) — CD16+ monocyte, ISG/IFN-high.** ADT FCGR3A(CD16)+ (pct 0.86 vs 0.05), IL3RA+,
  CD45RA+; RNA FCGR3A, VMO1, CDKN1C, MSR1, C1QA plus a clear type-I IFN program
  (OAS1, MX1, IFIT1/2/3, GBP4, IFITM1, RSAD2). The one genuine ISG cluster.
- **8 (170) — Non-classical monocyte (terminal, CD16++).** Strongest non-classical
  phenotype: ADT FCGR3A logFC 6.0 (pct 0.93), LILRB1 logFC 13.8, ITGAX(CD11c), CD86,
  PECAM1; RNA CDKN1C, VMO1, RHOC, IFITM2, CX3CR1, HLA-G. The canonical CD16++ CX3CR1+
  non-classical monocyte.

### 6 (574) — Classical monocyte, distinct (CCR2/IRF8-high) — **REVIEW**
Transcriptionally the most distinct cluster (dendrogram outgroup; correlation only
0.67–0.84 with all others), yet its top markers are sparse, non-canonical genes
(CAVIN2-AS1, CELF5, XYLB, GRK7, PLEKHA6, NUP210L, KCNMB4, SORCS2). Signatures: highest
Classical (0.73), high Inflammatory, with CCR2 (RNA), IRF8 and CIITA elevated; ADT is
monocyte-like (CD14, ITGAM, CD36, CD64) but weaker than cluster 0. QC is normal
(nFeature 2216, %mt 2.0). The atypical marker profile + outlier correlation is consistent
with either a genuine mono-DC/pre-DC-like state (CCR2/IRF8) or a technical (ambient/doublet)
cluster. **Action: review** — check donor/batch restriction and doublet score before
committing; provisional label "Classical monocyte (CCR2/IRF8-high)".

### 9 (121) — T-cell doublet / contamination — **EXCLUDE**
Unambiguous non-monocyte: ADT TCR-AB (logFC 8.9), CD5 (15.2), CD3D (10.3), CD2, CD7, CD8A,
CD4, CD28, CD27, ICOS, KLRG1; RNA dominated by mitochondrial genes (MT-ND6/ND4L/ND5/CYB/ATP8).
Monocyte–T doublets / T-cell contamination. Drop before downstream analysis.

---

## Proposed final annotation + combinations

| Cluster | n | Fine label | Merged group | Action |
|--------:|----:|------------|--------------|--------|
| 0 | 1435 | Classical monocyte (HSP/stress) | Classical monocyte | keep |
| 1 | 1328 | Classical monocyte (VCAN+/platelet-assoc) | Classical monocyte | keep (flag platelet) |
| 3 | 880 | Classical monocyte (hypoxic/glycolytic) | Classical monocyte | keep |
| 2 | 1026 | Inflammatory monocyte | Inflammatory monocyte | keep |
| 4 | 726 | MHC-II-high antigen-presenting monocyte | AP monocyte | keep |
| 7 | 313 | HLA-G+ monocyte | HLA-G+ monocyte | keep |
| 5 | 707 | CD16+ monocyte (ISG-high) | Non-classical/CD16+ monocyte | keep |
| 8 | 170 | Non-classical monocyte (CD16++) | Non-classical/CD16+ monocyte | keep |
| 6 | 574 | Classical monocyte (CCR2/IRF8-high) | review | review |
| 9 | 121 | T-cell doublet / contamination | — | EXCLUDE |

### Why these combinations

- **Merge 0 + 1 + 3 → "Classical monocyte".** Correlation 0.986–0.987, 115–137 DE genes;
  identical CD14+/CD16− core differing only by stress/metabolic programs (HSP, platelet,
  hypoxia). For a clean figure these collapse to one Classical population; keep the
  sub-states only if the stress/platelet programs are of specific interest. This matches the
  res-0.5 over-splitting we saw earlier — classical fragments into states, not identities.
- **Group 5 + 8 → "Non-classical / CD16+ monocyte"** (two states: 5 = ISG-high, 8 = terminal
  CD16++). They are the only CD16(FCGR3A)+ / Classical-negative clusters. Keep as two if the
  ISG program matters; collapse to one CD16+ population otherwise.
- **Keep 2, 4, 7 separate** — each is a distinct, marker-defined state (inflammatory;
  MHC-II-high AP; HLA-G+), well-separated on the dendrogram and by DE count.
- **Exclude 9** (T doublets). **Review 6** before assigning.

This yields a defensible 5-population scheme (Classical / Inflammatory / Antigen-presenting /
HLA-G+ / Non-classical-CD16+), with cluster 6 pending review and cluster 9 removed.

---

## Drop-in annotation map (for the next script)

```r
# fine label per subcluster (res 0.4)
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

# coarse merged group for figures
mono.group <- c(
  "0" = "Classical monocyte", "1" = "Classical monocyte", "3" = "Classical monocyte",
  "2" = "Inflammatory monocyte",
  "4" = "Antigen-presenting monocyte",
  "7" = "HLA-G+ monocyte",
  "5" = "Non-classical/CD16+ monocyte", "8" = "Non-classical/CD16+ monocyte",
  "6" = "REVIEW", "9" = "EXCLUDE"
)

# apply (drop the excluded cluster before downstream)
OPIS_MONO$Mono_fine  <- factor(unname(mono.fine [as.character(OPIS_MONO$Subcluster_ID)]))
OPIS_MONO$Mono_group <- factor(unname(mono.group[as.character(OPIS_MONO$Subcluster_ID)]))
OPIS_MONO <- subset(OPIS_MONO, subset = Subcluster_ID != "9")   # remove T-doublets
```

---

## Follow-ups worth doing before finalising
- **Cluster 6:** check donor/batch composition and a doublet score (scDblFinder); if it is
  donor-restricted or high-doublet, treat as technical, otherwise keep as a real state.
- **Cluster 1 platelet signal:** decide monocyte–platelet aggregates (biology) vs ambient
  (technical); a PPBP/PF4 vs monocyte-gene co-expression check on single cells settles it.
- **Optional automated cross-check:** flip `RUN_SINGLER <- TRUE` (celldex Monaco) to confirm
  classical / intermediate / non-classical calls independently.
