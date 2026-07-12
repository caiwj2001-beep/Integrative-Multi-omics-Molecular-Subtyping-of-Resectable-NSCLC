# Analysis Code — Integrative Multi-omics Molecular Subtyping of Resectable NSCLC

This repository contains the complete data-analysis pipeline for the study
*"Integrative Multi-omics Analysis Defines Molecular Subtypes of Resectable
Non–Small Cell Lung Cancer in the Adjuvant Setting."*

The code reproduces every quantitative result, statistic, and figure in the
manuscript from publicly available data. Manuscript/cover-letter assembly
scripts are intentionally excluded; only the research computation is included.

--------------------------------------------------------------------------------
## 1. Data sources (all public)

| Cohort | Role | Source |
|--------|------|--------|
| ALCHEMIST-ALCH | Discovery (multi-omics) | NCI Genomic Data Commons (project ALCHEMIST-ALCH); controlled clinical elements under dbGaP phs001140 |
| MSK-IMPACT (pan-cancer 2026) | Validation (genomic axis + OS) | cBioPortal |
| TCGA-LUAD / TCGA-LUSC (PanCancer Atlas) | Validation (transcriptomic axis, stage-adjusted OS/DFS) | cBioPortal studies `luad_tcga_pan_can_atlas_2018`, `lusc_tcga_pan_can_atlas_2018` |

Raw data are NOT redistributed here. See section 4 to obtain them.

--------------------------------------------------------------------------------
## 2. Software environment

- R >= 4.6.1
- Python >= 3.11 (data preparation, cBioPortal/GDC API calls, MOFA backend)
- R packages: data.table, survival, survminer, matrixStats, ggplot2, patchwork,
  ComplexHeatmap, circlize, maftools, ConsensusClusterPlus, GSVA, MOFA2,
  limma, edgeR, estimate, msigdbr, randomForest, irr, RColorBrewer
- Python packages: pandas, numpy, mofapy2 (MOFA2 backend via reticulate),
  requests/curl for API access

Install R packages with `00_setup_packages.R`. MOFA2 uses the Python `mofapy2`
backend through reticulate (see `05_mofa_integration.R`).

--------------------------------------------------------------------------------
## 3. Pipeline (run in order)

All scripts read paths from `config.R` (R) / `config.py` (Python). Set
`PROJECT_ROOT` there once; everything else is relative.

| Step | Script | Purpose |
|------|--------|---------|
| 00 | `00_setup_packages.R` | Install all R/Bioconductor packages |
| 01 | `01_prepare_data.py` | GDC file→patient UUID mapping; build expression / mutation / CNA / miRNA matrices + clinical master table (ALCHEMIST); build MSK-IMPACT NSCLC matrices |
| 02 | `02_fetch_tcga.py` | Retrieve TCGA-LUAD/LUSC genomic, expression, clinical (stage, OS, DFS) via cBioPortal API |
| 03 | `03_mutation_landscape.R` | Somatic mutation landscape, driver frequencies (maftools) |
| 04 | `04_cna_analysis.R` | Copy-number calls, fraction genome altered, focal events |
| 05 | `05_expression_clustering.R` | Consensus clustering of expression (ConsensusClusterPlus) |
| 06 | `06_immune_pathway.R` | ESTIMATE scores + Hallmark ssGSEA (GSVA) |
| 07 | `07_mofa_integration.R` | MOFA2 multi-omics integration → integrated subtypes (IMS1–4) |
| 08 | `08_subtype_definition.R` | Subtype biological characterization + statistics |
| 09 | `09_subtype_classifier.R` | Random-forest extension of subtypes to the 935-tumor core set |
| 10 | `10_msk_validation.R` | MSK-IMPACT genomic-axis transfer + OS (Cox) |
| 11 | `11_msk_treatment_sensitivity.R` | Treatment-era sensitivity analysis (panel proxy) |
| 12 | `12_tcga_build_signature.R` | Expression signature markers (limma) for transcriptomic transfer; writes `results/signature_genes.csv` |
| 12b | `02b_fetch_tcga_signature_expr.py` | Fetch TCGA mRNA expression for the signature genes → `tcga_signature_expr_matrix.csv` (run after step 12) |
| 13 | `13_tcga_validation.R` | TCGA genomic-axis transfer, stage-balance, stage-adjusted Cox |
| 14 | `14_tcga_transcriptomic_transfer.R` | TCGA transcriptomic-axis transfer + stage-adjusted OS/DFS |
| 15 | `15_dual_axis_concordance.R` | Genomic vs transcriptomic axis concordance (Cohen's kappa) |
| 16 | `16_figures.R` | Main Figures 1–6 (cohort, oncoplot, heatmap, biology, drivers, MSK KM) |
| 17 | `17_figures_tcga_concordance.R` | Figures 7 (TCGA), 8 (dual-axis), S1 (consensus diagnostics) |
| 18 | `18_tables.R` | Tables 1–3 and Supplementary Tables S1–S4 (statistics) |

--------------------------------------------------------------------------------
## 4. Reproduction notes

1. Download raw data from the sources in section 1 into the paths defined in
   `config.py` (`ALCHEMIST_DIR`, `MSK_DIR`). ALCHEMIST open-access files are
   fetched from GDC; MSK-IMPACT and TCGA from cBioPortal.
2. Run `00` (once), then `01`–`02` (data preparation), then `03`–`18` in order.
3. Intermediate objects are written to `results/`; figures to `figures/`;
   tables to `tables/`. These directories are created automatically.
4. Random seeds are fixed (`set.seed(20260707)` / `20260709`) for the
   stochastic steps (consensus clustering, MOFA, random forest).

--------------------------------------------------------------------------------
## 5. Key parameters

- Consensus clustering: PAM, Pearson distance, 500–1000 resamples, 80–90% subsampling
- MOFA2: 15 factors, Bernoulli likelihood for the mutation view, medium convergence
- Subtype number K = 4 (selected by proportion of ambiguous clustering + interpretability)
- Survival: Kaplan-Meier + log-rank; multivariable Cox (stage/histology adjusted in TCGA)
- Multiple testing: Benjamini-Hochberg FDR

Corresponding author: Wen-Jie Cai (caiwj2001@126.com).
