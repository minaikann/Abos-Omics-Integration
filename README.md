# Integrating Liver and Plasma Omics with Data‑Driven MASLD Clusters

## Aim of the Project

The aim of this project is to investigate how data‑driven clinical subtypes of metabolic dysfunction‑associated steatotic liver disease (MASLD) relate to molecular profiles in the liver and blood.  
The project combines liver and plasma metabolomics with liver transcriptomics to:
  - Characterise metabolic signatures of cardiometabolic and liver‑specific MASLD clusters,
  - Identify signals associated with steatohepatitis (MASH) and fibrosis, and 
  - Explore multi‑omics signatures as possible candidates for non‑invasive biomarkers.

## Project Description

This project builds on previously defined MASLD clusters from the ABOS bariatric surgery cohort and integrates three layers of data: clinical and histological features, untargeted liver and plasma metabolomics, and liver transcriptomics.  
It uses statistical modelling, differential expression analysis and multi‑omics integration to describe MASLD heterogeneity and to identify biologically meaningful signatures of MASH and significant fibrosis across omics layers.
it contains full analysis workflow, including both in my thesis work and additional exploratory analyses that were not all reported in the manuscript.
Some of the code used is adapted from the original data‑driven cluster analysis by Violeta Raverdy et al. \
available at [this GitLab repository](https://gitlab.com/bilille/2024-raverdy_et_al-masld_clusters/-/tree/main).

## Methods / Approach

### 1. Cohort and data

- Participants from the ABOS bariatric surgery cohort with:
  - Baseline clinical and biochemical data
  - Liver histology (steatosis, inflammation, ballooning, fibrosis, MASH)
  - Liver untargeted metabolomics
  - Plasma untargeted metabolomics
  - Liver transcriptomics (RNA‑seq or microarray, depending on your setup)
- Three subtypes analysed: cardiometabolic (CM), liver‑specific (LS), and control (CTRL).

### 2. Metabolomics analyses

- Untargeted liver and plasma metabolomics (UPLC–MS/MS).
- Preprocessing: log‑transformation, imputation of missing values, removal of xenobiotics.
- Differential analysis:
  - Linear models with empirical Bayes moderation to compare CM vs CTRL, LS vs CTRL, and CM vs LS.
  - False discovery rate (FDR) correction on p‑values and log2 fold‑change thresholds.
- Pathway analysis:
  - Over‑representation analysis on significant metabolites using curated pathway databases.
- Histology associations:
  - Spearman correlations between liver metabolites and binary histological features.
  - Selection of metabolite sets associated with MASH and with significant fibrosis.

### 3. Liver Transcriptomics (DEA) analyses

- Preprocessing:
  - Normalisation and quality control of liver transcriptomic data.
- Differential expression analysis (DEA):
  - Comparison of gene expression between MASLD clusters (e.g. CM vs CTRL, LS vs CTRL, CM vs LS).
  - Identification of genes differentially expressed in relation to MASH and fibrosis status.
  - FDR correction to control for multiple testing.
- Functional interpretation:
  - Gene‑set / pathway enrichment analysis on differentially expressed genes to highlight key biological pathways (e.g. lipid metabolism, inflammation, extracellular matrix, oxidative stress).

### 4. Multi‑omics integration

- Metabolite‑only integration:
  - Joint analysis of liver and plasma metabolites to identify cross‑tissue metabolic signatures.
- Metabolite + transcriptomics integration:
  - Integration of liver transcriptomics with liver and plasma metabolomics to link gene‑level changes to metabolic reprogramming.
  - Identification of multi‑omics components or modules associated with MASLD clusters, MASH and fibrosis.
  - Exploration of how sphingolipid, phospholipid and bile‑acid pathways at the metabolite level align with transcriptional signatures in the same pathways.

### 5. Modelling of MASH and fibrosis

- LASSO logistic regression and related models to build parsimonious signatures for:
  - MASH vs non‑MASH
  - Significant fibrosis (F2–F4) vs non‑significant fibrosis (F0–F1)

## Technologies and Tools

- R (4.x)
- Typical packages:
  - `limma` (metabolite DEA and for transcriptomics DEA)
  - `glmnet` (LASSO )
  - `mixOmics` (multi‑omics integration, e.g. DIABLO)
  - `RaMP`  (metabolite Over expression pathway enrichment analysis)

## Key Findings 

- MASLD subtypes show distinct but overlapping metabolic reprogramming in liver and plasma, centred on sphingolipid, phospholipid and bile‑acid pathways.  
- Transcriptomic DEA highlights complementary changes in genes involved in lipid metabolism, inflammation and fibrosis‑related pathways, supporting the metabolomics findings.  
- Multi‑omics integration suggests that liver sphingolipid and redox‑related signatures, together with their matching transcriptional programs and plasma markers, could form the basis of mechanistically grounded, subtype‑aware biomarker panels for MASH and fibrosis.
