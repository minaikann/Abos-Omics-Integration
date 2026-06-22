---
title: "ABOS Omics Analysis"
author: "Deborah Mina Ikann"
date: "`r Sys.Date()`"
output: html_document
  
---


```{r setup, include=FALSE}
knitr::opts_chunk$set(
  message = FALSE,
  warning = FALSE)
```

### Libaries used
```{r}
library(RaMP)
library(glmnet)
library(pROC)
library(caret)
library(rms)
library(janitor)
library(MASS)
library(readxl)
library(ggplot2)
library(sqldf)
library(cluster)
library(dplyr)
library(fmsb)
library(DT)
library(fpc)
library(limma)
library(data.table)
library(openxlsx)
library(ggpubr)
library(ggrepel)
library(tibble)
library(eulerr)
library(stringr)
library(tidyverse)
library(corrplot)
library(reshape2)
library(UpSetR)
library(patchwork)
library(tidyr)
library(clusterProfiler)
```


### Functions used 

```{r}
source("~/Documents/Clustering_ABOS/codes/functions/pie_chart_fct.R")
source("~/Documents/Clustering_ABOS/codes/functions/radar_chart_fct.R")
source("~/Documents/Clustering_ABOS/codes/functions/barplot_test_fct.R")
source("~/Documents/Clustering_ABOS/codes/functions/khiFish_KW_fcts.R")
source("~/Documents/Clustering_ABOS/codes/functions/functions_ liver_metabo.R")
source("~/Documents/Clustering_ABOS/codes/functions/metabolites_correlation_fct.R")
source("~/Documents/Clustering_ABOS/codes/functions/enrichment_analysis_fct.R")
source("~/Documents/Clustering_ABOS/codes/functions/lasso_model_fct.R")
```

### Data Preprocessing 

```{r merging_mash_&clinical_data}
# Clean and standardize column names from MASH study dataset
mash_data_study <- mash_data %>%
  rename(
    Abos_numeros_inclusion = IEX_PATIENT..Abo_numeros_inclusion,  # Patient inclusion ID
    sexe = IEX_PATIENT..ABOS_sexe,                                # Sex
    fibrosis = M000..Biopsie_Kleiner_Fibrose,                     # Fibrosis stage
    Steatosis = M000..Biopsie_SteatoseSeverite__lcn,              # Steatosis severity
    Inflammation = M000..Biopsie_Kleiner_Nb_Infla_Lobu,           # Lobular inflammation
    Ballooning = M000..Biopsie_Kleiner_Ballonisation,             # Hepatocyte ballooning
    NAS = M000..Biopsie_NAS                                        # NAFLD Activity Score
  ) %>%
  # Convert ballooning from categorical to numeric scale FIRST
  mutate(Ballooning = ifelse(Ballooning == "Quelques", 1,
                             ifelse(Ballooning == "Nombreux", 2, 0)),
         # Create binary MASH diagnosis based on histological triad
         Mash = ifelse(Steatosis >= 1 & Ballooning >= 1 & Inflammation >= 1, 1, 0))

# FIXED: Handle dates and create properly aligned data frame
date_data <- mash_data_study %>%
  mutate(InterventionDate_Clean = as.Date(ABOS_intervention_date, format = "%Y-%m-%d")) %>%
  group_by(Abos_numeros_inclusion) %>%
  arrange(InterventionDate_Clean) %>%
  slice(1) %>%
  ungroup() %>%
  select(Abos_numeros_inclusion, InterventionDate_Clean) %>%  # KEEP Abos_numeros_inclusion!
  as.data.frame() %>%
  column_to_rownames("Abos_numeros_inclusion")  # Now works!

# Select MASH-relevant columns with proper row alignment
mash_subset <- mash_data_study %>%
  select(Abos_numeros_inclusion, Steatosis, Ballooning, Inflammation, Mash, fibrosis) %>%
  # Take earliest biopsy per patient
  group_by(Abos_numeros_inclusion) %>%
  slice(1) %>%
  ungroup() %>%
  as.data.frame() %>%
  column_to_rownames("Abos_numeros_inclusion")

# Find common patients
common_rows <- intersect(rownames(clinical_data), rownames(mash_subset))
cat("Matched", length(common_rows), "patients out of", nrow(clinical_data), "clinical patients\n")
cat("Matched", length(common_rows), "patients out of", nrow(mash_subset), "MASH patients\n")

# Add MASH columns to clinical_data (NA for unmatched patients)
clinical_data$Steatosis <- mash_subset[match(rownames(clinical_data), rownames(mash_subset)), "Steatosis"]
clinical_data$Ballooning <- mash_subset[match(rownames(clinical_data), rownames(mash_subset)), "Ballooning"]
clinical_data$Inflammation <- mash_subset[match(rownames(clinical_data), rownames(mash_subset)), "Inflammation"]
clinical_data$Mash <- mash_subset[match(rownames(clinical_data), rownames(mash_subset)), "Mash"]
clinical_data$fibrosis <- mash_subset[match(rownames(clinical_data), rownames(mash_subset)), "fibrosis"]
clinical_data$InterventionDate_Clean <- date_data[match(rownames(clinical_data), rownames(date_data)), "InterventionDate_Clean"]


# Clean fibrosis to first character (if it's a string)
clinical_data$Biopsie_Kleiner_Fibrose <- substr(clinical_data$Biopsie_Kleiner_Fibrose, 1, 1)

# Save results
#write.xlsx(clinical_data, "clinical_data_Mash.xlsx", rowNames = TRUE)


# Matched 1545 patients out of 1545 clinical patients
# Matched 1545 patients out of 1777 MASH patients
```

```{r clean-qc}
# Define key clinical variables for analysis
var <- c("TGPUL",              # ALT (liver enzyme)
         "hba1cpc",            # HbA1c (%)
         "triglymmolL",        # Triglycerides (mmol/L)
         "BMI",                 # Body Mass Index
         "AgeJourIntervention", # Age at intervention (days)
         "cholest_LDLmmolL")    # LDL cholesterol (mmol/L)

# Total number of patients in cohort
no_row <- nrow(clinical_data)
cat("Total patients:", no_row, "\n")

# Sort clinical data by numeric sample ID for consistency
clinical_data <- clinical_data[order(as.numeric(rownames(clinical_data))), ]

# Alcohol  consumption filter 
clinical_data <- clinical_data[is.na(clinical_data$M0_alcool_quantite) | 
                              (clinical_data$sexe == "F" & clinical_data$M0_alcool_quantite <= 50) |
                              (clinical_data$sexe == "M" & clinical_data$M0_alcool_quantite <= 60), ]

exl_alch <- no_row - nrow(clinical_data)
cat("Excluded alcohol consumption:", exl_alch,"\n")

# BMI > 30 (obese cohort)
clinical_data <- clinical_data[clinical_data$BMI > 30, ]
no_row_after_alch <- no_row - exl_alch
exl_bmi <- no_row_after_alch - nrow(clinical_data)
cat("Excluded BMI ≤30:", exl_bmi,"\n")

# Cluster vars only
cluster_data <- clinical_data[, var]
rows_after_bmi <- nrow(cluster_data)
cat("After BMI filter:", rows_after_bmi,"\n")

# Remove rows with ANY NA in var
cluster_data <- na.omit(cluster_data)
exl_na <- rows_after_bmi - nrow(cluster_data)
cat("Excluded NA rows:", exl_na,"\n")


# Total patients: 1545 
# Excluded alcohol consumption: 54 
# Excluded BMI ≤30: 58 
# After BMI filter: 1433 
# Excluded NA rows: 27 
```

### Scaling and removal of outliers

```{r scaling}
# Scale all features to z-scores (mean = 0, SD = 1) for clustering
# Essential for fair comparison across metabolites with different scales
cluster_data_scale <- scale(cluster_data)

# Identify extreme outliers (> |5| SD from mean)
# These could unduly influence clustering results
outlier_mask <- apply(cluster_data_scale, 1, function(row) {
  any(abs(row) >= 5, na.rm = TRUE)
})

# Count number of samples flagged as outliers
exl_outliers <- sum(outlier_mask)


#Excluding outlier in the scaled data that will be used for the clustering 
cluster_data_scale <- cluster_data_scale[!outlier_mask, ]  

#Excluding outlier in the raw data that will be used for the Plotting of clinical variables  
cluster_data<- cluster_data[!outlier_mask, ]


cat("Clean rows for in cluster data :", nrow(cluster_data), "\n")
cat("Excluded outliers:", exl_outliers, "\n")
cat("Clean rows for clustering:", nrow(cluster_data_scale), "\n")
cat("Clean rows for plotting:", nrow(cluster_data), "\n")


# Clean rows for in cluster data : 1389 
# Excluded outliers: 17 
# Clean rows for clustering: 1389 
# Clean rows for plotting: 1389 
```

### K-Medoids

```{r kmediods_test}
nboot <- 100
K <- 3:10
silhouette <- numeric(length(K))

for(i in seq_along(K)) {
  pam_mod <- pam(cluster_data_scale, k = K[i], nstart = nboot)
  silhouette[i] <- mean(silhouette(pam_mod)[,3])
}

plot(K, silhouette, type = "o", pch = 19, xlab = "Number of clusters", ylab = "Mean Silhouette", main = "Optimal k")
abline(v = which.max(silhouette) + 2, col = "red", lty = 2)  # K starts at 3

```

### Jaccard Index

```{r jaccard_index}
set.seed(24)
jacc <- fpc::clusterboot(data = cluster_data_scale, B = 2000, clustermethod = pamkCBI,  krange = 6, count = FALSE)
meanJacc <- data.frame(meanJacc = round(apply(t(jacc$bootresult), 2, mean), 3))
meanJacc <- meanJacc$meanJacc

#The mean and SD
meanJaccAll <- round(mean(meanJacc),2)
sdJaccAll <- round(sd(meanJacc),2)
 
cat("The jaccard index number ", meanJaccAll," and the sd ", sdJaccAll)

# The jaccard index number  0.73  and the sd  0.07
``` 
 

 
```{r kmedoids}

# Ensure both datasets are ordered consistently by sample ID
# Critical for proper alignment between unscaled and scaled data
cluster_data <- cluster_data[order(as.numeric(rownames(cluster_data))), ]
cluster_data_scale <- cluster_data_scale[order(as.numeric(rownames(cluster_data_scale))), ]

# Run PAM clustering on scaled data with fixed parameters
set.seed(24)                    # Reproducible results
k <- 6                          # Number of desired clusters
pam_res <- pam(
  cluster_data_scale,           # Scaled metabolite/clinical features
  k = k,                        # 6 clusters
  nstart = 100,                 # Multiple random starts for stability
  metric = "euclidean"          # Distance metric
)


```

### Calinski-Harabasz Index

```{r calinski_h_factor}
cal_scale <- calinhara(cluster_data, pam_res$clustering, cn = k)

#Calinski-Harabasz index :
cat("The Cazlinski Score ", round(cal_scale))

# The Cazlinski Score  296
```

### VISUALIZATIONS OF THE CLUSTERS

```{r cluster_visualization processing}
# Extract cluster assignments from PAM (Partitioning Around Medoids) results
cluster <- data.frame(cluster = pam_res$clustering)


# Merge clustering results with existing cluster_data (by row names)
cluster_data <- merge(cluster_data, cluster, by = 0)
cluster_data <- cluster_data[order(as.numeric(rownames(cluster_data))), ]
# Restore original row names and clean up the Row.names column
row.names(cluster_data) <- cluster_data$Row.names
cluster_data$Row.names <- NULL

# Define consistent color palette for all visualizations (10 colors for flexibility)
def_colours <- c("forestgreen", "red", "darkblue", "orange", 
                 "cyan", "brown", "purple", "grey60", "gold", "green")

```

## PIE CHART 
### For the 6 clusters

```{r pie and radar chart for 6 clusters}
# Generate frequency table for the 6-cluster solution
freq_ABOS_6Clust <- data.frame(table(cluster_data$cluster))

# Create pie chart showing distribution across all 6 clusters
pieFct(
  dataFreq = freq_ABOS_6Clust,
  col = def_colours,        # Full 10-color palette (uses first 6)
  persLegend = TRUE,        # Include percentage labels in legend
  cohort = "ABOS cohort",
  labelSize = 4
)

# Generate radar plot comparing features across 6 clusters
radarInFct(
  dataName = "cluster_data",
  varClustName = "cluster", # Use full 6-cluster assignments
  color = def_colours[1:6]  # First 6 colors from palette
)

```

### For the 3 clusters

```{r Pie chart and radar plot for the 3 clusters}
# Create simplified 3-class grouping from 6-cluster results
# Maps clusters 2→CM (Cadiometabolic), 5→LS (Liver Specific), others→CTRL
cluster_data$class <- ifelse(cluster_data$cluster == 2, "CM",
                             ifelse(cluster_data$cluster == 5, "LS", "CTRL"))

# Generate frequency table for the 3-class grouping
freq_ABOS_3Clust <- data.frame(table(cluster_data$class))

# Define desired order for plotting (CM > LS > CTRL)
levels_3clust <- c("CM", "LS", "CTRL")

# Convert to factor with proper ordering
freq_ABOS_3Clust$Var1 <- factor(freq_ABOS_3Clust$Var1, levels = levels_3clust)

# Create pie chart showing class distribution in ABOS cohort
pieFct3cl(
  dataFreq = freq_ABOS_3Clust,
  col = def_colours[c(2,5,8)],      # Colors matching previous plots
  cohort = "ABOS cohort",
  labelSize = 4
)

# Generate radar plot showing metabolite/clinical features by 3-class
radarInFct(
  dataName = "cluster_data",
  varClustName = "class",           # Use the new 3-class variable
  color = def_colours[c(2,5,8)],
  clustOrder = levels_3clust
)

# Save final annotated clustering results if needed
#write.csv(cluster_data, "cluster_data_mash.csv", row.names = TRUE)

```

### MASH

```{r}
# Add MASH clinical phenotype to the clustering results
# Matches by row names to ensure proper sample alignment
cluster_data$Mash <- clinical_data$Mash[match(rownames(cluster_data),rownames(clinical_data))]

# Save intermediate results if needed
# write.xlsx(cluster_data, "cluster_data_mash.xlsx", row.names = TRUE)

# Remove any samples with missing MASH phenotype data
cluster_dataM <- cluster_data[!is.na(cluster_data$Mash), ]
nrow(cluster_dataM)

# Create bar plot for 6-cluster analysis
# Tests association between clusters and MASH phenotype
barPlotTestsFct(
  data = cluster_dataM,
  clust = "cluster",           # Column name for 6-cluster assignments
  outcome = "Mash",            # Binary outcome (MASH presence)
  position = c(0.5, 0, 0, 0.87, 0, 0.57, 0.67, 0, 0.95, 0, 0.77, 0, 0.5, 0, 0.57),
  color = def_colours[1:6],    # Colors for 6 clusters
  ymax = 1,
  labY = "MASH (%)"
)

# Check MASH distribution in the cleaned dataset
table(cluster_dataM$Mash)

# Create bar plot for 3-cluster analysis
barPlotTestsFct(
  data = cluster_dataM,
  clust = "class",             # Column name for 3-cluster assignments
  outcome = "Mash",
  position = c(0.75, 0, 0.55), # P-values/brackets positions for 3 clusters
  color = def_colours[c(2,5,8)], # Specific colors for 3 clusters
  ymax = 1,
  labY = "MASH (%)",
  clustOrder = levels_3clust    # Custom cluster ordering
)

nrow(cluster_dataM)

```

## =============================================================================
## DEA PLASMA METABOLITES 
## =============================================================================

### Preprocessing 
```{r Metabolomics }
# Check how many patients are matched between the expression and sample info sheets
cat("Matched", length(common_rows), "patients out of", nrow(expression_sheet), "\n")

# Subset sample_info so its rows correspond to the same order as expression_sheet
mah_subset <- sample_info[match(expression_sheet$PARENT_SAMPLE_NAME,
                                sample_info$PARENT_SAMPLE_NAME), ]

# Add the CLIENT_SAMPLE_ID (renamed acommon_samples <- intersect(colnames(liver_metabolite_DA), rownames(cluster_dataM))s Abos_id) to expression_sheet
# Backticks are used because the original column name contains spaces
expression_sheet$Abos_id <- mah_subset$CLIENT_SAMPLE_ID

#class(expression_sheet) Identify and remove any rows with missing IDs (unmatched samples)
sum(is.na(expression_sheet$Abos_id))
expression_sheet <- na.omit(expression_sheet)

# Convert the cleaned expression data into a proper dataframe
metabolite_data <- as.data.frame(expression_sheet)

# Assign unique row names using the sample identifier
rownames(metabolite_data) <- metabolite_data$Abos_id

# Remove columns no longer needed
metabolite_data$PARENT_SAMPLE_NAME <- NULL
metabolite_data$Abos_id <- NULL

# Order rows numerically by their sample IDs for consistency
metabolite_data <- metabolite_data[order(as.numeric(rownames(metabolite_data))), ]

# Match and rename columns (metabolites) using chemical_details
match_idx <- match(colnames(metabolite_data), chemical_details$CHEM_ID)
new_colnames <- chemical_details$COMP_ID[match_idx]

# Replace column names where a match exists
colnames(metabolite_data) <- ifelse(!is.na(match_idx),
                                    new_colnames,
                                    colnames(metabolite_data))

# Convert to a matrix and transpose it so metabolites are rows and samples are columns
metabolite_data <- as.matrix(metabolite_data)
plasma_matrix <- t(metabolite_data)

# Filter out Xenobiotics (non-endogenous compounds) from the chemical details
chemical_details <- chemical_details %>%
  filter(SUPER_PATHWAY != "Xenobiotics")

# Alternatively:
# chemical_details <- chemical_details[!is.na(chemical_details$SUPER_PATHWAY) &
#                                      chemical_details$SUPER_PATHWAY != "Xenobiotics", ]

# Get the list of valid metabolite IDs to keep
valid_ids <- chemical_details$COMP_ID

# Filter the plasma matrix to retain only valid metabolites
plasma_matrix <- plasma_matrix[rownames(plasma_matrix) %in% valid_ids, ]

# Save results if needed
# write.xlsx(plasma_matrix, "plasma_metabolite_matrix.xlsx", Row.names = TRUE)

# Summarize the final data dimensions
cat("Number of patients that have metabolite data: ", ncol(plasma_matrix), "\n")
cat("Number of metabolites used for the analysis: ", nrow(plasma_matrix))


# Matched 1545 patients out of 1558 
# [1] 14
# Number of patients that have metabolite data:  1544 
# Number of metabolites used for the analysis:  875

```


```{r}
# 1. Get the column names you want to match
cluster_metabo_name <- colnames(plasma_matrix)

# 2. Subset the cluster vector 
cluster_metaboplasm  <-setNames(cluster[colnames(plasma_matrix),],
                              colnames(plasma_matrix))

cluster_metaboplasm  <- na.omit(cluster_metaboplasm)

# 4. Check length before cleaning
cat("Length before cleaning:", length(cluster_metaboplasm), "\n")

# 5. Check final length
cat("Total number of patients used to Diffrential Analysis:", length(cluster_metaboplasm), "\n")

#6. Data to be used for Differential Analysis 
plasma_metabolite_DA <- plasma_matrix[,names(cluster_metaboplasm )]
identical(colnames(plasma_metabolite_DA), names(cluster_metaboplasm ))


cluster_metaboplasm <- ifelse(cluster_metaboplasm  == 2,"CM",
                             ifelse(cluster_metaboplasm  == 5,"LS", "CTRL"))


# 7. Check final length of the Data to be used for clustering
cat("Number of patients in the Data frame:", length(cluster_metaboplasm ), "\n")
cat("Number of metabolites in the Data frame:", nrow(plasma_metabolite_DA ), "\n")

table(cluster_metaboplasm)

# Length before cleaning: 1321 
# Total number of patients used to Diffrential Analysis: 1321 
# [1] TRUE
# Number of patients in the Data frame: 1321 
# Number of metabolites in the Data frame: 875 
# cluster_metaboplasm
#   CM CTRL   LS 
#  151 1075   95 
```

### Analysis

```{r}

map_chemical_names <- function(mstat, chemical_details) {
  indices <- match(row.names(mstat), chemical_details$COMP_ID)
  chemical_names <- chemical_details$CHEMICAL_NAME[indices]
  
  # Update row names safely (fallback to ID if name missing)
  row.names(mstat) <- ifelse(
    is.na(chemical_names) | chemical_names == "", 
    row.names(mstat), 
    chemical_names
  )
  return(mstat)
}
#  SETUP DESIGN & FIT
design <- model.matrix(~ -1 + cluster_metaboplasm)
fit <- lmFit(plasma_metabolite_DA, design)
# CM vs CTRL COMPARISON 
comp_CMvsCTRL <- comparisonsLimmaFct(fit, "cluster_metaboplasmCM - cluster_metaboplasmCTRL", design,  nrow(plasma_metabolite_DA))
comp_CMvsCTRL$mstat <- map_chemical_names(comp_CMvsCTRL$mstat, chemical_details)
results_CMvsCTRL <- printResultsLimmaFct(comp_CMvsCTRL, topPrintHist = FALSE,topPrintVolc = FALSE, color = c("red", "red", "grey70"), legendPos = "none")
down_CMvsCTRL <- nrow(results_CMvsCTRL$results[results_CMvsCTRL$results$logFC < 0,])
up_CMvsCTRL <- nrow(results_CMvsCTRL$results[results_CMvsCTRL$results$logFC >0,])
cat("Number of down regulated metabolites CMvsCTRL :", down_CMvsCTRL, "metabolites\n")
cat("Number of up regulated metabolites CMvsCTRL :", up_CMvsCTRL, "metabolites\n")

# --- 3. LS vs CTRL COMPARISON ---
comp_LSvsCTRL <- comparisonsLimmaFct(fit, "cluster_metaboplasmLS - cluster_metaboplasmCTRL", design, nbGenes = nrow(plasma_metabolite_DA))
comp_LSvsCTRL$mstat <- map_chemical_names(comp_LSvsCTRL$mstat, chemical_details)
results_LSvsCTRL <- printResultsLimmaFct(comp_LSvsCTRL, topPrintHist = FALSE,topPrintVolc = FALSE, color = c("red", "red", "grey70"), legendPos = "none")
down_LSvsCTRL<- nrow(results_LSvsCTRL$results[results_LSvsCTRL$results$logFC < 0,])
up_LSvsCTRL<- nrow(results_LSvsCTRL$results[results_LSvsCTRL$results$logFC >0,])
cat("Number of down regulated metabolites LSvsCTRL :", down_LSvsCTRL, "metabolites\n")
cat("Number of up regulated metabolites LSvsCTRL :", up_LSvsCTRL, "metabolites\n")


# --- 4. CM vs LS COMPARISON ---
comp_CMvsLS <- comparisonsLimmaFct(fit, "cluster_metaboplasmCM - cluster_metaboplasmLS", design, nbGenes = nrow(plasma_metabolite_DA))
comp_CMvsLS$mstat <- map_chemical_names(comp_CMvsLS$mstat, chemical_details)
results_CMvsLS <- printResultsLimmaFct(comp_CMvsLS, topPrintHist = FALSE,topPrintVolc = FALSE, color = c("red", "red", "grey70"), legendPos = "none")
down_CMvsLS <- nrow(results_CMvsLS$results[results_CMvsLS$results$logFC < 0,])
up_CMvsLS <-nrow(results_CMvsLS$results[results_CMvsLS$results$logFC >0,])
cat("Number of down regulated metabolites CMvsLS :", down_CMvsLS, "metabolites\n")
cat("Number of up regulated metabolites CMvsLS :", up_CMvsLS, "metabolites\n")




# There are 546 differentially expressed features considering threshold 0.05 on adjusted P-values, and 240 when adding the threshold 0.26 on absolute value of log2 Fold Change. 
# Number of down regulated metabolites CMvsCTRL : 41 metabolites
# Number of up regulated metabolites CMvsCTRL : 199 metabolites
# There are 123 differentially expressed features considering threshold 0.05 on adjusted P-values, and 39 when adding the threshold 0.26 on absolute value of log2 Fold Change. 
# Number of down regulated metabolites LSvsCTRL : 2 metabolites
# Number of up regulated metabolites LSvsCTRL : 37 metabolites
# There are 380 differentially expressed features considering threshold 0.05 on adjusted P-values, and 188 when adding the threshold 0.26 on absolute value of log2 Fold Change. 
# Number of down regulated metabolites CMvsLS : 37 metabolites
# Number of up regulated metabolites CMvsLS : 151 metabolites
```

### Visualizations

```{r}
# View tables and plots
results_CMvsCTRL$volcano
results_LSvsCTRL$volcano
results_CMvsLS$volcano


# Example datatable sorted by logFC
datatable(results_CMvsCTRL$results[order(results_CMvsCTRL$results$logFC, decreasing = TRUE),])
datatable(results_LSvsCTRL$results[order(results_LSvsCTRL$results$logFC, decreasing = TRUE),])
datatable(results_CMvsLS$results[order(results_CMvsLS$results$logFC, decreasing = TRUE),])

# Create Euler Diagram with Chemical Names
DEGs <- list(
  "CM versus LS" = rownames(results_CMvsLS$results),
  "CM versus CTRL" = rownames(results_CMvsCTRL$results),
  "LS versus CTRL" = rownames(results_LSvsCTRL$results)
)
# Your code + 1-line fix
DEGs <- euler(DEGs, 
              loss = "abs", 
              loss_aggregator = "max",
              shape = "circle")

# Plot (fills length matches your sets)
plot(DEGs, quantities = TRUE, fills = c("blueviolet", "red", "cyan"), alpha = 0.6)


# BAR PLOTS

create_top_plot(results_CMvsCTRL, "CM vs CTRL")
create_top_plot(results_LSvsCTRL, "LS vs CTRL")
create_top_plot(results_CMvsLS, "CM vs LS")



```

```{r}
genes<-df
genes2 <- genes
rownames(genes) <- genes$X    # Assign from column
genes$Id <- NULL 
genes$X <- NULL
genes$classeTot <- NULL# Remove column safely
#rownames(genes)  
genes <- t(genes)
```


```{r}
# 1. Get clusters matching matrix columns (FIX: add comma for columns)
cluster_gene_name <- colnames(genes)
cluster_transcript <- setNames(cluster[cluster_gene_name,], cluster_gene_name)

# Rest unchanged...
cluster_transcript <- cluster_transcript[!is.na(cluster_transcript)]
cat("Total patients for DA:", length(cluster_transcript), "\n")

plasma_transcript_DA <- genes[, names(cluster_transcript), drop = FALSE]
plasma_transcript_DA <- as.matrix(plasma_transcript_DA)
storage.mode(plasma_transcript_DA) <- "double"


# Verify (should be TRUE now)
all(colnames(plasma_transcript_DA) == names(cluster_transcript))
identical(colnames(plasma_transcript_DA), names(cluster_transcript))
sum(colnames(plasma_transcript_DA) != names(cluster_transcript))  # 0

# Recode
cluster_transcript <- ifelse(cluster_transcript == 2, "CM",
                             ifelse(cluster_transcript == 5, "LS",
                                    ifelse(cluster_transcript %in% c(1,3,4,6), "CTRL", "Pb")))

cat("Patients in DA:", length(cluster_transcript), "\n")
table(cluster_transcript)


ngenes<- nrow(plasma_transcript_DA)
cat("Total numebr of genes used " ,ngenes, "\n")

```

```{r}
design1 <- model.matrix(~ -1 + cluster_transcript)
row.names(design1) <-names(cluster_transcript)
colnames(design1)
fit <- lmFit(plasma_transcript_DA, design1)         # Fit linear model per gene


# CM vs CTRL COMPARISON 
comp1_CMvsCTRL <- comparisonsLimmaFct(fit, "cluster_transcriptCM - cluster_transcriptCTRL", design1,  ngenes)
results1_CMvsCTRL <- printResultsLimmaFct(comp1_CMvsCTRL, topPrintHist = FALSE,topPrintVolc = FALSE)
down_CMvsCTRL <- nrow(results1_CMvsCTRL$results[results1_CMvsCTRL$results$logFC < 0,])
up_CMvsCTRL <- nrow(results1_CMvsCTRL$results[results1_CMvsCTRL$results$logFC >0,])
cat("Number of down regulated genes CMvsCTRL :", down_CMvsCTRL, "genes\n")
cat("Number of up regulated genes CMvsCTRL :", up_CMvsCTRL, "genes\n")



# LS vs CTRL COMPARISON 
comp1_LSvsCTRL <- comparisonsLimmaFct(fit, "cluster_transcriptLS - cluster_transcriptCTRL", design1,  ngenes)
results1_LSvsCTRL <- printResultsLimmaFct(comp1_LSvsCTRL, topPrintHist = FALSE,topPrintVolc = FALSE)
down_LSvsCTRL <- nrow(results1_LSvsCTRL$results[results1_LSvsCTRL$results$logFC < 0,])
up_LSvsCTRL <- nrow(results1_LSvsCTRL$results[results1_LSvsCTRL$results$logFC >0,])
cat("Number of down regulated genes LSvsCTRL :", down_LSvsCTRL, "genes\n")
cat("Number of up regulated genes LSvsCTRL :", up_LSvsCTRL, "genes\n")


# --- 4. CM vs LS COMPARISON ---
comp1_CMvsLS <- comparisonsLimmaFct(fit, "cluster_transcriptCM - cluster_transcriptLS", design1, ngenes)
results1_CMvsLS <- printResultsLimmaFct(comp1_CMvsLS, topPrintHist = FALSE,topPrintVolc = FALSE)
down_CMvsLS <- nrow(results1_CMvsLS$results[results1_CMvsLS$results$logFC < 0,])
up_CMvsLS <-nrow(results1_CMvsLS$results[results1_CMvsLS$results$logFC >0,])
cat("Number of down regulated genes CMvsLS :", down_CMvsLS, "genes\n")
cat("Number of up regulated genes CMvsLS :", up_CMvsLS, "genes\n")


results1_CMvsCTRL$volcano
results1_LSvsCTRL$volcano
results1_CMvsLS$volcano

datatable(results1_CMvsCTRL$results[order(results1_CMvsCTRL$results$logFC, decreasing = TRUE),])
datatable(results1_LSvsCTRL$results[order(results1_LSvsCTRL$results$logFC, decreasing = TRUE),])
datatable(results1_CMvsLS$results[order(results1_CMvsLS$results$logFC, decreasing = TRUE),])

# There are 482 differentially expressed features considering threshold 0.05 on adjusted P-values, and 371 when adding the threshold 0.26 on absolute value of log2 Fold Change. 
# Number of down regulated genes CMvsCTRL : 189 genes
# Number of up regulated genes CMvsCTRL : 182 genes
# There are 59 differentially expressed features considering threshold 0.05 on adjusted P-values, and 57 when adding the threshold 0.26 on absolute value of log2 Fold Change. 
# Number of down regulated genes LSvsCTRL : 15 genes
# Number of up regulated genes LSvsCTRL : 42 genes
# There are 35 differentially expressed features considering threshold 0.05 on adjusted P-values, and 35 when adding the threshold 0.26 on absolute value of log2 Fold Change. 
# Number of down regulated genes CMvsLS : 10 genes
# Number of up regulated genes CMvsLS : 25 genes
```




## =============================================================================
## DEA LIVER METABOLITES 
## =============================================================================
```{r}
# Filter out Xenobiotics pathway metabolites
sample_info_liver <- sample_info_liver %>% filter(sample_info_liver$SUPER.PATHWAY != "Xenobiotics")

# Keep only matching rows between data and sample info
keep_rows <- sample_info_liver$COMP.ID

liver_matrix <- liver_data[rownames(liver_data) %in% keep_rows, , drop = FALSE]

```

```{r}
# Align cluster data with metabolite matrix columns
cluster_liver_metaboplasm <- setNames(cluster[colnames(liver_matrix), ], colnames(liver_matrix))
cluster_liver_metaboplasm <- na.omit(cluster_liver_metaboplasm)  # Remove unmatched samples

cat("Length before cleaning:", length(cluster_liver_metaboplasm), "\n")

# Filter matrix to cluster-matched samples only
liver_metabolite_DA <- as.data.frame(liver_matrix[, names(cluster_liver_metaboplasm)])

# Convert numeric clusters to labels: 2=CM, 5=LS, others=CTRL
cluster_liver_metaboplasm <- ifelse(cluster_liver_metaboplasm == 2, "CM",
                                    ifelse(cluster_liver_metaboplasm == 5, "LS", "CTRL"))

cat("Number of liver metabolite samples for analysis:", ncol(liver_metabolite_DA), "\n")
cat("Number of liver metabolites:", nrow(liver_metabolite_DA), "\n")



# Length before cleaning: 295 
# Number of liver metabolite samples for analysis: 295 
# Number of liver metabolites: 817 
```

### Analysis

```{r}
map_metabolite_names <- function(mstat_object, sample_info){ 
  # Map row names to chemical names safely
  sample_info$COMP.ID <- as.character(sample_info$COMP.ID)
  indices <- match(row.names(mstat_object), sample_info$COMP.ID)
  chemical_names <- sample_info$BIOCHEMICAL[indices]
  
  # Replace row names (fallback to ID if missing)
  row.names(mstat_object) <- ifelse(is.na(chemical_names) | chemical_names == "", 
                                    row.names(mstat_object), chemical_names)
  return(mstat_object)
}

# --- 1. SETUP ---
 design <- model.matrix(~ -1 + cluster_liver_metaboplasm)
fit <- lmFit(liver_metabolite_DA, design)
# --- 2. CM vs CTRL ---
comp_CMvsCTRL <- comparisonsLimmaFct(fit, "cluster_liver_metaboplasmCM - cluster_liver_metaboplasmCTRL", design, nrow(liver_metabolite_DA))
comp_CMvsCTRL$mstat <- map_metabolite_names(comp_CMvsCTRL$mstat, sample_info_liver)
res_CMvsCTRL <- printResultsLimmaFct(comp_CMvsCTRL, topPrintHist = FALSE, topPrintVolc = FALSE)

down_CMvsCTRL <- nrow(res_CMvsCTRL$results[res_CMvsCTRL$results$logFC < 0,])
up_CMvsCTRL <- nrow(res_CMvsCTRL$results[res_CMvsCTRL$results$logFC >0,])
cat("Number of down regulated metabolites CMvsCTRL :", down_CMvsCTRL, "metabolites\n")
cat("Number of up regulated metabolites CMvsCTRL :", up_CMvsCTRL, "metabolites\n")

# --- 3. LS vs CTRL ---
comp_LSvsCTRL <- comparisonsLimmaFct(fit, "cluster_liver_metaboplasmLS - cluster_liver_metaboplasmCTRL", design, nrow(liver_metabolite_DA))
comp_LSvsCTRL$mstat <- map_metabolite_names(comp_LSvsCTRL$mstat, sample_info_liver)
res_LSvsCTRL <- printResultsLimmaFct(comp_LSvsCTRL, topPrintHist = FALSE, topPrintVolc = FALSE)

down_LSvsCTRL<- nrow(res_LSvsCTRL$results[res_LSvsCTRL$results$logFC < 0,])
up_LSvsCTRL<- nrow(res_LSvsCTRL$results[res_LSvsCTRL$results$logFC >0,])
cat("Number of down regulated metabolites LSvsCTRL :", down_LSvsCTRL, "metabolites\n")
cat("Number of up regulated metabolites LSvsCTRL :", up_LSvsCTRL, "metabolites\n")

# --- 4. CM vs LS ---
comp_CMvsLS <- comparisonsLimmaFct(fit, "cluster_liver_metaboplasmCM - cluster_liver_metaboplasmLS", design, nrow(liver_metabolite_DA))
comp_CMvsLS$mstat <- map_metabolite_names(comp_CMvsLS$mstat, sample_info_liver)
res_CMvsLS <- printResultsLimmaFct(comp_CMvsLS, topPrintHist = FALSE, topPrintVolc = FALSE)

down_CMvsLS <- nrow(res_CMvsLS $results[res_CMvsLS $results$logFC < 0,])
up_CMvsLS <-nrow(res_CMvsLS $results[res_CMvsLS $results$logFC >0,])
cat("Number of down regulated metabolites CMvsLS :", down_CMvsLS, "metabolites\n")
cat("Number of up regulated metabolites CMvsLS :", up_CMvsLS, "metabolites\n")


# There are 184 differentially expressed features considering threshold 0.05 on adjusted P-values, and 178 when adding the threshold 0.26 on absolute value of log2 Fold Change. 
# Number of down regulated metabolites CMvsCTRL : 89 metabolites
# Number of up regulated metabolites CMvsCTRL : 89 metabolites
# There are 122 differentially expressed features considering threshold 0.05 on adjusted P-values, and 120 when adding the threshold 0.26 on absolute value of log2 Fold Change. 
# Number of down regulated metabolites LSvsCTRL : 49 metabolites
# Number of up regulated metabolites LSvsCTRL : 71 metabolites
# There are 24 differentially expressed features considering threshold 0.05 on adjusted P-values, and 24 when adding the threshold 0.26 on absolute value of log2 Fold Change. 
# Number of down regulated metabolites CMvsLS : 12 metabolites
# Number of up regulated metabolites CMvsLS : 12 metabolites
```

### Visualizations

```{r}
 res_CMvsCTRL$volcano
res_LSvsCTRL$volcano
res_CMvsLS$volcano
# ggsave("volcano_CMvsCTRL.pdf",plot = res_CMvsCTRL$volcano,width = 10, height = 8, units = "in")

# ggsave("volcano_LSvsCTRL.pdf",plot = res_LSvsCTRL$volcano, width = 10, height = 8, units = "in")

# ggsave("volcano_CMvsLS.pdf", plot = res_CMvsLS$volcano, width = 10, height = 8, units = "in")


# Tables
datatable(res_CMvsCTRL$results[order(res_CMvsCTRL$results$logFC, decreasing = TRUE),])
datatable(res_LSvsCTRL$results[order(res_LSvsCTRL$results$logFC, decreasing = TRUE),])
datatable(res_CMvsLS$results[order(res_CMvsLS$results$logFC, decreasing = TRUE),])


#png("top10_CMvsCTRL.png", width=11, height=8, units="in", res=300, bg="white")
create_top_plot(res_CMvsCTRL, "CM vs CTRL")
#dev.off()

#png("top10_LSvsCTRL.png", width=11, height=8, units="in", res=300, bg="white")
create_top_plot(res_LSvsCTRL, "LS vs CTRL")
#dev.off()

#png("top10_CMvsLS.png", width=11, height=8, units="in", res=300, bg="white")
create_top_plot(res_CMvsLS, "CM vs LS")
#dev.off()


#Euler plot h smaller text
euler_data <- list(
  "CM versus CTRL" = rownames(res_CMvsCTRL$results),
  "LS versus CTRL" = rownames(res_LSvsCTRL$results),
  "CM versus LS"   = rownames(res_CMvsLS$results)
)

deg <- euler(euler_data, 
             loss = "abs", 
             loss_aggregator = "max",
             shape = "circle")
#png("enhanced_deg_plot", width=11, height=8, units="in", res=300, bg="white")
plot(deg, 
     quantities = TRUE, fills = c("blueviolet", "red", "cyan"), alpha = 0.6)
#dev.off()


# write.csv(res_CMvsCTRL$results[order(res_CMvsCTRL$results$logFC, decreasing = TRUE), ],
#           "res_liver_CMvsCTRL_sorted.csv", row.names = TRUE)
# 
# write.csv(res_LSvsCTRL$results[order(res_LSvsCTRL$results$logFC, decreasing = TRUE), ],
#           "res_liver_LSvsCTRL_sorted.csv", row.names = TRUE)
# 
#  write.csv(res_CMvsLS$results[order(res_CMvsLS$results$logFC, decreasing = TRUE), ],
#           "res_liver_CMvsLS_sorted.csv", row.names = TRUE)



```

```{r}
A <- unique(euler_data[["CM versus CTRL"]])
B <- unique(euler_data[["LS versus CTRL"]])
C <- unique(euler_data[["CM versus LS"]])

only_C <- setdiff(C, union(A, B))

cat("Metabolites that are only in CM vs LS compartment:\n")
print(only_C)
print(res_CMvsLS$results[only_C,])
```

## =============================================================================
## ENRICHMENT ANALYSIS 
## =============================================================================

## PLASMA
```{r}
# Map the plasma metabolite lists to HMDB identifiers
LSvsCTRL_enrichp <- map_to_hmdb(rownames(results_LSvsCTRL$results), chemical_details, name_col = "CHEMICAL_NAME")
CMvsCTRL_enrichp <- map_to_hmdb(rownames(results_CMvsCTRL$results), chemical_details, name_col = "CHEMICAL_NAME")
CMvsLS_enrichp   <- map_to_hmdb(rownames(results_CMvsLS$results), chemical_details, name_col = "CHEMICAL_NAME")

# Run pathway enrichment analysis for each plasma comparison
LSvsCTRL_resultsp <- run_pathway_enrichment(LSvsCTRL_enrichp)
CMvsCTRL_resultsp <- run_pathway_enrichment(CMvsCTRL_enrichp)
CMvsLS_resultsp   <- run_pathway_enrichment(CMvsLS_enrichp)


cols_wanted <- c(
  "pathwayName",
  "pathwayId",
  "pathwaySource",
  "Num_In_Path",
  "Total_In_Path",
  "Pval_FDR",
  "analytes"
)

# Select only these columns for each table
CMvsCTRL_resultsp <- CMvsCTRL_resultsp[, cols_wanted]
LSvsCTRL_resultsp <- LSvsCTRL_resultsp[, cols_wanted]
CMvsLS_resultsp   <- CMvsLS_resultsp[, cols_wanted]

# Display the enrichment results
datatable(CMvsCTRL_resultsp)
datatable(LSvsCTRL_resultsp)
datatable(CMvsLS_resultsp)

# Generate bar plots for the enriched pathways in each plasma comparison
plot_pathway_enrichment(CMvsCTRL_resultsp, "PLASMA: CM vs CTRL")
plot_pathway_enrichment(LSvsCTRL_resultsp, "PLASMA: LS vs CTRL")
plot_pathway_enrichment(CMvsLS_resultsp, "PLASMA: CM vs LS")
```

## LIVER

```{r}
# Map the metabolite lists to HMDB identifiers for the liver comparisons
LSvsCTRL_enrich <- map_to_hmdb(rownames(res_LSvsCTRL$results), sample_info_liver)
CMvsCTRL_enrich <- map_to_hmdb(rownames(res_CMvsCTRL$results), sample_info_liver)
CMvsLS_enrich   <- map_to_hmdb(rownames(res_CMvsLS$results), sample_info_liver)

# Run pathway enrichment analysis for each liver comparison
LSvsCTRL_results <- run_pathway_enrichment(LSvsCTRL_enrich)
CMvsCTRL_results <- run_pathway_enrichment(CMvsCTRL_enrich)
CMvsLS_results   <- run_pathway_enrichment(CMvsLS_enrich)

# Define the columns you want, in this exact order
cols_wanted <- c(
  "pathwayName",
  "pathwayId",
  "pathwaySource",
  "Num_In_Path",
  "Total_In_Path",
  "Pval_FDR",
  "analytes"
)

# Select only these columns for each table
CMvsCTRL_results <- CMvsCTRL_results[, cols_wanted]
LSvsCTRL_results <- LSvsCTRL_results[, cols_wanted]
CMvsLS_results   <- CMvsLS_results[, cols_wanted]

# View the results
datatable(CMvsCTRL_results)
datatable(LSvsCTRL_results)
datatable(CMvsLS_results)



# Generate bar plots for the enriched pathways in each liver comparison
p1 <- plot_pathway_enrichment(CMvsCTRL_results, "LIVER: CM vs CTRL")
p2 <- plot_pathway_enrichment(LSvsCTRL_results, "LIVER: LS vs CTRL")
p3 <- plot_pathway_enrichment(CMvsLS_results, "LIVER: CM vs LS")

# ggsave("pathway_enrichment_CMvsCTRL.pdf", plot = p1, width = 10, height = 8, units = "in")
# ggsave("pathway_enrichment_LSvsCTRL.pdf", plot = p2, width = 10, height = 8, units = "in")
# ggsave("pathway_enrichment_CMvsLS.pdf", plot = p3, width = 10, height = 8, units = "in")

```


## Clusters Validation
```{r}
# Find common samples between liver metabolite data & clinical clusters
common_samples <- intersect(colnames(liver_matrix), rownames(cluster_dataM))

#  Filter cluster_data to EXACTLY matching samples 
cluster_aligned <- cluster_dataM[common_samples, , drop = FALSE]

# Save ALL THREE cluster analysis plots as thesis-ready PNGs
freq_3Clust <- data.frame(table(cluster_aligned$class))

# 3. Pie chart - 3-cluster frequencies
#png("pie_3cluster_frequencies.png", width=10, height=9, units="in", res=300, bg="white")
pieFct3cl(dataFreq = freq_3Clust,
          col= def_colours[c(2,8,5)],  
          cohort = "ABOS patients with liver metabolite", 
          labelSize = 4)
#dev.off()

# Load ggplot2 if not already loaded


# 4. Radar plot - Clinical variables comparison
#pdf(" Radar plot - Clinical variables 3cluster.pdf", width = 10, height = 8, bg = "white")
radarInFct(dataName = "cluster_aligned", 
                      varClustName = "class",           
                      color = def_colours[c(2,5,8)],   
                      clustOrder = levels_3clust)


#dev.off()



# 5. Bar plot - MASH prevalence per cluster

#pdf("mash_prevalence_3cluster.pdf", width = 10, height = 8, bg = "white")
 barPlotTestsFct(data = cluster_aligned,
                         clust = "class",              
                         outcome = "Mash",             
                         position = c(0.95, 0, 0.75),  
                         color = def_colours[c(2,5,8)],
                         ymax = 1,                     
                         labY = "MASH (%)",
                         clustOrder = levels_3clust)
#dev.off()


       

# Calculate MASH percentages per cluster
mash_table <- table(cluster_aligned$class, cluster_aligned$Mash)
mash_percent <- prop.table(mash_table, margin = 1) * 100

# View results
print(mash_table)        # Raw counts
print(round(mash_percent, 1)) 



cluster_aligned2 <- clinical_data[common_samples, , drop = FALSE]

# Absolute frequencies by sex
sex_table <- table(cluster_aligned2$sexe) 

# Percentages for the entire dataset (total = 100%)
sex_percent <- prop.table(sex_table) * 100

# Display results
print("Absolute frequencies:")
print(sex_table)
print("Overall percentages:")
print(round(sex_percent, 1))
```

## =============================================================================
## CORRELATION OF METABOLITES WITH HISTOLOGY VARIABLES
## =============================================================================
```{r}

# Create a data frame for histological scores from the clinical dataset
histo <- data.frame(
  row.names = rownames(clinical_data),
  Steatosis = as.numeric(clinical_data$M0_steatose_score),
  Inflammation = as.numeric(clinical_data$Inflammation),
  Ballooning = as.numeric(clinical_data$Ballooning),
  Mash = as.numeric(clinical_data$Mash),
  Fibrosis = as.numeric(clinical_data$Biopsie_Kleiner_Fibrose)
)

clinical_var <- data.frame(
  row.names = rownames(clinical_data),
  Age  = clinical_data$AgeJourIntervention,
  BMI  = clinical_data$BMI,
  SEX  = clinical_data$sexe
)

# Match cluster class by row names, leaving NA where no match
clinical_var$CLASS <- cluster_data$class[ match(
  rownames(clinical_var),
  rownames(cluster_data)
) ]


# Keep only rows whose rownames exist in cluster_data, in the same order
common.ids <- intersect( rownames(histo), rownames(cluster_data) )
histo      <- histo[ common.ids, , drop = FALSE ]
# Remove rows with any missing values from both tables
clinical_var <- na.omit(clinical_var)
histo <- na.omit(histo)

```

## PLASMA 

```{r}
# Find common samples between histological data and plasma metabolite data
common_stage_samples_plasma <- intersect(rownames(histo), colnames(plasma_metabolite_DA))
cat("Plasma samples:", length(common_stage_samples_plasma), "\n")

# Map metabolite annotations for plasma samples that match histology data
plasma_stages <- map_plasma_bio(plasma_metabolite_DA[, common_stage_samples_plasma, drop = FALSE], chemical_details)

# Transpose histology data: histological scores (rows) × samples (columns)
stages_t_plasma <- t(histo[common_stage_samples_plasma, , drop = FALSE])

# Calculate correlations between plasma metabolites and histological scores
plasma_results <- cor_metab_histology(plasma_stages, stages_t_plasma, "Plasma")

# Generate heatmap visualization of plasma metabolite-histology correlations

# Example call
ht <- make_stage_heatmap(
  sig_results = plasma_results,        # or liver_results
  name        = "Plasma vs Histology",
  histo_cols  = c("Steatosis", "Inflammation", "Ballooning", "Mash", "Fibrosis"),
  priority_mets =NULL,
  n_total     = 50,
  fdr_cutoff  = 0.05
)

draw(ht)

# How many cells are above fdr_cutoff?

significant_plasma <-plasma_results[plasma_results$Significant == "TRUE",]

#table(significant_plasma$Histology_St)
#write.csv(significant_plasma, "correlation_liver.csv", row.names = FALSE)



# Plasma samples: 1172 
# ===  Plasma ===
# Input metab dims: 875 1172 
# Input histo dims: 5 1172 
# Common samples: 1172 
# Sample names match: 47 48 49 50 53 54 
# Scaled metab dims: 875 1172 
# Processing stage: Steatosis 
# Processing stage: Inflammation 
# Processing stage: Ballooning 
# Processing stage: Mash 
# Processing stage: Fibrosis 
# Plasma sig metabolite-stage pairs (|ρ|>0.3 FDR<0.05): 2032 
# 
# 
#   Ballooning     Fibrosis Inflammation         Mash    Steatosis 
#          397          414          278          408          535 
```

## LIVER 

```{r}
# Find common samples between histological data and liver metabolite data
common_stage_samples_liver <- intersect(rownames(histo), colnames(liver_metabolite_DA))
cat("Liver samples:", length(common_stage_samples_liver), "\n")

# Map metabolite annotations for liver samples that match histology data
liver_stages <- map_liver_bio(liver_metabolite_DA[, common_stage_samples_liver, drop = FALSE], sample_info_liver)

# Transpose histology data: histological scores (rows) × samples (columns)
stages_t_liver <- t(histo[common_stage_samples_liver, , drop = FALSE])

# Calculate correlations between liver metabolites and histological scores
liver_results <- cor_metab_histology(liver_stages, stages_t_liver, "Liver")

# Use heatmap with FDR gating
ht <- make_stage_heatmap(
  sig_results     = liver_results,
  name            = NULL,
  histo_cols      = c("Steatosis", "Inflammation", "Ballooning", "Mash", "Fibrosis"),
  priority_mets   = priority_mets,
  n_total         = 30
)

#png("Liver vs Histology.png",width=11, height=10, units="in", res=300, bg="white")
draw(ht)
#dev.off()


significant_liver <-liver_results[liver_results$Significant == "TRUE",]

sig_mash <- significant_liver[
  significant_liver$Histology_Stages == "Mash" & 
    (significant_liver$rho <= -0. | significant_liver$rho >= 0.3),
]

sig_fibrosis<- significant_liver[
  significant_liver$Histology_Stages == "Fibrosis" & 
    (significant_liver$rho <= -0.1 | significant_liver$rho >= 0.1),
]


table(significant_liver$Histology_Stages)
#write.csv(sig_mash, "sig_mash.csv", row.names = FALSE)

# Liver samples: 287 
# ===  Liver ===
# Input metab dims: 817 287 
# Input histo dims: 5 287 
# Common samples: 287 
# Sample names match: 70 76 108 109 128 133 
# Scaled metab dims: 817 287 
# Processing stage: Steatosis 
# Processing stage: Inflammation 
# Processing stage: Ballooning 
# Processing stage: Mash 
# Processing stage: Fibrosis 
# Liver sig metabolite-stage pairs (|ρ|>0.3 FDR<0.05): 1534 
# 
# 
#   Ballooning     Fibrosis Inflammation         Mash    Steatosis 
#          251          221          355          245          462 
```

## =============================================================================
## METABOLITES ASSOCIATED WITH MASH AND FIBROSIS 
## =============================================================================
## MASH
```{R}
# STEP 1 & 2: Your code
common_samples <- intersect(rownames(histo), colnames(liver_stages))
liver_lasso <- liver_stages[rownames(liver_stages), common_samples]
liver_lasso <- t(liver_lasso)

liver_data <- data.frame(row.names = rownames(liver_lasso))
liver_data$Mash <- histo[rownames(liver_lasso), "Mash"]

# Add metabolites
met_names_original <- colnames(liver_lasso)
for (i in 1:length(met_names_original)) {
  liver_data[[met_names_original[i]]] <- liver_lasso[, i]
}

# Add covariates
covar_rows <- rownames(liver_data)
liver_data$age <- clinical_var[covar_rows, "Age"]
liver_data$sex <- clinical_var[covar_rows, "SEX"]
liver_data$BMI <- clinical_var[covar_rows, "BMI"]

# Remove rows with any NA
liver_data <- na.omit(liver_data)

# Convert sex to binary: 0 = female/other, 1 = male
liver_data$sex <- as.character(liver_data$sex)
liver_data$sex_bin <- ifelse(
  liver_data$sex %in% c("M", "Male", "MALE", "m", "male", 1, "1"),
  1,
  0
)

cat("Mash distribution:\n")
print(table(liver_data$Mash))
cat("Final N =", nrow(liver_data), "\n")
cat("Metabolites =", length(met_names_original), "\n")


# Verify binary outcome AFTER na.omit
table(liver_data$Mash)  # Must show exactly 2 levels
if(length(unique(liver_data$Mash)) != 2) {
  stop("Mash must be binary after na.omit!")
}

met_cols_mash <- setdiff(colnames(liver_data),
                         c("Mash", "age", "sex", "sex_bin", "BMI"))


res_mash_A <- run_lasso_nested_cv_modelA(
  data = liver_data,
  outcome_bin_col = "Mash",
  metabolite_cols = met_cols_mash,
  alpha = 1,
  n_outer = 10,
  n_inner = 10
)

res_mash <- run_lasso_nested_cv(
  data            = liver_data,
  outcome_bin_col = "Mash",
  covariate_cols  = c("age", "sex_bin", "BMI"),
  metabolite_cols = met_cols_mash)
```



```{r}
plot(res_mash_A$final_cv)
abline(v = log(res_mash_A$final_cv$lambda.min), col = "blue", lty = 2)
abline(v = log(res_mash_A$final_cv$lambda.1se), col = "red", lty = 2)
legend("topright",
       legend = c("lambda.min", "lambda.1se"),
       col = c("blue", "red"),
       lty = 2,
       bty = "n")



plot(res_mash$final_cv)
abline(v = log(res_mash$final_cv$lambda.min), col = "blue", lty = 2)
abline(v = log(res_mash$final_cv$lambda.1se), col = "red", lty = 2)
legend("topright",
       legend = c("lambda.min", "lambda.1se"),
       col = c("blue", "red"),
       lty = 2,
       bty = "n")

```


```{r}

# Extract nonzero coefficients from a glmnet/cv.glmnet object
make_coef_table <- function(fit, s = "lambda.1se", drop_intercept = TRUE) {
  coefs <- coef(fit$final_cv, s = s)
  coefs_mat <- as.matrix(coefs)

  out <- data.frame(
    variable = rownames(coefs_mat),
    coefficient = as.numeric(coefs_mat[, 1]),
    row.names = NULL
  )

  out <- out[out$coefficient != 0, ]

  if (drop_intercept) {
    out <- out[out$variable != "(Intercept)", ]
  }

  out
}


# Table 1: with AGE, BMI, SEX
final_coef_table <- make_coef_table(res_mash, s = "lambda.1se")

# Table 2: without AGE, BMI, SEX
final_coef_table_A <- make_coef_table(res_mash_A, s = "lambda.1se")

# Print tables
print("Signatures without covariates")
print(final_coef_table_A)

print("Signatures with covariates")
print(final_coef_table)
```


##  FIBROSIS
```{r}
# STEP 1: Fibrosis-specific selection
common_samples_fib <- intersect(rownames(histo), colnames(liver_stages))
liver_lasso_fib <- liver_stages[rownames(liver_stages), common_samples_fib]
liver_lasso_fib <- t(liver_lasso_fib)

# STEP 2: Fibrosis data prep
liver_data_fib <- data.frame(row.names = rownames(liver_lasso_fib))
liver_data_fib$Fibrosis <- histo[rownames(liver_lasso_fib), "Fibrosis"]
liver_data_fib$Fibrosis_F2 <- ifelse(liver_data_fib$Fibrosis >= 2, 1, 0)

met_names_fib <- colnames(liver_lasso_fib)
for (i in 1:length(met_names_fib)) {
  liver_data_fib[[met_names_fib[i]]] <- liver_lasso_fib[, i]
}

# Add covariates
liver_data_fib$age <- clinical_var[rownames(liver_lasso_fib), "Age"]
liver_data_fib$sex <- clinical_var[rownames(liver_lasso_fib), "SEX"]
liver_data_fib$BMI <- clinical_var[rownames(liver_lasso_fib), "BMI"]

# Clean data
liver_data_fib <- na.omit(liver_data_fib)

# Convert sex to binary here as well
liver_data_fib$sex <- as.character(liver_data_fib$sex)
liver_data_fib$sex_bin <- ifelse(
  liver_data_fib$sex %in% c("M", "Male", "MALE", "m", "male", 1, "1"),
  1,
  0
)

# DIAGNOSTIC: Check binary outcome & sample size
cat("Fibrosis_F2 distribution:\n")
print(table(liver_data_fib$Fibrosis_F2))
cat("Final N =", nrow(liver_data_fib), "\n")
cat("Metabolites =", length(met_names_fib), "\n")

# Run analysis
met_cols_fib <- setdiff(colnames(liver_data_fib),
                        c("Fibrosis_F2", "Fibrosis", "age", "sex", "sex_bin", "BMI"))

res_fib_A <- run_lasso_nested_cv_modelA(
  data = liver_data_fib,
  outcome_bin_col = "Fibrosis_F2",
  metabolite_cols = met_cols_fib,
  alpha = 1,
  n_outer = 10,
  n_inner = 10
)


res_fib <- run_lasso_nested_cv(
  data            = liver_data_fib,
  outcome_bin_col = "Fibrosis_F2",
  covariate_cols  = c("age", "sex_bin", "BMI"),
  metabolite_cols = met_cols_fib,
  alpha           = 1,
  n_outer         = 10,
  n_inner         = 10
)

final_coef <- coef(res_mash $final_cv, s = "lambda.1se")
final_coefA <- coef(res_mash_A $final_cv, s = "lambda.1se")
```

                                                                           
```{r}

plot(res_fib$final_cv)
abline(v = log(res_fib$final_cv$lambda.min), col = "blue", lty = 2)
abline(v = log(res_fib$final_cv$lambda.1se), col = "red", lty = 2)
legend("topright",
       legend = c("lambda.min", "lambda.1se"),
       col = c("blue", "red"),
       lty = 2,
       bty = "n")



plot(res_fib_A$final_cv)
abline(v = log(res_fib_A$final_cv$lambda.min), col = "blue", lty = 2)
abline(v = log(res_fib_A$final_cv$lambda.1se), col = "red", lty = 2)
legend("topright",
       legend = c("lambda.min", "lambda.1se"),
       col = c("blue", "red"),
       lty = 2,
       bty = "n")


# Table 1: with AGE, BMI, SEX
final_coef_table <- make_coef_table(res_fib, s = "lambda.1se")

# Table 2: without AGE, BMI, SEX
final_coef_table_A <- make_coef_table(res_fib_A, s = "lambda.1se")

# Print tables
print("Table 2: without AGE, BMI, SEX")
print(final_coef_table_A)

print("Table 2: with AGE, BMI, SEX")
print(final_coef_table)

mash_met <- res_mash$final_selected_metabolites
fib_met <- res_fib$final_selected_metabolites
```

## =============================================================================
## DIABLO (LIVER METABOLITES SIGNATURES OF MASH AND FIBROSIS WITH PLASMA METABOLITES) 
## =============================================================================

### Preprocessing 
```{r}
library(mixOmics)

# 1. Define common samples and selected liver metabolites


common_samples <- intersect(colnames(liver_stages), colnames(plasma_stages))
common_samples <- intersect(rownames(histo), common_samples)

mets_names <- unique(c(mash_met, fib_met))
mets_names <- gsub("`", "", mets_names)

cat("Total number of liver metabolites predicted by lasso associated with MASH and Fibrosis:", length(mets_names), "\n")
cat("Total number of common samples:", length(common_samples), "\n")

mets_names_used <- intersect(mets_names, rownames(liver_stages))

liver_sub  <- liver_stages[mets_names, common_samples, drop = FALSE]
plasma_sub <- plasma_stages[, common_samples, drop = FALSE]


# 2. Transpose to samples x variables


liver_sub  <- t(liver_sub)
plasma_sub <- t(plasma_sub)


# 3. Remove plasma columns that are all NA


all_na_plasma <- apply(plasma_sub, 2, function(x) all(is.na(x)))
plasma_sub <- plasma_sub[, !all_na_plasma, drop = FALSE]

# Optional: remove all-NA liver columns too, just in case
all_na_liver <- apply(liver_sub, 2, function(x) all(is.na(x)))
liver_sub <- liver_sub[, !all_na_liver, drop = FALSE]


# 4. Global scaling


liver_scaled  <- scale(liver_sub, center = TRUE, scale = TRUE)
plasma_scaled <- scale(plasma_sub, center = TRUE, scale = TRUE)


# 5. Remove zero-variance columns

liver_var  <- apply(liver_scaled, 2, var, na.rm = TRUE)
plasma_var <- apply(plasma_scaled, 2, var, na.rm = TRUE)

liver_final  <- liver_scaled[, liver_var > 0, drop = FALSE]
plasma_final <- plasma_scaled[, plasma_var > 0, drop = FALSE]


# 6. Make names block-specific and unique


colnames(liver_final)  <- paste0("liv_", colnames(liver_final))
colnames(plasma_final) <- paste0("pls_", colnames(plasma_final))

colnames(liver_final)  <- make.unique(colnames(liver_final))
colnames(plasma_final) <- make.unique(colnames(plasma_final))

# Check uniqueness
stopifnot(length(unique(colnames(liver_final)))  == ncol(liver_final))
stopifnot(length(unique(colnames(plasma_final))) == ncol(plasma_final))


# 7. Build outcome vector aligned to sample IDs


common_ids <- rownames(liver_final)

# For MASH
liver_outcome <- histo[common_ids, "Mash", drop = TRUE]

# If instead you want fibrosis:
# liver_outcome <- ifelse(histo[common_ids, "Fibrosis", drop = TRUE] >= 2, 1, 0)

liver_outcome <- factor(liver_outcome, levels = c(0, 1), labels = c("No", "Yes"))

# Sanity checks
stopifnot(nrow(liver_final) == length(liver_outcome))
stopifnot(nrow(plasma_final) == nrow(liver_final))
```

## DIABLO BLock 
```{r}

# 8. Build DIABLO input


X <- list(
  liver  = liver_final,
  plasma = plasma_final
)

design <- matrix(
  c(0, 0.5,
    0.5, 0),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(names(X), names(X))
)


# 9. Tune DIABLO


 test.keep <- list(
   
   liver = c(3, 5),
   plasma = c(5,10)
 )

set.seed(123)

tune <- tune.block.splsda(
  X = X,
  Y = liver_outcome,
  ncomp = 2,
  design = design,
  folds = 10,
  nrepeat = 10,
  dist = "centroids.dist",
  test.keepX = test.keep
)

print(tune)


# 10. Fit final DIABLO model

diablo <- block.splsda(
  X = X,
  Y = liver_outcome,
  ncomp = 2,
  design = design,
  keepX = tune$choice.keepX,
  near.zero.var = TRUE 
)

# DIABLO's own stability (no FDR needed)
perf_diablo <- perf(diablo, validation = "Mfold", folds = 10, nrepeat = 10)
plot(perf_diablo)
```

### Optional plots

```{r}
#pdf("sample_seperation_plot.pdf", width = 16, height = 10)
plotIndiv(
  diablo, 
  comp = 1:2, 
  group = liver_outcome, 
  ind.names = FALSE,
  legend = TRUE
)
#dev.off ()


plotVar(diablo, comp = 1:2, blocks = c("liver", "plasma"))

#pdf("circos_diablo.pdf", width = 20, height = 16)
circosPlot(
  diablo,
  comp = 1:2,
  cutoff = 0.5,
  size.variables = 0.78,
  size.legend = 0.7,
  ncol.legend = 1
)
#dev.off()

```



## =============================================================================
## DIABLO (LIVER, PLASMA METABOLITES WITH LIVER TRANSCRIPTOMICS) 
## =============================================================================

```{r}
transcripto <- plasma_transcript_DA
metabolite <- liver_matrix
plasma_metabolite <- plasma_matrix

common <- intersect(colnames(transcripto), colnames(metabolite))
common <- intersect(colnames(plasma_metabolite), common)
length(common)
 
metabolite <- map_liver_bio(metabolite[, common, drop= FALSE], sample_info_liver)
plasma_metabolite <- map_chemical_names(plasma_metabolite[,common, drop = FALSE], chemical_details)

transcripto <- transcripto[,common, drop=FALSE]

transcripto <- t(transcripto)
metabolite<- t(metabolite)
plasma_metabolite <- t(plasma_metabolite)

plasma_metabolite_scaled <- scale(plasma_metabolite, center = TRUE, scale = TRUE)
transcripto_scaled  <- scale(transcripto, center = TRUE, scale = TRUE)
metabolite_scaled <- scale(metabolite, center = TRUE, scale = TRUE)


# For MASH
cluster_id <- cluster_data[common, "class" , drop = TRUE]


```

```{r}
# -------------------------------------------------------------------
# 8. Build DIABLO input
# -------------------------------------------------------------------

X <- list(
  Liver = metabolite_scaled,
  Transcripto = transcripto_scaled,
  Plasma = plasma_metabolite_scaled
)

design <- matrix(
  c(0,   0.5, 0.5,
    0.5, 0,   0.5,
    0.5, 0.5, 0),
  nrow = 3,
  byrow = TRUE,
  dimnames = list(names(X), names(X))
)
# -------------------------------------------------------------------
# 9. Tune DIABLO
# -------------------------------------------------------------------

 test.keep <- list(
   
   Liver = c(5, 10),
   Transcripto  = c(5,10),
   Plasma = c(5,10)
 )

set.seed(123)

tune <- tune.block.splsda(
  X = X,
  Y = cluster_id,
  ncomp = 2,
  design = design,
  folds = 10,
  nrepeat = 10,
  dist = "centroids.dist",
  test.keepX = test.keep
)

print(tune)

# -------------------------------------------------------------------
# 10. Fit final DIABLO model
# -------------------------------------------------------------------

diablo <- block.splsda(
  X = X,
  Y = cluster_id,
  ncomp = 2,
  design = design,
  keepX = tune$choice.keepX,
  near.zero.var = TRUE 
)

perf_diablo <- perf(diablo, validation = "Mfold", folds = 10, nrepeat = 10)
plot(perf_diablo)

circosPlot(
  diablo,
  comp = 1:2,
  cutoff = 0.5,
  size.variables = 0.7)

plotIndiv(
  diablo, 
  comp = 1:2, 
  group = cluster_id, 
  ind.names = FALSE,
  legend = TRUE
)
```



```{r}

#sessionInfo()
#sink("sessionInfo.txt")
```




