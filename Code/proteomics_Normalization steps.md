---
title: "Proteomics_data Normalization"
author: "Deborah Mina Ikann"
output: html_document
---
```{r}
library(openxlsx)
library(dplyr)

raw <- read.xlsx("~/Documents/Clustering_ABOS/codes/Data/proteiniomics_baseline.xlsx", rowNames = FALSE)

```
<img width="261" height="240" alt="baseline_data_preprocessing" src="https://github.com/user-attachments/assets/fd429c62-d80d-4724-bd9f-6283b7d27537" />


```{r}

meta_path <- "~/Documents/Clustering_ABOS/..."
sheet_names <- getSheetNames(meta_path)

all_sheets <- lapply(2:17, function(s) {
  dat <- read.xlsx(meta_path, sheet = s)     # headers ABOS_ID / OK read directly
  names(dat) <- trimws(names(dat))
  dat$Boite <- sheet_names[s]                # keep the box name (e.g. "Boîte 1")
  dat
})

combined <- bind_rows(all_sheets)

# Drop empty rows
combined <- combined[!is.na(combined$ABOS_ID), , drop = FALSE]

head(combined)
dim(combined)
```

```{r}
# If it's Excel (.xlsx):
library(readxl)

raw <- as.data.frame(raw)

# ALWAYS look at the corner first:
raw[1:4, 1:6]
```
```{r}
n_id_cols <- 3

# Reference table
ref <- raw[, 1:n_id_cols]
colnames(ref) <- c("Uniprot", "Protein_Name", "Gene")
ref[] <- lapply(ref, trimws)

# Numeric matrix — samples start at column 4
mat <- as.matrix(raw[, 4:ncol(raw)])
mode(mat) <- "numeric"
colnames(mat) <- colnames(raw)[4:ncol(raw)]   # "LFQ intensity 1", ... unchanged

# Index by Uniprot FOR NOW (needed for the drop/sanity checks below)
rownames(mat) <- ref$Uniprot

# ---- Remove the stray "Proteomics ID" row ----
# That row's Uniprot cell is "Proteomics ID" (or NA/blank), not a real accession
drop_row <- is.na(ref$Uniprot) | ref$Uniprot == "" | ref$Uniprot == "NA." | ref$Uniprot == "Uniprot" 
mat <- mat[!drop_row, , drop = FALSE]
ref <- ref[!drop_row, , drop = FALSE]

keep_cols <- grepl("^[0-9]", colnames(mat))
mat <- mat[, keep_cols, drop = FALSE]

# Sanity checks (still on Uniprot rownames)
stopifnot(nrow(ref) == nrow(mat))
stopifnot(all(ref$Uniprot == rownames(mat)))
cat("Matrix:", nrow(mat), "proteins x", ncol(mat), "columns\n")
colnames(mat)[1:5]    # LFQ intensity 1 ... 5
head(ref)

# ------------------------------------------------------------
# Switch matrix rownames from Uniprot -> Protein_Name
# Fall back to Uniprot where Protein_Name is missing/blank,
# and make.unique() to avoid duplicate row.names errors.
# ------------------------------------------------------------
# new_names <- ref$Protein_Name
# new_names <- ifelse(is.na(new_names) | new_names == "", ref$Uniprot, new_names)
# rownames(mat) <- make.unique(new_names)
# 
# # Record the exact rowname used, so ref still maps name <-> Uniprot <-> Gene
# ref$RowName <- rownames(mat)
# 
# # ref and mat stay aligned by position from here on
# stopifnot(nrow(ref) == nrow(mat))
# stopifnot(all(ref$RowName == rownames(mat)))

#write.xlsx(ref, "protein_chemical_details.xlsx")

# Safe mapper: Uniprot -> Protein_Name with dedup + NA fallback
map_protein <- function(data_matrix, ref) {
  match_idx <- match(rownames(data_matrix), ref$Uniprot)
  bio_names <- ref$Protein_Name[match_idx]
  new_names <- ifelse(!is.na(match_idx) & !is.na(bio_names) & bio_names != "",
                      bio_names, rownames(data_matrix))
  rownames(data_matrix) <- make.unique(new_names)
  return(data_matrix)
}
```

```{r}
# Look up a protein by Uniprot
ref[ref$Uniprot == "P62258", ]

# Find a gene
ref[ref$Gene == "YWHAE", ]

# Save it for your records
#write.csv(ref, "protein_reference_table.csv", row.names = FALSE)
```





```{r}
library(matrixStats)

cat("Total proteins before processing", nrow(mat))
# ------------------------------------------------------------
# Step 1: Remove proteins where > 30% of samples are MISSING
# Missing = NA. A 0 is a real measured value and is NOT missing,
# so it is counted as observed here (kept).
# ------------------------------------------------------------
frac_na <- rowMeans(is.na(mat))
keep_rows <- frac_na <= 0.30
mat <- mat[keep_rows, , drop = FALSE]
ref <- ref[keep_rows, , drop = FALSE]
cat("Proteins kept after >30% NA filter:", nrow(mat), "\n")

# ------------------------------------------------------------
# Step 2: Log2 transformation
# +1 pseudocount keeps real zeros finite (log2(0+1)=0).
# NA values stay NA through the transform.
# ------------------------------------------------------------
mat <- log2(mat+1)

# ------------------------------------------------------------
# Step 3: Remove outlier samples (ID 714, 715 & first 24)
# ------------------------------------------------------------
print(colnames(mat))
outlier_ids <- c("714", "715")
drop_cols   <- which(colnames(mat) %in% outlier_ids)
drop_cols   <- union(drop_cols, 1:24)
cat("Dropping these columns:\n")
print(colnames(mat)[drop_cols])
if (length(drop_cols) > 0) {
  mat <- mat[, -drop_cols, drop = FALSE]
}
cat("Samples kept after outlier removal:", ncol(mat), "\n")

# ------------------------------------------------------------
# Map each remaining sample to its box (needed for impute + batch)
# ------------------------------------------------------------
batch <- combined$Boite[match(colnames(mat), combined$ABOS_ID)]
stopifnot(!any(is.na(batch)))   # this is the box-lookup check, not data NAs

# ------------------------------------------------------------
# Step 4: Downshift imputation BY BOX
# Only NA values are imputed. Real zeros are observed values
# and contribute to mu/sigma like any other measurement.
# ------------------------------------------------------------
downshift_impute_box <- function(mat, batch, shift = 1.8, width = 0.3) {
  for (b in unique(batch)) {
    cols <- which(batch == b)
    block <- mat[, cols, drop = FALSE]
    obs   <- block[!is.na(block)]           # observed = non-NA (zeros included)
    mu    <- mean(obs)
    sigma <- sd(obs)
    is_miss <- is.na(block)
    n_miss  <- sum(is_miss)
    if (n_miss > 0) {
      block[is_miss] <- rnorm(n_miss, mean = mu - shift * sigma, sd = width * sigma)
    }
    mat[, cols] <- block
  }
  mat
}

set.seed(42)
mat <- downshift_impute_box(mat, batch)
cat("Any NAs left after imputation?", any(is.na(mat)), "\n")

# ------------------------------------------------------------
# Step 5: Batch correction by box (median centering)
# ------------------------------------------------------------
median_center_batch <- function(mat, batch) {
  grand_median <- median(mat, na.rm = TRUE)
  for (b in unique(batch)) {
    cols <- which(batch == b)
    batch_median <- median(mat[, cols], na.rm = TRUE)
    mat[, cols] <- mat[, cols] - batch_median + grand_median
  }
  mat
}
mat <- median_center_batch(mat, batch)

cat("Done. Final matrix:", nrow(mat), "proteins x", ncol(mat), "samples\n")

# Rownames are now Protein_Name (deduplicated) -> keep them on write
#write.xlsx(mat, "proteomics_matrix.xlsx", rowNames = TRUE)
```


```{r}
# Overall distribution of all intensity values
hist(mat, breaks = 100, main = "Distribution of log2 intensities",
     xlab = "log2 intensity")

# Per-sample distribution — boxplot, one box per sample column
# (this is the key plot to see batch/box effects before correction)
boxplot(mat[,1:5], las = 2, outline = FALSE,
        main = "Per-sample intensity distribution",
        ylab = "log2 intensity", cex.axis = 0.5)
```
<img width="700" height="432" alt="image" src="https://github.com/user-attachments/assets/c650a0c4-d0c7-4c54-9161-0bd4a8e9b4f3" />

```{r}
batch <- combined$Boite[match(colnames(mat), combined$ABOS_ID)]

plot(density(mat[,1], na.rm = TRUE), ylim = c(0, 0.15),
     main = "Per-sample density by box", xlab = "log2 intensity",
     col = as.integer(factor(batch))[1])
for (i in 2:ncol(mat)) {
  lines(density(mat[,i], na.rm = TRUE),
        col = as.integer(factor(batch))[i])
}
```
