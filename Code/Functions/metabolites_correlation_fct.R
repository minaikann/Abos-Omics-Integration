# Mapping functions (unchanged)
map_plasma_bio <- function(data_matrix, chemical_details) {
  match_idx <- match(rownames(data_matrix), chemical_details$COMP_ID)
  bio_names <- chemical_details$CHEMICAL_NAME[match_idx]
  new_matrix <- data_matrix
  rownames(new_matrix) <- ifelse(!is.na(match_idx), bio_names, rownames(new_matrix))
  return(new_matrix)
}

map_liver_bio <- function(data_matrix, sample_info_liver) {
  match_idx <- match(rownames(data_matrix), sample_info_liver$COMP.ID)
  bio_names <- sample_info_liver$BIOCHEMICAL[match_idx]
  new_matrix <- data_matrix
  rownames(new_matrix) <- ifelse(!is.na(match_idx), bio_names, rownames(new_matrix))
  return(new_matrix)
}

# FILTERING FUNCTION (unchanged)
filter_significant_mets <- function(results_list, logfc_threshold = 0.26, pval_col = "adj.P.Val") {
  sig_mets <- character(0)
  
  for(comp in names(results_list)) {
    df <- results_list[[comp]]
    logfc <- suppressWarnings(as.numeric(df[, "logFC"]))
    adj_p <- suppressWarnings(as.numeric(df[, pval_col]))
    
    significant <- !is.na(logfc) & !is.na(adj_p) & 
      abs(logfc) > logfc_threshold & adj_p < 0.05
    
    sig_mets <- unique(c(sig_mets, rownames(df)[significant]))
  }
  
  return(sig_mets)
}

plasma_liver_correlation <- function(
    plasma_bio, liver_bio, common_samples,
    top_plasma, top_liver, 
    fdr_threshold = 0.05
) {
  
  cat("Testing:", length(top_plasma), "×", length(top_liver), "correlations \n")  
  
  if(length(top_plasma) == 0 || length(top_liver) == 0) {
    return(list(message = "No significant metabolites", cor_df = data.frame(), sig_correlations = data.frame()))
  }
  
  # 0. ALIGN & SCALE both matrices to common samples (like cor_metab_histology)
  plasma_aligned <- plasma_bio[top_plasma, common_samples, drop = FALSE]
  liver_aligned <- liver_bio[top_liver, common_samples, drop = FALSE]
  
  cat("Samples used:", length(common_samples), "\n")
  cat("Plasma mets:", nrow(plasma_aligned), "×", ncol(plasma_aligned), "\n")
  cat("Liver mets:", nrow(liver_aligned), "×", ncol(liver_aligned), "\n")
  
  # Scale each metabolite (row) - SAME as cor_metab_histology
  plasma_scaled <- t(scale(t(plasma_aligned)))
  liver_scaled <- t(scale(t(liver_aligned)))
  
  cat("Scaled plasma dims:", dim(plasma_scaled), "\n")
  cat("Scaled liver dims:", dim(liver_scaled), "\n")
  
  # 1. Full correlation matrices
  cor_mat <- matrix(NA, nrow = length(top_plasma), ncol = length(top_liver))
  pval_mat <- matrix(NA, nrow = length(top_plasma), ncol = length(top_liver))
  rownames(cor_mat) <- top_plasma
  colnames(cor_mat) <- top_liver
  rownames(pval_mat) <- top_plasma
  colnames(pval_mat) <- top_liver
  
  counter <- 0
  for(i in seq_along(top_plasma)) {
    for(j in seq_along(top_liver)) {
      counter <- counter + 1
      
      p_met <- plasma_scaled[i, ]
      l_met <- liver_scaled[j, ]
      
      valid <- !is.na(p_met) & !is.na(l_met) & is.finite(p_met) & is.finite(l_met)
      
      if(sum(valid) >= 3) {
        test <- suppressWarnings(cor.test(p_met[valid], l_met[valid], method = "spearman", exact = FALSE))
        cor_mat[i, j] <- test$estimate
        pval_mat[i, j] <- test$p.value
      }
    }
    if(counter %% 50 == 0) cat("Processed", counter, "pairs\n")
  }
  
  # 2. FDR adjustment
  pvals_vec <- as.vector(pval_mat)
  fdr_vec <- p.adjust(pvals_vec, "BH")
  fdr_mat <- matrix(fdr_vec, nrow=nrow(pval_mat), ncol=ncol(pval_mat))
  rownames(fdr_mat) <- rownames(pval_mat)
  colnames(fdr_mat) <- colnames(pval_mat)
  
  # 3. Correlation data frame
  cor_df <- data.frame(
    Plasma_Metabolite = rep(top_plasma, each = length(top_liver)),
    Liver_Metabolite = rep(top_liver, times = length(top_plasma)),
    Spearman_R = as.vector(cor_mat),
    P_value = as.vector(pval_mat),
    FDR = as.vector(fdr_mat),
    N_samples = length(common_samples),
    stringsAsFactors = FALSE
  )
  
  sig_correlations <- cor_df[cor_df$FDR < fdr_threshold & !is.na(cor_df$FDR), ]
  
  cat("Significant correlations (FDR <", fdr_threshold, "):", nrow(sig_correlations), "\n")
  
  # 4. Return
  list(
    all_correlations = cor_df,
    sig_correlations = sig_correlations,
    top10 = head(sig_correlations[order(abs(sig_correlations$Spearman_R), decreasing = TRUE), ], 10),
    cor_matrix = cor_mat,
    pval_matrix = pval_mat,
    fdr_matrix = fdr_mat
  )
}


# === FUNCTION: Histology metabolite correlation 
cor_metab_histology <- function(metab_matrix, histology_t_matrix, name) {
  cat("=== ", name, "===\n")
  cat("Input metab dims:", dim(metab_matrix), "\n")
  cat("Input histo dims:", dim(histology_t_matrix), "\n")
  
  common_samples <- intersect(colnames(metab_matrix), colnames(histology_t_matrix))
  cat("Common samples:", length(common_samples), "\n")
  cat("Sample names match:", head(common_samples), "\n")
  
  metab_aligned <- metab_matrix[, common_samples, drop = FALSE]
  histo_aligned <- histology_t_matrix[, common_samples, drop = FALSE]
  
  metab_scaled <- t(scale(t(metab_aligned)))
  cat("Scaled metab dims:", dim(metab_scaled), "\n")
  
  n_stages <- nrow(histo_aligned)
  n_mets   <- nrow(metab_scaled)
  
  results <- data.frame(stringsAsFactors = FALSE)
  
  for (stage_i in 1:n_stages) {
    stage_name <- rownames(histo_aligned)[stage_i]
    cat("Processing stage:", stage_name, "\n")
    
    rho_vec    <- numeric(n_mets)
    pvals_vec  <- numeric(n_mets)
    
    suppressWarnings({
      for (i in 1:n_mets) {
        test <- cor.test(metab_scaled[i, ], histo_aligned[stage_i, ], method = "spearman")
        rho_vec[i]   <- test$estimate
        pvals_vec[i] <- test$p.value
      }
    })
    
    pvals_fdr <- p.adjust(pvals_vec, method = "BH")
    
    stage_results <- data.frame(
      Metabolite       = rownames(metab_matrix),
      Histology_Stages = stage_name,
      Samples_n        = length(common_samples),
      rho              = rho_vec,
      p_value          = pvals_vec,
      p_fdr            = pvals_fdr,
      Significant      = pvals_fdr < 0.05,
      stringsAsFactors = FALSE
    )
    
    results <- rbind(results, stage_results)
  }
  
  n_sig <- sum(results$Significant)
  cat(name, "sig metabolite-stage pairs (|ρ|>0.3 FDR<0.05):", n_sig, "\n")
  write.csv(results, paste0(name, "_vs_histology.csv"), row.names = FALSE)
  
  top_sig <- results[results$Significant | results$p_fdr < 0.1, ]
  cat("\n")
  
  return(results)
}

# === FUNCTION: Single heatmap===

priority_mets <- c(
  c(
    "1,5-anhydroglucitol (1,5-AG)",
    "3-(3-amino-3-carboxypropyl)uridine",
    "4-cholesten-3-one",
    "gamma-carboxyglutamate",
    "glycerophosphoinositol",
    "glycochenodeoxycholate",
    "glycosyl ceramide (d18:1/20:0, d16:1/22:0)",
    "glycosyl-N-stearoyl-sphingosine (d18:1/18:0)",
    "lactosyl-N-arachidoyl-sphingosine (d18:1/20:0)",
    "lactosyl-N-behenoyl-sphingosine (d18:1/22:0)",
    "lactosyl-N-nervonoyl-sphingosine (d18:1/24:1)",
    "lactosyl-N-stearoyl-sphingosine (d18:1/18:0)",
    "N-acetylglutamate",
    "N6-carbamoylthreonyladenosine",
    "noradrenaline",
    "stearoyl-arachidonoyl-glycerol (18:0/20:4) [2]",
    "stearoyl-arachidonoyl-glycerol (18:0/20:4) [2]*",
    "turine"
  )
  
)


make_stage_heatmap <- function(
    sig_results, name, histo_cols,
    priority_mets = NULL,
    n_total = 30,
    fdr_cutoff = 0.05
) {
  library(dplyr)
  library(reshape2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  
  # 1. Prepare data (keep p_fdr column, ensure numeric, no duplicates)
  data <- sig_results %>%
    filter(Histology_Stages %in% histo_cols) %>%
    mutate(
      rho   = suppressWarnings(as.numeric(rho)),
      p_fdr = suppressWarnings(as.numeric(p_fdr))
    ) %>%
    distinct(Metabolite, Histology_Stages, .keep_all = TRUE)
  
  if (nrow(data) == 0) {
    stop("No rows found for selected histology columns.")
  }
  
  # 2. Choose top metabolites (with priority)
  if (!is.null(priority_mets)) {
    priority_present <- intersect(priority_mets, unique(data$Metabolite))
    remaining_n <- max(0, n_total - length(priority_present))
    
    other_mets <- data %>%
      filter(!Metabolite %in% priority_present) %>%
      group_by(Metabolite) %>%
      summarise(max_abs_rho = max(abs(rho), na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(max_abs_rho)) %>%
      slice_head(n = remaining_n) %>%
      pull(Metabolite)
    
    top_mets <- unique(c(priority_present, other_mets))
  } else {
    top_mets <- data %>%
      group_by(Metabolite) %>%
      summarise(max_abs_rho = max(abs(rho), na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(max_abs_rho)) %>%
      slice_head(n = n_total) %>%
      pull(Metabolite)
  }
  
  # 3. Build rho and p_fdr matrices
  rho_matrix <- acast(
    data,
    Metabolite ~ Histology_Stages,
    value.var = "rho",
    fill = NA_real_
  )
  
  p_fdr_matrix <- acast(
    data,
    Metabolite ~ Histology_Stages,
    value.var = "p_fdr",
    fill = NA_real_
  )
  
  # 4. Ensure all histo_cols exist
  for (mc in setdiff(histo_cols, colnames(rho_matrix))) {
    rho_matrix <- cbind(rho_matrix, NA_real_)
    colnames(rho_matrix)[ncol(rho_matrix)] <- mc
  }
  
  for (mc in setdiff(histo_cols, colnames(p_fdr_matrix))) {
    p_fdr_matrix <- cbind(p_fdr_matrix, NA_real_)
    colnames(p_fdr_matrix)[ncol(p_fdr_matrix)] <- mc
  }
  
  # 5. Subset to top metabolites and histo_cols
  heatmap_data <- rho_matrix[top_mets, histo_cols, drop = FALSE]
  p_fdr_matrix  <- p_fdr_matrix[top_mets, histo_cols, drop = FALSE]
  
  # 6. Per‑cell masking: only blur if p_fdr > fdr_cutoff or NA
  mask <- is.na(p_fdr_matrix) | p_fdr_matrix > fdr_cutoff
  heatmap_data[mask] <- NA_real_
  
  # 7. Colormap
  col_fun <- colorRamp2(
    c(-1, 0, 1),
    c("darkblue", "white", "darkred")
  )
  
  # 8. Cell text: only where heatmap_data is not NA
  cell_fun <- function(j, i, x, y, width, height, fill) {
    val <- heatmap_data[i, j]
    if (!is.na(val)) {
      grid.text(
        sprintf("%.2f", val),
        x, y,
        gp = gpar(fontsize = 6, fontface = "plain")
      )
    }
  }
  
  # 9. Create heatmap
  ht <- Heatmap(
    heatmap_data,
    name = "Spearman ρ",
    col = col_fun,
    na_col = "grey95",           # grey for non‑significant cells
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    show_row_dend = TRUE,
    show_column_dend = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_side = "right",
    row_names_gp = gpar(fontsize = 10),
    row_names_max_width = unit(12, "cm"),
    column_names_side = "bottom",
    column_names_gp = gpar(fontsize = 12),
    column_names_rot = 90,
    column_names_max_height = unit(8, "cm"),
    width = ncol(heatmap_data) * unit(15, "mm"),
    height = nrow(heatmap_data) * unit(6, "mm"),
    cell_fun = cell_fun,
    heatmap_legend_param = list(
      title = "Spearman ρ",
      at = c(-1, -0.5, 0, 0.5, 1),
      legend_height = unit(5, "cm")
    ),
    column_title = name,
    column_title_gp = gpar(fontsize = 14, fontface = "bold")
  )
  
  return(ht)
}

find_exact_matches <- function(liver_results, plasma_results, histo_cols, 
                                       rho_threshold = 0.2, fdr_threshold = 0.05) {
  
  stage_order <- c("Steatosis", "Inflammation", "Ballooning", "Mash", "Fibrosis")
  
  # 1. FILTER significant
  liver_sig <- liver_results[liver_results$p_fdr < fdr_threshold & abs(liver_results$rho) > rho_threshold, ]
  plasma_sig <- plasma_results[plasma_results$p_fdr < fdr_threshold & abs(plasma_results$rho) > rho_threshold, ]
  
  exact_matches <- intersect(unique(liver_sig$Metabolite), unique(plasma_sig$Metabolite))
  
  if(length(exact_matches) == 0) return(NULL)
  
  # 2. Filter to exact matches only
  liver_complete <- liver_results[liver_results$Metabolite %in% exact_matches, ]
  plasma_complete <- plasma_results[plasma_results$Metabolite %in% exact_matches, ]
  
  final_matches <- intersect(unique(liver_complete$Metabolite), unique(plasma_complete$Metabolite))
  cat("Similar metabolites found in both tissues with significant correlaton :", length(final_matches), "\n")
  
  
  if(length(final_matches) == 0) return(NULL)
  
  
  #cat("Exact matches:\n")
  #print(final_matches)
  
  return(final_matches)
}


metabs<- c("1-myristoyl-2-palmitoyl-GPC (14:0/16:0)", "1-myristoylglycerol (14:0)", 
             "1-palmitoleoylglycerol (16:1)*", "1-palmitoyl-2-arachidonoyl-GPE (16:0/20:4)*",
             "1-palmitoyl-2-arachidonoyl-GPI (16:0/20:4)*", "1-palmitoyl-2-dihomo-linolenoyl-GPC (16:0/20:3n3 or 6)*",
             "1-palmitoyl-2-linoleoyl-GPE (16:0/18:2)", "1-palmitoyl-2-linoleoyl-GPI (16:0/18:2)",
             "1-palmitoyl-2-oleoyl-GPC (16:0/18:1)", "1-palmitoyl-2-oleoyl-GPE (16:0/18:1)",
             "1-palmitoyl-2-oleoyl-GPI (16:0/18:1)*", "1-palmitoyl-2-palmitoleoyl-GPC (16:0/16:1)*",
             "1-palmitoyl-GPG (16:0)*", "1-stearoyl-2-arachidonoyl-GPE (18:0/20:4)",
             "1-stearoyl-2-docosahexaenoyl-GPE (18:0/22:6)*", "1-stearoyl-2-linoleoyl-GPE (18:0/18:2)*",
             "1-stearoyl-2-linoleoyl-GPI (18:0/18:2)", "1-stearoyl-2-oleoyl-GPC (18:0/18:1)",
             "1-stearoyl-2-oleoyl-GPE (18:0/18:1)", "1-stearoyl-GPE (18:0)", "1-stearoyl-GPG (18:0)",
             "diacylglycerol (16:1/18:2 [2], 16:0/18:3 [1])*", "oleoyl-oleoyl-glycerol (18:1/18:1) [2]*",
             "palmitoyl-arachidonoyl-glycerol (16:0/20:4) [2]*", "palmitoyl-docosahexaenoyl-glycerol (16:0/22:6) [1]*",
             "palmitoyl-linoleoyl-glycerol (16:0/18:2) [2]*", "palmitoyl-oleoyl-glycerol (16:0/18:1) [2]*",
             "lactosyl-N-behenoyl-sphingosine (d18:1/22:0)*", "lactosyl-N-palmitoyl-sphingosine (d18:1/16:0)")




