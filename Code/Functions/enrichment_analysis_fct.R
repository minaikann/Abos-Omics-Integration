map_to_hmdb <- function(metabolite_names, info_table, name_col = "BIOCHEMICAL") {
  # 1. Clean the reference HMDB column
  info_table$HMDB_clean <- sapply(strsplit(as.character(info_table$HMDB), ","), `[`, 1)
  
  # 2. Filter out rows with no ID
  valid_table <- info_table[!is.na(info_table$HMDB_clean) & info_table$HMDB_clean != "", ]
  
  # 3. Create lookup vector (Name -> ID)
  # Uses the name_col argument to select the correct column
  lookup <- setNames(valid_table$HMDB_clean, valid_table[[name_col]])
  
  # 4. Map the input names to IDs
  mapped_ids <- lookup[metabolite_names]
  
  # 5. Remove NAs, add 'hmdb:' prefix, and return the vector
  final_vector <- paste0("hmdb:", na.omit(mapped_ids))
  
  return(final_vector)
}

run_pathway_enrichment <- function(hmdb_ids, pathway_sources = c("kegg", "reactome" )) {
  library(RaMP)
  
  # 1. Run Fisher's Test
  d <- runEnrichPathways(analytes = hmdb_ids)
  
  # 2. Filter by p-value
  d2 <- filterEnrichResults(d)
  
  # 3. Filter by source (dynamic based on input argument)
  d2$fishresults <- d2$fishresults[d2$fishresults$pathwaySource %in% pathway_sources, ]
  
  
  return(d2$fishresults)
}

# Filter & rank top 10 pathways (example: plasma subset)
plot_pathway_enrichment <- function(results_df, plot_title = "Enriched Pathways") {
  library(ggplot2)
  library(dplyr)
  
  # 1. Filter, sort, and slice for top 15
  plot_data <- results_df %>%
    filter(Pval_FDR <= 0.05) %>%
    arrange(Pval_FDR) %>%
    slice(1:15) %>%
    mutate(log10FDR = -log10(Pval_FDR))
  
  # 2. Create the plot
  p <- ggplot(plot_data, aes(x = log10FDR, y = reorder(pathwayName, log10FDR), fill = pathwaySource)) +
    geom_col(width = 0.7) +
    scale_x_continuous(name = expression(-log[10]~"(FDR-adjusted p-value)")) +
    scale_fill_brewer(name = "Pathway source", palette = "Set1") +
    theme_minimal() +
    theme(
      axis.title.y = element_blank(),
      axis.text.y  = element_text(size = 10),
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(title = plot_title)
  
  return(p)
}

