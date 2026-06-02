
#-------------------------------------------------------------------------------
# Function to perform khi-2 or fisher exact tests + print barplots
#-------------------------------------------------------------------------------

#The barPlotTestsFct function represent the part of outcomes in each cluster, with confidence intervals and khi-2 or Fisher comparisons, for the patients for whom we have the information concerning outcomes.

# Following is the description of the parameters used : 
# data : Table containing clinical variables, clusters, and outcomes variables
# clust : Cluster variable name
# outcome : Outcome variable name
# position : Brace positions on the graph
# color : indicates the color palette
# ymax = 1 : ylim max 
# labY : y-axix label
# cohort = cohort name


barPlotTestsFct = function(data, clust, outcome, position, color, ymax = 1, labY = NULL, cohort = "Others", clustOrder = NULL){
  
  ## Calculate confidence intervals
  # Retrieve the number of patients with information regarding the cluster and the outcome
  cnt = cbind(table(data[!is.na(data[,clust]) & !is.na(data[,outcome]),clust]),
              data.frame(tapply(data[!is.na(data[,clust]) & !is.na(data[,outcome]), outcome],
                                data[!is.na(data[,clust]) & !is.na(data[,outcome]), clust],
                                function(r) sum(as.numeric(r)/length(r))))
  )
  names(cnt) = c("clust", "nb", "tx")
  
  # if the order of clusters is specified, use it
  if (!is.null(clustOrder))
    cnt <- cnt[order(match(cnt[,"clust"], clustOrder)),]
  
  # Calculate 
  cnt[,"ec"] = 1.96*sqrt(cnt[,"tx"]*(1-cnt[,"tx"])/cnt[,"nb"])
  # Lower bound
  cnt[,"borneInf"] = round(cnt[,"tx"]-cnt[,"ec"], 4)
  # Upper bound
  cnt[,"borneSup"] = round(cnt[,"tx"]+cnt[,"ec"], 4)
  
  # Special case
  if (cohort == "MAFALDA cohort"){
    cnt$clust = as.numeric(as.character(cnt$clust))
    cnt[nrow(cnt)+1,] = c(3, 0, 0, 0, 0, 0)
    cnt <- cnt[order(cnt$clust),]
  }
  # Factor to keep the right order
  cnt$clust = factor(cnt$clust, levels = cnt$clust)
  
  
  ## tests 2-2
  
  # Prepare the table with all combinations
  chiFishTests = data.frame(t(data.frame(combn(sort(unique(data[!is.na(data[,clust]),clust])), 2))))
  names(chiFishTests) = c("group1", "group2")
  row.names(chiFishTests) = c(1:nrow(chiFishTests))
  
  # Perform the tests and retrieve the p-values
  for (i in 1:nrow(chiFishTests)){
    chiFishTests[i, "p"] = khiFishFct(data, clust, c(chiFishTests[i,1], chiFishTests[i,2]), outcome)$pval
  }
  
  # Correction for multiple tests (bonferroni)
  chiFishTests$Adj.P.Value = p.adjust(chiFishTests$p, method = "bonferroni")
  
  # Number of stars based on the adjusted p-value
  chiFishTests$signif <- ifelse(chiFishTests$Adj.P.Value < 0.001, "***",
                                ifelse(chiFishTests$Adj.P.Value < 0.01, "$",
                                       ifelse(chiFishTests$Adj.P.Value < 0.025, "@",
                                              ifelse(chiFishTests$Adj.P.Value < 0.05, "*", "NS"))))
  
  
  # Adjust the braces on the graph
  chiFishTests = cbind(y.position = position, 
                       chiFishTests)
  
  # Graphe
  bar = ggplot(cnt, aes(x = clust, y = tx)) +
    geom_bar(stat = "identity", position = "dodge", fill = color) + 
    geom_errorbar(aes(ymin = tx-ec, ymax = tx+ec), width = 0.2) +  # Intervalles de confiance
    theme_bw() +  # Fond blanc
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_text(size = 15),
          axis.text.x =   element_text(size = 12),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 15)) +  # Format des étiquettes et titres des axes
    ylab(labY) +
    scale_y_continuous(labels = function(x) paste0(x*100, "%"), limits = c(0,ymax)) +
    
    # Tests (Stars corresponding to 2-2 tests with significant adjusted p-values)
    stat_pvalue_manual(chiFishTests[chiFishTests$signif != "NS",], label = "signif", label.size = 7)
  
  return(list(bar = bar, chiFishTests = chiFishTests))
}

#-------------------------------------------------------------------------------
