
#-------------------------------------------------------------------------------
# Function to perform khi-2 or fisher exact test and return p-value
#-------------------------------------------------------------------------------

# Following is the description of the parameters used : 
# data : Table containing variables
# varClustName : Name of the variable containing the clusters
# grpes : List containing the labels of the 2 clusters we want to compare
# var : Name of the variable compared between the 2 groups specified in grpes

khiFishFct = function(data, varClustName, grpes, var){
  
  freqTableToCompare <- table(data[data[,varClustName] %in% grpes, varClustName], data[data[,varClustName] %in% grpes,var])
  test = chisq.test(freqTableToCompare)
  
  # If the hypothesis is not verified, we use fisher exact test instead of khi-2 test
  if (length(which(test$expected < 5)) > 0){
    
    # If more than 2x2 dimension, use Monte-Carlo method to simulate the p-value
    if (any(dim(freqTableToCompare) > 2)){
      pval <- fisher.test(freqTableToCompare, simulate.p.value = TRUE)$p.value
      testUsed <- "Fisher with MonteCarlo"
    } else {
      pval <- fisher.test(freqTableToCompare, simulate.p.value = FALSE)$p.value
      testUsed <- "Fisher Exact Test"
    }
  } else {
    pval <- test$p.value
    testUsed <- "Khi-2 Test"
  }
  
  return(list(pval = pval, testUsed = testUsed))
}


#-------------------------------------------------------------------------------
# Function to perform Kruskal Wallis test and return p-value
#-------------------------------------------------------------------------------

# Kruslal tests
kwFct = function(data, varTest, varCible){
  kwPval = kruskal.test(formula = as.formula(paste(varTest, varCible, sep = " ~ ")),
                        data = data[!is.na(data[,varTest]),])
  akwPval = kwPval$p.value
  return(akwPval)
}
#-------------------------------------------------------------------------------