 
# ------------------------------------------------------------------------------
# Radarplot :
# ------------------------------------------------------------------------------

# The radarInFct function is used to produce the radarplot of clinical variables and represents the median values of age, BMI, HbA1c, LDL, triglycerides and ALT for each cluster. 

# Following is the description of the parameters used : 

# dataName : Table name containing the numeric variables we want to display in the radar plot, where each variable corresponds to a specific axis, along with the cluster variable
# varClustName : Cluster variable name
# dataAll : Table containing the numeric variables we want to use in the radar plot to calculate the percentiles 5 and 95, which we will name min and max
# color : indicates the color palette


radarInFct = function(dataName, varClustName, dataAll = get(dataName), color, clustOrder = NULL){
  
  # Calculate the medians of clinical variables for each cluster
  dataRadarN = sqldf(paste0('select ',varClustName ,', median(TGPUL) as ALT, median(hba1cpc) as HbA1c, median(triglymmolL) as Triglycerides, median(BMI) as BMI, median(AgeJourIntervention) as Age, median(cholest_LDLmmolL) as LDL from ', dataName,' group by ', varClustName))
  
  # if the order of clusters is specified, use it
  if (!is.null(clustOrder))
    dataRadarN <- dataRadarN[order(match(dataRadarN[,varClustName], clustOrder)),]
  
  # Calculate the percentiles 5 and 95 of clinical variables, which we will name min and max
  tempMinMax = apply(dataAll[,c("TGPUL", "hba1cpc", "triglymmolL", "BMI", "AgeJourIntervention", "cholest_LDLmmolL")], 2, 
                     function(r) quantile(r, probs = seq(0,1,0.05), na.rm = TRUE))
  dataRadarN[nrow(dataRadarN) + 1, varClustName] = "Min"
  dataRadarN[nrow(dataRadarN) + 1, varClustName] = "Max"
  dataRadarN[nrow(dataRadarN) - 1, 2:ncol(dataRadarN)] = tempMinMax[2,]
  dataRadarN[nrow(dataRadarN), 2:ncol(dataRadarN)] = tempMinMax[nrow(tempMinMax)-1,]
  
  graph = radarchart(dataRadarN[,-1],
                     
                     # Values displayed on the axes
                     axistype = 0,
                     palcex = 0.7,
                     
                     # Polygon options
                     pcol = c(color, scales::alpha("grey60", 0), "grey60"),
                     pfcol = c(scales::alpha(color, 0.4), scales::alpha(c("grey60", "grey60"), 0)),
                     plty = 1,
                     plwd = 2,
                     
                     # Customize the grid
                     cglty = 1,
                     cglwd = 0.8,
                     cglcol = "grey40",
                     
                     # Variable labels
                     vlcex = 1,
                     
                     # Axis options
                     maxmin = FALSE,
                     axislabcol = "grey20",
                     seg = 3,
                     
                     # Point options
                     pty = 20)
}

# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# Radarplots outcomes
# ------------------------------------------------------------------------------

# The radarOutFct function is used to calculate the proportion of patients by outcome for each cluster. 

# Following is the description of the parameters used : 

# data : Table containing outcomes variables
# varCibleBin : Outcome variable names
# varClustName : Cluster variable name


radarOutFct = function(data, varCibleBin, varClustName){
  
  dataRadar = data.frame(clust = sort(unique(data[,varClustName])))
  
  for (i in dataRadar$clust){
    for (v in varCibleBin){
      dataRadar[dataRadar$clust == i,v] = round(sum(as.numeric(data[!is.na(data[,v]) & data[,varClustName] == i, v]))/nrow(data[!is.na(data[,v]) & data[,varClustName] == i,]), 3)
    }
  }
  return(dataRadar)
}

# ------------------------------------------------------------------------------
