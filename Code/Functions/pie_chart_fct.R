
# The pieFct and pieFct3cl functions are used to produce the pie chart, which is just a bar plot transformed into polar coordinates.

# Following is the description of the parameters used : 

# dataFreq : Frequency table
# col : indicates the color palette
# persLegend : Legend parameter
# cohort : Cohort name
# labelSize : label Size 



# ------------------------------------------------------------------------------
# Pie chart 6 clusters : 
# ------------------------------------------------------------------------------


pieFct = function(dataFreq, col, persLegend = FALSE, cohort = "Other", labelSize = 6){
  
  # Input table : Frequency table with first column indicating the cluster and the second one indicating its frequency
  names(dataFreq) = c("cluster", "Freq")
  
  
  # Calculate the positions for the labels
  if (cohort == "MAFALDA cohort"){    # Special case
    gradEt = c(dataFreq$Freq[5]/2,
               dataFreq$Freq[5] + dataFreq$Freq[4]/2,
               dataFreq$Freq[5] + dataFreq$Freq[4] + dataFreq$Freq[3]/2,
               dataFreq$Freq[5] + dataFreq$Freq[4] + dataFreq$Freq[3] + dataFreq$Freq[2]/2,
               dataFreq$Freq[5] + dataFreq$Freq[4] + dataFreq$Freq[3] + dataFreq$Freq[2] + dataFreq$Freq[1]/2)
  } else {
    gradEt = c(dataFreq$Freq[6]/2,
               dataFreq$Freq[6] + dataFreq$Freq[5]/2,
               dataFreq$Freq[6] + dataFreq$Freq[5] + dataFreq$Freq[4]/2,
               dataFreq$Freq[6] + dataFreq$Freq[5] + dataFreq$Freq[4] + dataFreq$Freq[3]/2,
               dataFreq$Freq[6] + dataFreq$Freq[5] + dataFreq$Freq[4] + dataFreq$Freq[3] + dataFreq$Freq[2]/2,
               dataFreq$Freq[6] + dataFreq$Freq[5] + dataFreq$Freq[4] + dataFreq$Freq[3] + dataFreq$Freq[2] + dataFreq$Freq[1]/2)
  }
  
  options(repr.plot.width = 2, repr.plot.height =2)
  
  # Add the % to the legend
  if (persLegend){
    
    dataFreq$partInd = prop.table(dataFreq$Freq)
    dataFreq$partRound = round(dataFreq$partInd*100)
    dataFreq$partRound = paste(dataFreq$partRound, "%")
    
    # Calculate IC 95%
    dataFreq$ec = 1.96 * sqrt(dataFreq$partInd*(1-dataFreq$partInd)/sum(dataFreq$Freq))
    dataFreq$IC = paste0("(", round((dataFreq$partInd-dataFreq$ec)*100, 1), "%-", round((dataFreq$partInd+dataFreq$ec)*100, 1), "%)")
    
    # Legend parameter
    nameLeg = paste0(cohort, " (n=", sum(dataFreq$Freq), ")")
    dataFreq[,"nameLeg"] = paste0("cluster ", dataFreq$cluster, " n=", dataFreq$Freq)
    
  }
  
  ## pie chart
 
  if (persLegend){
    
    bar <- ggplot(dataFreq, aes(x="", y=Freq, fill=nameLeg)) + 
      geom_bar(stat="identity") +
      
   
      coord_polar("y", start=0) +
      
      # Rename the legend
      labs(fill = nameLeg) +
      
      # Add the labels with the number of individuals, using the scale defined above
      geom_text(aes(y = gradEt,
                    label = rev(partRound)),
                color = "white", size = labelSize)
    
  } else {
    
    bar <- ggplot(dataFreq, aes(x="", y=Freq, fill=cluster)) + 
      geom_bar(stat="identity") +
      
 
      coord_polar("y", start=0)
    
  }
  
  bar = bar +
    scale_fill_manual(values = col[as.numeric(as.character(unique(dataFreq$cluster)))]) +
    
    # Theme 
    theme_minimal()+
    theme(
      axis.title.x = element_blank(),
      axis.text.x=element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_blank(),
      panel.grid=element_blank(),
      axis.ticks = element_blank(),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = labelSize*2+3),
      legend.text=element_text(size=labelSize*2)
    )
  
  return(bar)
}
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Pie chart 3 clusters 
# ------------------------------------------------------------------------------

pieFct3cl = function(dataFreq, col, cohort = "Other", labelSize = 6){
  
  # Input table : Frequency table
  names(dataFreq) = c("cluster", "Freq")
  dataFreq <- dataFreq[order(dataFreq$cluster),]
  nameLeg = paste0(cohort, " (n=", sum(dataFreq$Freq), ")")
  dataFreq[,"nameLeg"] = paste0("cluster ", dataFreq$cluster, " n=", dataFreq$Freq)
  dataFreq$nameLeg <- factor(dataFreq$nameLeg, levels = dataFreq$nameLeg) 
  
  dataFreq$partInd = prop.table(dataFreq$Freq)
  dataFreq$partRound = round(dataFreq$partInd*100)
  dataFreq$partRound = paste(dataFreq$partRound, "%")
  
  # Calculate the positions for the labels
  gradEt = c(dataFreq$Freq[3]/2,
             dataFreq$Freq[3] + dataFreq$Freq[2]/2,
             dataFreq$Freq[3] + dataFreq$Freq[2] + dataFreq$Freq[1]/2)
  
  options(repr.plot.width = 2, repr.plot.height =2)
  
  ## Pie chart
  bar <- ggplot(dataFreq, aes(x="", y=Freq, fill=nameLeg)) + 
    geom_bar(stat="identity") +
    
    coord_polar("y", start=0) +
    
    # Rename the legend
    labs(fill = nameLeg) +
    
    # Add the labels with the number of individuals, using the scale defined above
    geom_text_repel(aes(y = gradEt,
                        label = rev(partRound)), size = labelSize, color = "black")
  
  bar = bar +
    scale_fill_manual(values = col) +
    
    # Theme 
    theme_minimal()+
    theme(
      axis.title.x = element_blank(),
      axis.text.x=element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_blank(),
      panel.grid=element_blank(),
      axis.ticks = element_blank(),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = labelSize*2+3),
      legend.text=element_text(size=labelSize*2)
    )
  
  return(bar)
}
# ------------------------------------------------------------------------------

