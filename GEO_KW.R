#this function servers to compare expression values of normal,tumor and metastatic gene chip samples using Kruskal-Wallis test
tnm_geo_KW <- function(tissueID, geneName){
  
  ################################data preparation#######################################
  
  #opening clinical table
  clinical_geo <- read_xlsx(path = "GEO_clinical.xlsx")
  clinical_geo <- as.data.frame(clinical_geo)
  NameExpTable <- paste0(tissueID,"-clean-ann.RData")
  
  #aim folder of expression table
  expWD <- paste0("TNMdata/GEO_Expression-Tables/", NameExpTable)
  
  #opening expression table
  load(expWD)
  expTable <- txtFile
  rm(txtFile)
  
  #screening for selected tissue in clinical table
  usedTissue <- clinical_geo[grepl(pattern = tissueID, clinical_geo$Origin),]
  sampNames <- usedTissue[,c(1,4,5,6)]
  sampNames$type <- NA
  
  #converting normal tumor meta columns to one column
  for (i in 1:nrow(sampNames)){
    if (sampNames$Normal[i] == 1){
      sampNames$type[i] <- "Normal"
    } else if(sampNames$Tumor[i] == 1){
      sampNames$type[i] <- "Tumor"
    } else if(sampNames$Meta[i] == 1){
      sampNames$type[i] <- "Metastatic"
    }
  }
  sampNames <-  sampNames[,-(2:4)]
 
  #creating expression table for the selected gene
  expvalue <- as.data.frame(t(expTable[expTable$gene_sym == geneName, colnames(expTable) %in% sampNames$GSM]))
  expvalue$GSM <- rownames(expvalue) 
  rownames(expvalue) <- c()
  
  #merging the clinical data and expression table - this data frame will be used for further investigations
  sampExp <- sampNames
  sampExp <- merge(x = sampExp, y = expvalue, by = "GSM")
  colnames(sampExp)[3] <- "Geneexp"
  sampExp$type <- factor(x = sampExp$type, levels =  c("Normal", "Tumor", "Metastatic"))
    
  ###################################statistical analysis##############################
  
  nNorm <- nrow(sampExp[sampExp$type == "Normal",])
  nTum <- nrow(sampExp[sampExp$type == "Tumor",])
  nMet <- nrow(sampExp[sampExp$type == "Metastatic",])
  nAll <- nNorm + nTum +nMet
  
  meanNorm <- mean(sampExp[sampExp$type == "Normal",3])
  sdnorm <- sd(sampExp[sampExp$type == "Normal",3])
  meanTum <- mean(sampExp[sampExp$type == "Tumor",3])
  sdtum <- sd(sampExp[sampExp$type == "Tumor",3])
  meanMet <- mean(sampExp[sampExp$type == "Metastatic",3])
  sdmet <- sd(sampExp[sampExp$type == "Metastatic",3])
  
  FCnt <- meanTum/meanNorm
  FCtm <- meanMet/meanTum
    
  kwt <- kruskal.test(x = list("Normal" = sampExp[sampExp$type == "Normal",3],
                               "Tumor" = sampExp[sampExp$type == "Tumor",3],
                               "Metastatic" = sampExp[sampExp$type == "Metastatic",3]))
  pvalue <- format(kwt$p.value, digits = 3, scientific = T)
    
  #identification of outliers
  normOuts <- boxplot.stats(sampExp[sampExp$type == "Normal",3])$stats
  tumOuts <- boxplot.stats(sampExp[sampExp$type == "Tumor",3])$stats
  metOuts <- boxplot.stats(sampExp[sampExp$type == "Metastatic",3])$stats
  outs <- range(normOuts, tumOuts, metOuts)
  
  
  ##############################data visualization with boxplots######################
  
  library(ggplot2)
  boxPlot <- ggplot(sampExp, aes(x = type, y = Geneexp, color = type, fill = type))  + 
    geom_boxplot(outlier.size = -1, notch = F, width = 0.55) + ylim(outs[1], outs[2]) +
    geom_jitter(width = 0.25) + 
    scale_color_manual( values =  c("darkgreen", "red", "darkorange3")) + 
    scale_fill_manual(values =   c("darkolivegreen1", "coral1", "goldenrod1")) +
    labs(x = NULL, y = paste0(geneName, " gene expression")) + 
    geom_label(label = paste0("P = ", pvalue), fill = "white", color = "black",x =3.2, y = outs[2] , label.size = 0, size = 7) +
    theme(legend.position = "none", 
          panel.background = element_rect(fill = NA, colour = "grey50"),  
          panel.grid = element_blank(),
          axis.text.x = element_text(size = 26, colour = "black"),
          axis.text.y = element_text(size = 26, colour = "black"),
          axis.title.y = element_text(size = 26, face = "bold",  vjust = 2))
  
  
  vioPlot <- ggplot(sampExp, aes(x = type, y = Geneexp, color = type, fill = type))  + 
    geom_violin(lwd = 1) + ylim(outs[1], outs[2]) +
    geom_boxplot(lwd = 1.2, outlier.size = -1, width = 0.1, color = c("darkgreen", "red", "darkorange3")) +
    scale_color_manual( values =  c("darkgreen", "red", "darkorange3")) + 
    scale_fill_manual(values = c("white", "white", "white")) +
    labs(x = NULL, y = paste0(geneName, " gene expression")) + 
    geom_label(label = paste0("P = ", pvalue), fill = "white", color = "black",x =3.2, y = outs[2] , label.size = 0, size = 7) +
    theme(legend.position = "none", 
          panel.background = element_rect(fill = NA, colour = "grey50"),  
          panel.grid = element_blank(),
          axis.text.x = element_text(size = 26, colour = "black"),
          axis.text.y = element_text(size = 26, colour = "black"),
          axis.title.y = element_text(size = 26, face = "bold",  vjust = 2))
  
  res <- list(boxPlot = boxPlot, vioPlot = vioPlot, kwp = pvalue, FCnt = FCnt, FCtm = FCtm, nNormal = nNorm, nTumor = nTum, nMeta = nMet)
 
  
}

