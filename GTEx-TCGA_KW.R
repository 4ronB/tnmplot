#this function servers to compare expression values of normal,tumor and metastatic RNA-Seq samples from GTEx and GDC databases using Kruskal-Wallis test
tnmk_GTEx_TCGA <- function(tissueID, geneName){
  
  # open co-normalized expression table 
  expWD <- paste0("TNMdata/GTEx_Expression-Tables/",tissueID[[1]],"-raw-clean-norm.RData")
  load(expWD)
  expTable <- txtFile
  rm(txtFile)
  
  #opening clinical table
  clinical_TCGA <- read_xlsx(path = "GDC_clinical.xlsx")
  clinical_TCGA <- as.data.frame(clinical_TCGA)
  
  #screening for selected tissue in clinical table
  TCGAID <- unlist(strsplit(tissueID[[1]], split = ".", fixed = T))[2]
  usedTissue <- clinical_TCGA[grepl(pattern = TCGAID, clinical_TCGA$`Project ID`),]
  sampNames <- usedTissue[,1:4]
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
  colnames(sampNames)[1] <- "SampleID"
  expValue <- as.data.frame(t(expTable[expTable$geneid == geneName,]))
  expValue <- as.data.frame(expValue[2:nrow(expValue),])
  colnames(expValue)[1] <- "exp"
  expValue$SampleID <- row.names(expValue)
  
  gtexData <-  data.frame("SampleID" = row.names(expValue)[grep(pattern = "GTEX", x = row.names(expValue))], "type" = "Normal")
  sampNames <- rbind(sampNames, gtexData)
  
  sampExp <- sampNames
  sampExp <- merge(x = sampExp, y = expValue, by = "SampleID")
  colnames(sampExp) <- c("SampID", "type", "Geneexp")
  sampExp$type <- factor(x = sampExp$type, levels =  c("Normal", "Tumor", "Metastatic"))
  
  #removing unnecessary data
  rm(list = c("clinical_TCGA", "expTable", "expValue", "sampNames", "usedTissue", "i", "gtexData"))
  
  ###################################statistical analysis##############################
  sampExp[,3] <- as.numeric(as.vector(sampExp[,3]))
  
  nNorm <- nrow(sampExp[sampExp$type == "Normal",])
  nTum <- nrow(sampExp[sampExp$type == "Tumor",])
  nMet <- nrow(sampExp[sampExp$type == "Metastatic",])
  nAll <- nNorm + nTum + nMet
  
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
  
  library(ggplot2)
  boxPlot <- ggplot(sampExp, aes(x = type, y = Geneexp, color = type, fill = type))  + 
    geom_boxplot(outlier.size = -1, notch = F, width = 0.55) + 
    ylim(outs[1], outs[2]) +
    geom_jitter(width = 0.25) + 
    scale_color_manual( values =  c("darkgreen", "red", "darkorange3")) + 
    scale_fill_manual(values =  c("darkolivegreen1", "coral1", "goldenrod1")) +
    labs(x = NULL, y = paste0(geneName, " gene expression")) + 
    geom_label(label = paste0("P = ", pvalue), fill = "white", color = "black",x =3.2, y = outs[2] , label.size = 0, size = 7) +
    theme(legend.position = "none", 
          panel.background = element_rect(fill = NA, colour = "grey50"),  
          panel.grid = element_blank(),
          axis.text.x = element_text(size = 26, colour = "black"),
          axis.text.y = element_text(size = 26, colour = "black"),
          axis.title.y = element_text(size = 27, face = "bold",  vjust = 2))
  
  
  vioPlot <- ggplot(sampExp, aes(x = type, y = Geneexp, color = type, fill = type))  + 
    geom_violin(lwd = 1) + ylim(outs[1], outs[2]) +
    geom_boxplot(lwd = 1.2, outlier.size = -1, width = 0.1, color = c("darkgreen", "red", "darkorange3")) +
    scale_color_manual( values =  c("darkgreen", "red", "darkorange3")) + 
    scale_fill_manual(values = c("white", "white", "white"))+
    labs(x = NULL, y = paste0(geneName, " gene expression")) + 
    geom_label(label = paste0("P = ", pvalue), fill = "white", color = "black",x =3.2, y = outs[2] , label.size = 0, size = 7) +
    theme(legend.position = "none", 
          panel.background = element_rect(fill = NA, colour = "grey50"),  
          panel.grid = element_blank(),
          axis.text.x = element_text(size = 26, colour = "black"),
          axis.text.y = element_text(size = 26, colour = "black"),
          axis.title.y = element_text(size = 26, face = "bold",  vjust = 2))
  
  res <- list(boxPlot = boxPlot, vioPlot = vioPlot, kwp = pvalue, FCnt = FCnt, FCtm = FCtm, nNormal = nNorm, nTumor = nTum, nMeta = nMet, DpValues = dunnPvalues)
  
}
