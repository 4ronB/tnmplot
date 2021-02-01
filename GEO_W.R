#this function servers to compare expression values of paired normal and tumor gene chip samples using Wilcoxon test
tnm_geo_w <- function(tissueID, geneName){
  
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
  sampNames <- usedTissue[!usedTissue$Pair == 0,]#screening for paired data
  sampNames <- sampNames[!sampNames$Meta == 1,]
  norms <-  as.vector(unlist(sampNames$Pair))
  normSamps <- usedTissue[usedTissue$GSM %in% norms,]
  sampNames <- merge(sampNames, normSamps, all = T, sort = F)
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
  sampNames <-  sampNames[,-(2:9)]
  
  #creating expression table for the selected gene
  expvalue <- as.data.frame(t(expTable[expTable$gene_sym == geneName, colnames(expTable) %in% sampNames$GSM]))
  expvalue$GSM <- rownames(expvalue) 
  rownames(expvalue) <- c()
  #merging the clinical data and expression table - this data frame will be used for further investigations
  sampExp <- sampNames
  sampExp <- merge(x = sampExp, y = expvalue, by = "GSM", sort = F)
  colnames(sampExp)[3] <- "Geneexp"
    
  ###################################statistical analysis##############################
  
  nNorm <- nrow(sampExp[sampExp$type == "Normal",])
  nTum <- nrow(sampExp[sampExp$type == "Tumor",])
  nAll <- nNorm + nTum
  
  meanNorm <- mean(sampExp[sampExp$type == "Normal",3])
  sdnorm <- sd(sampExp[sampExp$type == "Normal",3])
  meanTum <- mean(sampExp[sampExp$type == "Tumor",3])
  sdtum <- sd(sampExp[sampExp$type == "Tumor",3])
  
  FC <- meanTum/meanNorm
  
  mwt <- wilcox.test(x = sampExp[sampExp$type == "Normal",3],
                     y = sampExp[sampExp$type == "Tumor",3],
                     paired = T,
                     exact = F,
                     conf.int = T)
  pvalue <- format(mwt$p.value, digits = 3, scientific = T)
  
  
  #identification of outliers
  normOuts <- boxplot.stats(sampExp[sampExp$type == "Normal",3])$stats
  tumOuts <- boxplot.stats(sampExp[sampExp$type == "Tumor",3])$stats
  outs <- range(normOuts, tumOuts)
  
  #barchart data
  
  tumExp <- sampExp[sampExp$type == "Tumor", 3]
  normExp <- sampExp[sampExp$type == "Normal", 3]
  
  Maxpercent <- round((length(tumExp[tumExp > max(normExp)])/nTum)*100, 1)
  Q3percent <- round((length(tumExp[tumExp > normOuts[4]])/nTum)*100, 1)
  Medpercent <- round((length(tumExp[tumExp > normOuts[3]])/nTum)*100, 1)
  Q1percent <- round((length(tumExp[tumExp > normOuts[2]])/nTum)*100, 1)
  Minpercent <- round((length(tumExp[tumExp > normOuts[1]])/nTum)*100, 1)
  
  #preparing values for specificity examination
  if(FC>1){
    sMax <- round((length(tumExp[tumExp > max(normExp)])/(length(normExp[normExp > max(normExp)]) + length(tumExp[tumExp > max(normExp)])))*100, 1)
  } else if(FC<1){
    sMax <- round((length(tumExp[tumExp < max(normExp)])/(length(normExp[normExp < max(normExp)]) + length(tumExp[tumExp < max(normExp)])))*100, 1)
  }
  if(FC>1){
    sQ3 <- round((length(tumExp[tumExp > normOuts[4]])/(length(normExp[normExp > normOuts[4]]) + length(tumExp[tumExp > normOuts[4]])))*100, 1)
  } else if(FC<1){
    sQ3 <- round((length(tumExp[tumExp < normOuts[4]])/(length(normExp[normExp < normOuts[4]]) + length(tumExp[tumExp < normOuts[4]])))*100, 1)
  }
  if(FC>1){
    sMed <- round((length(tumExp[tumExp > normOuts[3]])/(length(normExp[normExp > normOuts[3]]) + length(tumExp[tumExp > normOuts[3]])))*100, 1)
  } else if(FC<1){
    sMed <- round((length(tumExp[tumExp < normOuts[3]])/(length(normExp[normExp < normOuts[3]]) + length(tumExp[tumExp < normOuts[3]])))*100, 1)
  }
  if(FC>1){
    sQ1 <- round((length(tumExp[tumExp > normOuts[2]])/(length(normExp[normExp > normOuts[2]]) + length(tumExp[tumExp > normOuts[2]])))*100, 1)
  } else if(FC<1){
    sQ1 <- round((length(tumExp[tumExp < normOuts[2]])/(length(normExp[normExp < normOuts[2]]) + length(tumExp[tumExp < normOuts[2]])))*100, 1)
  }
  if(FC>1){
    sMin <- round((length(tumExp[tumExp > normOuts[1]])/(length(normExp[normExp > normOuts[1]]) + length(tumExp[tumExp > normOuts[1]])))*100, 1)
  } else if(FC<1){
    sMin <- round((length(tumExp[tumExp < normOuts[1]])/(length(normExp[normExp < normOuts[1]]) + length(tumExp[tumExp < normOuts[1]])))*100, 1)
  }
  
  barchart <- data.frame("Q" = c("Max", "Q3", "Med", "Q1", "Min") ,
                         "perc" = c(Maxpercent, Q3percent, Medpercent, Q1percent, Minpercent), 
                         "spec" = c(sMax, sQ3, sMed, sQ1, sMin))
  barchart$Q <- factor(barchart$Q, levels = rev(barchart$Q))
  
  
  ##############################data visualization with boxplots######################
  
  library(ggplot2)
  boxPlot <- ggplot(sampExp, aes(x = type, y = Geneexp, color = type, fill = type))  + 
    geom_boxplot(outlier.size = -1, notch = F, width = 0.55) + ylim(outs[1], outs[2]) +
    geom_jitter(width = 0.25) + 
    scale_color_manual( values =  c("darkgreen", "red")) + 
    scale_fill_manual(values =  c("darkolivegreen1", "coral1")) +
    labs(x = NULL, y = paste0(geneName, " gene expression")) + 
    geom_label(label = paste0("P = ", pvalue), fill = "white", color = "black",x =2.32, y = outs[2] , label.size = 0, size = 7) +
    theme(legend.position = "none", 
          panel.background = element_rect(fill = NA, colour = "grey50"),  
          panel.grid = element_blank(),
          axis.text.x = element_text(size = 26, colour = "black"),
          axis.text.y = element_text(size = 26, colour = "black"),
          axis.title.y = element_text(size = 26, face = "bold",  vjust = 2))
  
  
  vioPlot <- ggplot(sampExp, aes(x = type, y = Geneexp, color = type, fill = type))  + 
    geom_violin(lwd = 1) + ylim(outs[1], outs[2]) +
    geom_boxplot(lwd = 1.2, outlier.size = -1, width = 0.1, color = c("darkgreen", "red")) +
    scale_color_manual( values =  c("darkgreen", "red")) + 
    scale_fill_manual(values = c("white", "white")) +
    labs(x = NULL, y = paste0(geneName, " gene expression")) + 
    geom_label(label = paste0("P = ", pvalue), fill = "white", color = "black",x =2.32, y = outs[2] , label.size = 0, size = 7) +
    theme(legend.position = "none", 
          panel.background = element_rect(fill = NA, colour = "grey50"),  
          panel.grid = element_blank(),
          axis.text.x = element_text(size = 26, colour = "black"),
          axis.text.y = element_text(size = 26, colour = "black"),
          axis.title.y = element_text(size = 26, face = "bold",  vjust = 2))
  
  barPlot <- if (FC > 1){
    ggplot(barchart) +
      geom_col(aes(x = Q, y = perc, fill = Q), width = 0.9) +
      geom_line(aes(x = Q, y = spec), size = 1, color="red", group = 1) +
      geom_point(aes(x = Q, y = spec), shape = 21, color = "red", fill = "red", size = 6) +
      ylim(c(0,100)) +
      xlab("Cut-off in normal tissues") + ylab("% of tumors higher than normal") +
      scale_y_continuous(sec.axis = sec_axis(~ . /100, name = "Specificity (T/T+N over cutoff)")) +
      scale_fill_manual(values =  c("darkslategray1", "darkslategray2",  "cyan2",  "darkslategray3", "darkslategray4")) +
      geom_text(aes(x = Q, y = perc -2.2, label = paste0(perc, "%")), size = 9)+
      theme(legend.position = "none", 
            panel.background = element_rect(fill = NA, colour = "grey50"),  
            panel.grid = element_blank(),
            axis.text.x = element_text(size = 25, colour = "black"),
            axis.text.y = element_text(size = 25, colour = "blue"),
            axis.title.x = element_text(size = 25, face = "bold",  vjust = 2, margin = margin(t = 25, r = 0, b = 0, l = 0)),
            axis.title.y = element_text(size = 25, face = "bold",  vjust = 2,, colour = "blue"),
            axis.text.y.right = element_text(size = 26, colour = "red"),
            axis.title.y.right = element_text(size = 25, face = "bold",  vjust = 2, colour = "red"))
  } else if (FC <= 1){
    ggplot(barchart) +
      geom_col(aes(x = Q, y = perc, fill = Q), width = 0.9) +
      geom_line(aes(x = Q, y = spec), size = 1, color="red", group = 1) +
      geom_point(aes(x = Q, y = spec), shape = 21, color = "red", fill = "red", size = 6) +
      ylim(c(0,100)) +
      xlab("Cut-off in normal tissues") + ylab("% of tumors over normal") +
      scale_y_continuous(sec.axis = sec_axis(~ . /100, name = "Specificity (T/T+N below cutoff)")) +
      scale_fill_manual(values =  c("darkslategray1", "darkslategray2",  "cyan2",  "darkslategray3", "darkslategray4")) +
      geom_text(aes(x = Q, y = perc -2.2, label = paste0(perc, "%")), size = 9)+
      theme(legend.position = "none", 
            panel.background = element_rect(fill = NA, colour = "grey50"),  
            panel.grid = element_blank(),
            axis.text.x = element_text(size = 25, colour = "black"),
            axis.text.y = element_text(size = 25, colour = "blue"),
            axis.title.x = element_text(size = 25, face = "bold",  vjust = 2, margin = margin(t = 25, r = 0, b = 0, l = 0)),
            axis.title.y = element_text(size = 25, face = "bold",  vjust = 2,, colour = "blue"),
            axis.text.y.right = element_text(size = 26, colour = "red"),
            axis.title.y.right = element_text(size = 25, face = "bold",  vjust = 2, colour = "red"))
  }
  res <- list(boxPlot = boxPlot, vioPlot = vioPlot, barPlot = barPlot, mwp = pvalue, FC = FC, nNormal = nNorm, nTumor = nTum)
  
  
}

