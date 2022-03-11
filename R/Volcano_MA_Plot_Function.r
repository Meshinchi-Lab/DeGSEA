#Jenny Smith 

#3/5/18 


#purpose: function for volcano plots and MA plots 


volcano_plot <- function(twoGroups_DEGs.res, cut.off=2.5, label.offset=0.5){
  
  library(limma)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  
  #Only prints out - cannot be saved
  # vplot1 <- volcanoplot(twoGroups_DEGs.res$DE$eBayesFit, 
  #             highlight = 5, 
  #             names=rownames(twoGroups_DEGs.res$DE$eBayesFit),
  #             pch=16, cex=1.75)

  df <- data.frame(logFC=twoGroups_DEGs.res$DE$eBayesFit$coefficients[,1],
                   pValue=twoGroups_DEGs.res$DE$eBayesFit$p.value[,1],
                   FDR=p.adjust(twoGroups_DEGs.res$DE$eBayesFit$p.value[,1], method="BH"),
                   MeanExpression=twoGroups_DEGs.res$DE$eBayesFit$Amean) %>%
      rownames_to_column("Gene") %>%
      mutate(Neg.Log10.P= -log10(pValue),
             Neg.Log10.FDR= -log10(FDR),
             DEGs.Groups=case_when(
                  logFC > 1.0 & FDR < 0.05 ~ "FC Greater than 2",
                  logFC < -1.0 & FDR < 0.05 ~ "FC Less than 2",
                  TRUE ~ "Not Significant FC"))

  
  #Select differentially expressed genes to highlight in the plot. 
  idx <- which(abs(df$logFC) > cut.off & df$FDR < 0.05)
  
 
  vplot2 <- ggplot(df, aes(x=logFC, y=Neg.Log10.FDR)) + 
    geom_point(data = filter(df, DEGs.Groups == "Not Significant FC"), 
               mapping = aes(x=logFC, y=Neg.Log10.FDR, color=DEGs.Groups), alpha=0.65)  +
    
    geom_point(data= filter(df, grepl("2", DEGs.Groups)), 
               mapping = aes(x=logFC, y=Neg.Log10.FDR, color=DEGs.Groups)) +
    
    geom_vline(xintercept=c(-1,1)) +
    geom_hline(yintercept = -log10(0.05)) +
    
    scale_color_manual(values=c("FC Greater than 2"="red", 
                                "FC Less than 2"="blue",
                                "Not Significant FC"="lightgrey")) +
    
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          panel.background = element_rect(fill="white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", fill=NA),
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(angle = 0,hjust=0.5,vjust = 0.5, size = 26),
          axis.text.y = element_text(size = 25),
          axis.title = element_text(size = 30),
          plot.margin = margin(2,2,2,2, unit = "mm")) +
      
    geom_text_repel(aes(x=logFC, y=Neg.Log10.FDR, label=Gene),
                    size=3.5,
                    data=df[idx, ])
 
  res <- list("gg.volcano"=vplot2, "data"=df)
  return(res)
  
}



MA_plot <- function(twoGroups_DEGs.res, cut.off=2.5, label.offset=0.5){

  df <- data.frame(logFC=twoGroups_DEGs.res$DE$eBayesFit$coefficients[,1],
                   pValue=twoGroups_DEGs.res$DE$eBayesFit$p.value[,1],
                   FDR=p.adjust(twoGroups_DEGs.res$DE$eBayesFit$p.value[,1], method="BH"),
                   MeanExpression=DS.DE$DE$eBayesFit$Amean) %>%
    rownames_to_column("Gene") %>%
    mutate(Neg.Log10.P= -log10(pValue),
           DEGs.Groups=case_when(
             logFC > 1.0 & pValue < 0.05 ~ "FC Greater than 2",
             logFC < -1.0 & pValue < 0.05 ~ "FC Less than 2",
             TRUE ~ "Not Significant FC"))
  
  
    MAplot <- ggplot(df, aes(x=MeanExpression, y=logFC, color=DEGs.Groups)) + 
      
      geom_point(data = filter(df, DEGs.Groups == "Not Significant FC"), 
                 mapping = aes(x=MeanExpression, y=logFC, color=DEGs.Groups), alpha=0.65)  +
      
      geom_point(data= filter(df, grepl("2", DEGs.Groups)), 
                 mapping = aes(x=MeanExpression, y=logFC, color=DEGs.Groups)) +
      
      scale_color_manual(values=c("FC Greater than 2"="red", 
                                  "FC Less than 2"="blue",
                                  "Not Significant FC"="lightgrey")) + 
      
      theme(plot.title = element_text(hjust = 0.5, size = 20),
            panel.background = element_rect(fill="white"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "black", fill=NA),
            axis.text = element_text(color = "black"),
            axis.text.x = element_text(angle = 0,hjust=0.5,vjust = 0.5, size = 26),
            axis.text.y = element_text(size = 25),
            axis.title = element_text(size = 30),
            plot.margin = margin(2,2,2,2, unit = "mm"))
    
    return(MAplot)
}















