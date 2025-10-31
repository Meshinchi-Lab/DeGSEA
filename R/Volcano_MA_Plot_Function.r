#Jenny Smith

#3/5/18



#' function for volcano plot
#'
#' @param limma_eBayesFit a limma eBayesFit object
#' @param cut.off a numeric value of log2 fold-change for gene labels on points
#' @param label.offset a numeric value for nudge_y in ggrepel::geom_text_repel
#' @param contrast a numeric value for the contrast (column) from the coefficients to plot
#'
#' @returns a ggplot object
#' @export
#'
#' @examples
#' output <- c('TBD')
volcano_plot <- function(limma_eBayesFit, contrast = 1, cut.off=2.5, label.offset=0.5){


  df <- data.frame(Gene = names(limma_eBayesFit$coefficients[,contrast]),
                   logFC=limma_eBayesFit$coefficients[,contrast],
                   pValue=limma_eBayesFit$p.value[,contrast],
                   FDR=p.adjust(limma_eBayesFit$p.value[,contrast], method="BH"),
                   MeanExpression=limma_eBayesFit$Amean) %>%
    dplyr::mutate(Neg.Log10.P= -log10(pValue),
                  Neg.Log10.FDR= -log10(FDR),
                  gene_group=case_when(
                    logFC > 1.0 & FDR < 0.05 ~ "FC Greater than 2",
                    logFC < -1.0 & FDR < 0.05 ~ "FC Less than 2",
                    TRUE ~ "Not Significant FC")
    )


  #Select differentially expressed genes to highlight in the plot.
  idx <- which(abs(df$logFC) > cut.off & df$FDR < 0.05)


  vplot2 <- ggplot(df, aes(x=logFC, y=Neg.Log10.P)) +
    geom_point(data = filter(df, gene_group == "Not Significant FC"),
               mapping = aes(x=logFC, y=Neg.Log10.P, color=gene_group),
               alpha=0.5)  +

    geom_point(data= filter(df, grepl("2", gene_group)),
               mapping = aes(x=logFC, y=Neg.Log10.P, color=gene_group),
               alpha=0.75) +
    labs(y = "-log10(p-value)") +
    scale_color_manual(values=c("FC Greater than 2"="red",
                                "FC Less than 2"="deepskyblue",
                                "Not Significant FC"="lightgrey")) +

    theme(plot.title = element_text(hjust = 0.5, size = 14),
          panel.background = element_rect(fill="white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", fill=NA),
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(angle = 0,hjust=0.5,vjust = 0.5, size = 14),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(size = 16),
          plot.margin = margin(2,2,2,2, unit = "mm")) +
    ggrepel::geom_text_repel(aes(x=logFC, y=Neg.Log10.P, label=Gene),
                             size=3.5,
                             nudge_y = label.offset,
                             data=df[idx, ])


  return(vplot2)

}


#' MA plot
#'
#' @param limma_eBayesFit a limma eBayesFit object
#' @param cut.off a numeric value of log2 fold-change for gene labels on points
#' @param label.offset a numeric value for nudge_y in ggrepel::geom_text_repel
#'
#' @returns a ggplot object
#' @export
#'
#' @examples
#' output <- c('TBD')
MA_plot <- function(limma_eBayesFit, cut.off=2.5, label.offset=0.5){

  df <- data.frame(logFC=limma_eBayesFit$coefficients[,1],
                   pValue=limma_eBayesFit$p.value[,1],
                   FDR=p.adjust(limma_eBayesFit$p.value[,1], method="BH"),
                   MeanExpression=limma_eBayesFit$Amean) %>%
    tibble::rownames_to_column("Gene") %>%
    dplyr::mutate(Neg.Log10.P= -log10(pValue),
           gene_group=case_when(
             logFC > 1.0 & pValue < 0.05 ~ "FC Greater than 2",
             logFC < -1.0 & pValue < 0.05 ~ "FC Less than 2",
             TRUE ~ "Not Significant FC"))


    MAplot <- ggplot(df, aes(x=MeanExpression, y=logFC, color=gene_group)) +

      geom_point(data = filter(df, gene_group == "Not Significant FC"),
                 mapping = aes(x=MeanExpression, y=logFC, color=gene_group), alpha=0.65)  +

      geom_point(data= filter(df, grepl("2", gene_group)),
                 mapping = aes(x=MeanExpression, y=logFC, color=gene_group)) +

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


    #Select differentially expressed genes to highlight in the plot.
    idx <- which(abs(df$logFC) > cut.off & df$FDR < 0.05)

    MAplot <- MAplot +
      ggrepel::geom_text_repel(aes(x=logFC, y=Neg.Log10.FDR, label=Gene),
                             size=3.5,
                             nudge_y = label.offset,
                             data=df[idx, ])


    return(MAplot)
}















