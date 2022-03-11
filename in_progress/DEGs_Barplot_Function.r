#Jenny Smith

#July 11, 2017 

#purpose: plot the highest genes fold-changes

barplot <- function(limma.res, decile=NULL) {
  
  dec <- quantile(limma.res$logFC, probs = seq(0,1,length= 11), type=5)
  limma.res$gene <- rownames(limma.res)
  
  if (is.null(decile)){
    limma.res <- limma.res
  }else{
    decile.90 <- dec[10]
    decile.10 <- dec[2]
    limma.res <- limma.res[which(limma.res$logFC >= decile.90 | limma.res$logFC <= decile.10), ]
  }
  limma.res$cut <- ifelse(limma.res$logFC > 0, "up", "dn")
  
  ggplot(limma.res, aes(x=reorder(limma.res$gene, limma.res$logFC), y=logFC, fill=cut)) + 
    geom_bar(stat = "identity") + theme_JS + labs(title="Differentially Regulated Genes: Largest FC",x="") +
    scale_fill_manual(values = c(up="red", dn="darkgreen"))
}


