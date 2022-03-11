#Jenny Smith 


#Feb. 20, 2018

#Puprpose: to create a correlation plot with samples, and provide a color bar for cytogenetics

get.CPM.Df <- function(twoGroups_DEGs.res, trend=TRUE){
  
  if(trend){
    CPM.DEGs <- twoGroups_DEGs.res$DE$dge %>%
      data.frame() %>%
      rownames_to_column("gene") %>%
      filter(gene %in% rownames(twoGroups_DEGs.res$DE$DE)) %>%
      column_to_rownames("gene")
  }else{
    CPM.DEGs <- twoGroups_DEGs.res$DE$Voom$E %>%
      data.frame() %>%
      rownames_to_column("gene") %>%
      filter(gene %in% rownames(twoGroups_DEGs.res$DE$DE)) %>%
      column_to_rownames("gene")
  }

  return(CPM.DEGs)
}


corPlot_withColorBar <- function(expnData,CDE,col2=NULL, meanctr=TRUE,log2=TRUE,
                                 lab.size=0,diag.point.size=15, offset=15){
  #expnData is normalized expresion data, log2 scale, with genes as rownames
  #CDE is the clinical data with a column called,"Status" and "Num.Status" for the factor levels of each cyto group as factor/and numeric,
  #and CDE has a column called "USI"for merging. 
  #col2 is an optional column 
  #data will be mean centered (geometric mean on log2 scale) in order to visualize the sample correlations best. 
  
  #lab.size is a single numeric.  if > 0 will indicate to keep patient labels on the x and y axis. 
  #diag.point.size is the size of the point for the patient CDE groups.
  #offset a single  numeric on the x axis to move the points for the diaganol. 
  
  library(dplyr)
  library(RColorBrewer)
  library(cowplot)
  library(ggplot2)
  library(dplyr)
  library(reshape2)
  
  #Function to reorder the cor matrix using correlation as distance and ward.D2 clustering. 
  reorder_cormat <- function(cormat,vector=FALSE){
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd, method="ward.D2")
    cormat <-cormat[hc$order, hc$order]
    if (vector){
      order <- rownames(cormat)
      return(order)
    }else{
      return(cormat)
    }}
  
  if(!log2){
    expnData <- log2(expnData+1)
  }
  
  #this is geometric mean centering because the antilog of the arithmatic mean of logged values == geomtric mean of non-logged values 
  #mean center the data/scale to make z-scores. 
  expnData.ctr <- scale(t(expnData), center = TRUE, scale = FALSE)
  expnData.ctr <- t(expnData.ctr)
  cormat.ctr <- cor(expnData.ctr)
  
  #find the order of the correlation matrix using hclust with ward.D2 algorithm. 
  cormat.order <- reorder_cormat(cormat.ctr, vector=TRUE)
  cormat.final <- reorder_cormat(cormat.ctr)
  cormat.final[upper.tri(cormat.final)] <- NA
  
  
  #Melt the correlation matrix and add the clinical data
  cormat.melt <- cormat.final %>%
    melt(.,na.rm=TRUE, as.is=TRUE, varnames=c("Var1","Var2")) %>%
    # melt(cormat.final, na.rm = TRUE)  #original code. The behavior of melt() changed? Now defauls to X1 and X2.
    inner_join(., CDE, by=c("Var1"="USI")) %>%
    mutate_at(c("Var1", "Var2"),funs(factor(., levels=cormat.order))) %>% #ensure factor levels are correct before plotting
    mutate(Num.USI=as.numeric(Var1),
           Num.pair=as.numeric(Var2)) %>% #convert to numeric ordering
    select(Var1,Var2,Status,Num.USI,Num.pair,Num.Status, everything())
  
  #Create a subset of the cormat.melt to contain just the diaganol. 
  diag <- subset(cormat.melt, Num.USI==Num.pair)
  n.samp <- nrow(diag)
  
  #Use geom_raster to create a correlation heatmap. 
  ggheatmap <- ggplot(cormat.melt, aes(x=Num.USI, y=Num.pair, fill = value))+
    geom_raster() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, 
                         limit = c(-1,1), space = "Lab",   name="Pearson\nCorrelation") +
    theme_minimal() + # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1, size = lab.size, hjust = 1),
          axis.text.y = element_text(angle = 45, vjust = 1, size = lab.size, hjust = 1),
          panel.border = element_blank(),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_blank()) +
    coord_fixed(ratio = 1, xlim = c(1,n.samp), ylim=c(1,n.samp))

  if(n.samp < 50){
    heatmap2 <-  ggheatmap +
      geom_point(data = diag,
                 mapping = aes(x=Num.USI-offset, y=Num.pair, color=Num.Status), size=diag.point.size, shape=15, show.legend = FALSE) 
  }else{
    heatmap2 <-  ggheatmap +
      geom_jitter(data = diag,
                  mapping = aes(x=Num.USI-offset, y=Num.pair, color=Num.Status), size=diag.point.size, shape=95, show.legend = FALSE) #shape 124 is a line or pipe shape, 95 is a dash
  }

  
  if (!is.null(col2)){
    #Note: this has not been tested. 
    offset2 <- offset*3
    
    heatmap2 <- heatmap2 + 
      geom_jitter(data = diag,
                  mapping = aes_(x=Num.USI-offset2, y=Num.pair, color=as.name(col2)), 
                  size=15, shape=95, show.legend = FALSE)
    
  }
  
  
  n <- length(unique(CDE$Num.Status))
  
  if( n <= 2){
    col <- c("deepskyblue","navy")
  }else{
    pal <- brewer.pal(min(n,9), "Set1") %>% 
      gsub("#F781BF", "burlywood2", .) %>%
      c(., "turquoise3", "blue", "black","seagreen2", "maroon", 
        "orchid", "cornflowerblue",  "darkblue", "azure4", "chartreuse1", 
        "darkmagenta","orange1", "deeppink", "darkslategray1",
        "green4", "navajowhite2","brown3", "darkgoldenrod3", "deepskyblue1", 
        "lightcoral", "mediumorchid", "saddlebrown")
    col <- pal[1:n]
  }
  
  heatmap2 <- heatmap2 + 
      scale_color_gradientn(colors = col )  +
      labs(x="AML Samples", y="AML Samples") +
      theme(axis.title  = element_text(size=24),
            # axis.text.x = element_blank(),
            # axis.text.y = element_blank(),
            legend.key.size = unit(30, units = "points"),
            legend.text = element_text(size=22),
            legend.title = element_text(size=24)) 
  
  #If the label point size is not zero, add the patient labels. 
  if (lab.size > 0){
    heatmap2<- heatmap2 + 
      scale_x_reverse(breaks=seq(1,n.samp,by=1), labels = diag$Var1) +
      scale_y_continuous(breaks=seq(1,n.samp,by=1), labels = diag$Var1)
  }else{
    heatmap2 <- heatmap2 +
      scale_x_reverse()
  }
  
  
  #create a legend for the heatmap
  subset <- CDE %>%
    select(Status, Num.Status) %>%
    arrange(Num.Status) %>%
    unique() 
  
  color.Legend <- ggplot(subset, aes(x=Status, y=0.25, fill=Status)) + 
    geom_bar(stat = "identity") + 
    scale_fill_manual(values=col) +
    theme(legend.text = element_text(size=15))
  
  l <- cowplot::get_legend(color.Legend)
  
    
  list <- list(cormat.melt,  ggheatmap, heatmap2,l)
  names(list) <- c("cormat.melt",  "ggheatmap", "heatmap2", "legend")
  return(list)
}




#Code for the selection of highly correlated pairs of genes
maxCorrPairs <- function(expnData,cutoff){
  library(dplyr)
  library(psych)
  #expnData is genes as rows and patients as columns. must be a data frame class 
  
  #find absolute correlation values 
  gene.cor.abs <- abs(cor(t(expnData)))
  p <- corr.test(t(expnData), adjust = "BH",ci = FALSE)
  
  #function to pullout the highest genes correlation pair
  maxAbsCor <- function(x){
    x <- x[x < 1] #remove the correlation to itself
    x <- x[x == max(x)]
    res <- data.frame(cor=x,
                      pair=names(x), 
                      stringsAsFactors = FALSE)
    
    return(res)
  }
  
  getPval <- function(gene, pair, pvalues){
    p <- pvalues[gene,pair]
    return(p)
  }
  
  #apply the function to each row or column (doesn't matter its a diagnol upper and lower matrix still)
  maxcors <- bind_rows(apply(gene.cor.abs,2, maxAbsCor)) %>%
    mutate(gene=colnames(gene.cor.abs)) %>%
    rowwise() %>%
    mutate(pval.BH=getPval(gene,pair,p$p)) %>%
    ungroup() %>%
    filter(!duplicated(cor)) %>%
    arrange(desc(cor)) %>%
    filter(cor >= quantile(cor, c(cutoff)) & pval.BH < 0.01) #select correlations in the 80th percentile or greater
  
  
  subset <- expnData[unique(c(maxcors$pair, maxcors$gene)),]
  
  list <- list(gene.cor.abs, maxcors, subset)
  names(list) <- c("gene.cor.abs", "maxcors","subset.expn")
  
  return(list)
  
}














