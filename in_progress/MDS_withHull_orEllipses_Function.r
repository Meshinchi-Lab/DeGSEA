#Jenny Smith 


#Feb 28, 2018 


#purpose: to create a classical multidimensional scaling plot with convex hulls based on Kmeans or kmedians. 


MDSplot.Hull <- function(df,cols, method="kmeans",ellipseType="convex", k){
  #df has the clinical data and expression data. 
  #patients as rows and genes as columns 
  #enusre ONLY the gene expression data is numeric!! 
  #cols is a character vector of colnames for the clinical data to be used for fills or shapes 
  
  library(dplyr)
  library(vegan)
  library(flexclust)
  library(ggpubr)
  
  MDS <- capscale(select_if(df,is.numeric) ~ 1, distance = "bray", add=TRUE)
  
  scores.mds <- scores(MDS, display="sites") %>%
    data.frame() 
  
  if (method == "kmeans"){
    kmean <- kmeans(scores[,1:2],k)
    groups <- as.factor(kmean$cluster)
    
  }else if (method == "kmedians"){
    require(flexclust)
    kmed <- flexclust::kcca(scores[,1:2],k=k, family = kccaFamily("kmedians"), save.data=TRUE)
    groups <- as.factor(kmed@cluster)
  }
  
  scores <- scores %>%
    mutate(k.groups=groups) %>%
    bind_cols(df[,cols]) #same order as in the input df to the capscale() function. 
  
  
  mds.plot <- ggscatter(scores, x="MDS1", y="MDS2", 
                        size=3.0, 
                        color=cols[1], #just use the first column by default, will change later to have more than 1 variable 
                        palette = "jco",
                        ellipse =FALSE,
                        ellipse.type = "norm",
                        repel = TRUE,
                        add.params = list(alpha=0.5))
  
  
  if (ellipseType=="convex"){
    find_hull <- function(df) df[chull(df$MDS1, df$MDS2), ]
    hulls <- plyr::ddply(scores,"k.groups",find_hull)
    
    mds.plot <- mds.plot + 
      geom_polygon(data=hulls, aes(fill=k.groups), alpha=0.2) +
      theme(text = element_text(size=20)) +
      labs(x="MDS1", y="MDS2")
    
  }else if(ellipseType=="t"){
    mds.plot <- mds.plot + 
      stat_ellipse(data=scores, mapping = aes(group=k.groups,fill=k.groups), 
                   geom="polygon", alpha=0.15, color="black", level=0.95) +
      theme(text = element_text(size=20)) +
      labs(x="MDS1", y="MDS2")
    
  }
  
  res <- list(MDS, scores, mds.plot)
  names(res) <- c("mds.obj","mds.scores","mds.plot")
  
  return(res)
}


