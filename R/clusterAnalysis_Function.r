#Jenny Smith

#May 8, 2017

#Purpose: Create a PCA, and MDS ordination plots given expression data and genes of interest.


#' Multidimensional Scaling plot
#'
#' @param expnData normalized counts, typically log2 scale.
#' @param phenovector named character vector of patient ID = groupID
#' @param title plot title
#' @param colorCodes named vector like groupID = 'red'
#' @param geneList character vector gene_names or gene_ids
#' @param size point size for the plot
#'
#' @return list
#' @export
#'
#' @examples
#' ex <- c('TBD')
plotPCoA <- function(expnData,phenovector, title="",colorCodes=NULL, geneList=NULL, size=3){
  #expnData is normalized counts, typically log2 scale.

  #factor is the name of the factor column
  # suppressPackageStartupMessages(library(vegan))
  # suppressPackageStartupMessages(library(ggplot2))
  # suppressPackageStartupMessages(library(dplyr))

  #Ensure correct order of patients in both datasets
  expnData <- expnData[ ,intersect(names(phenovector), colnames(expnData))]
  phenovector <- phenovector[intersect(names(phenovector),colnames(expnData))]

  if (! is.null(geneList)){
    expnData <- t(expnData[geneList, ])
  }else{
    expnData <- t(expnData) #note: must remove all zero count genes or will  fail on an error
  }

  PCoA <- capscale(expnData ~ 1, distance = "bray", add=TRUE)
  scores <- data.frame(scores(PCoA, display="sites"),
                       group=phenovector) %>%
            dplyr::arrange(desc(group))


  p <- ggplot(scores, aes(x=MDS1, MDS2)) +
    geom_point(aes(color=scores[,"group"]), size=size, alpha=0.75) +
    theme_numX +
    labs(title=title)

  if(!is.null(colorCodes)){
    p <- p +
      scale_color_manual(values=colorCodes)
  }
  #If wanted to add ellipses and convex hulls use the code below

  # k <- 4
  # clust <- kmeans(scores[,1:2], k)
  # groups <- as.factor(clust$cluster)
  #
  # scores <- scores %>%
  #   mutate(group=groups,
  #          Status=clinData[rownames(expnData),factor])
  # head(scores)
  # mds.plot <- ggscatter(scores, x="MDS1", y="MDS2",
  #                       size=3.5,
  #                       color="Status",
  #                       group="group",
  #                       palette = "Set1",
  #                       ellipse =FALSE,
  #                       ellipse.type = "norm",
  #                       repel = TRUE)
  # #https://stats.stackexchange.com/questions/22805/how-to-draw-neat-polygons-around-scatterplot-regions-in-ggplot2/22855
  # find_hull <- function(df) df[chull(df$MDS1, df$MDS2), ]
  # hulls <- ddply(scores,"group",find_hull )
  #
  # #add hulls or ellipses
  # mds.plot <- mds.plot +
  #   # stat_ellipse(data = scores,
  #   # mapping = aes(group=group, fill=group), geom = "polygon", alpha=0.15, color="black", level=0.99)
  #   geom_polygon(data=hulls, aes(fill=group), alpha=0.15) +
  #   scale_fill_manual(values = rainbow_hcl(4)) +
  #   theme(text = element_text(size=20))

  # aov <- aov(scores$MDS1 ~ df[,factor]) #NOTE: This is only valid for balanced experimental designs! Equal # of obs in each factor level.

  list <- list(PCoA,scores, p)
  names(list) <- c("PCoA","scores","plot")

  return(list)
}


#from https://github.com/mikelove/DESeq2/blob/master/R/plots.R
#Want to return the whole scores matrix so can examine 3d pca plots.
plotPCA.DESeq.mod <- function(object, intgroup="condition", ntop=500, returnData=FALSE, PC3=FALSE)
{
  # library(matrixStats)
  # calculate the variance for each gene
  rv <- rowVars(assay(object))

  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))

  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }

  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])

  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    colData(object)[[intgroup]]
  }

  # assembly the data for the plot - first 10 PCs
  # d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=group, intgroup.df, name=colnames(object))
  n <- min(10, ncol(as.data.frame(pca$x)))
  d <- data.frame(as.data.frame(pca$x)[,1:n], group=group, intgroup.df, name=colnames(object))

  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:10]
    rot <- pca$rotation[,1:10] #for first 10 PCs
    dat <- list("scores"=d,"rotation"=rot)
    return(dat)
  }

  if(PC3){
    ggplot(data=d, ggplot2::aes_string(x="PC1", y="PC3", color="group")) +
      geom_point(size=3, alpha=0.75) +
      xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
      ylab(paste0("PC3: ",round(percentVar[3] * 100),"% variance"))

  }else{
    ggplot(data=d, ggplot2::aes_string(x="PC1", y="PC2", color="group")) +
      geom_point(size=3, alpha=0.75) +
      xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
      ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance"))
    # coord_fixed()
  }
}




#' Generate PCA plot
#'
#' @param expnData raw counts (not normalized) has patient IDs as colnames and genes as rownames.
#' @param phenovector named character vector with patient ID = groupID
#' @param title character vector
#' @param round boolean
#' @param colorCodes named vector like c(groupID = 'red')
#' @param ntop number of most varied genes to use in PCA
#' @param PC3 plot PC1 vs PC3 as well as PC1 vs PC2
#' @param GOI character vector of gene names/ gene id that match the expnData.
#'
#' @return list
#' @export
#'
#' @examples
#' ex <- c('TBD')
PCA <- function(expnData,phenovector,title="",round=TRUE,colorCodes=NULL,
                ntop=500,PC3=FALSE, GOI=NULL){

  # suppressPackageStartupMessages(library(DESeq2))
  # library(ggplot2)
  #expnData is the raw counts (not normalized) has patient IDs as colnames and genes as rownames.

  # countData <- expnData[,match(names(phenovector), colnames(expnData))]
  samples <- intersect(names(phenovector), colnames(expnData))
  countData <- expnData[,samples]
  phenovector <- phenovector[samples]

  countData <- round(countData, digits = 0)
  colData <- as.data.frame(phenovector)

  if(length(unique(phenovector)) < 2){
    dds <- DESeqDataSetFromMatrix(countData = countData,
                                  colData = colData,
                                  design = ~ 1)
  }else{
    dds <- DESeqDataSetFromMatrix(countData = countData,
                                  colData = colData,
                                  design = ~ phenovector)
  }

  #create deseq2 dataset and perform variance stabilized transformation for sample to sample comparisons
  dds <- dds[ rowSums(counts(dds)) > 10, ]
  varianceStab <- vst(dds, blind = TRUE)

  #if given a list of genes of interest
  if (! is.null(GOI)){
    GOI <- intersect(GOI, rownames(assay(varianceStab)))
  }else{
    GOI <- 1:nrow(varianceStab)
  }

  # Create a PCA plot
  plot.1 <- plotPCA.DESeq.mod(varianceStab[GOI,], intgroup = "phenovector", ntop = ntop, PC3=FALSE) +
    theme_numX +
    labs(title=title)

  if(PC3){
    plot.2 <- plotPCA.DESeq.mod(varianceStab[GOI,], intgroup = "phenovector", ntop = ntop, PC3=TRUE) +
      theme_numX +
      labs(title=title)
  }

  #change points to custom colors if provided
  if (!is.null(colorCodes)){
    plot.1 <- plot.1 +
      scale_color_manual(values=colorCodes)

    if(exists("plot.2")){
      plot.2 <- plot.2 +
        scale_color_manual(values=colorCodes)
    }

  }

  #PCA data frame with the wieghts/loadings and eigen vectors
  pca.dat <- plotPCA.DESeq.mod(varianceStab[GOI,], intgroup = "phenovector", ntop = ntop,
                               returnData=TRUE)

  #Final Results object
  res <- list(dds, varianceStab,pca.dat$scores,pca.dat$rotation, plot.1)
  names(res) <- c("dds", "vst","pca_data","pca_loadings", "pca_plot")

  if(is.character(GOI)){
    res[["GOI"]] <- GOI
  }

  if(PC3){
    res[["pca_plot2"]] <- plot.2
  }

  return(res)
}


#' Custom PCA plot
#'
#' @param expnData expn data is log2, normalized counts
#' @param CDE dataframe - has patients as rownames
#' @param fillCol character string of column name for fill colors
#' @param colorCol character string of column name for border colors
#' @param colorCode an option named vector of colors
#' @param PC3 boolean to generate PC1 vs PC3 plot
#' @param single.col.outline boolean whether to susbet to dataframe for the border color
#' @param toHighlight if single.col.outline is true, Gives the factor to highlight with borders from the colorColumn
#' @param ellipse boolean to draw an ellipse based on fillCol column
#'
#' @returns list of plots
#' @export
#'
#' @examples
#' mat <- sapply(1:10, function(i){ runif(100,min = 0, max = 15) })
#' colnames(mat) <- paste0("s",1:10)
#' metadata <- data.frame(sample_id = paste0("s",1:10), mutation = rep("yes","no", length.out = 10))
#' rownames(metadata) <- paste0("s",1:10)
#' pca_custom(expnData = mat,  CDE = metadata, fillCol = "mutation", colorCol = "mutation")
pca_custom <- function(expnData,CDE,fillCol, colorCol, colorCode=NULL, PC3=FALSE,
                       single.col.outline=FALSE, toHighlight=NULL, ellipse=FALSE){
  # library(tibble)
  # library(dplyr)
  #expn data is log2, normalized counts
  #CDE has patients as rownames.
  #fillCol == character string of column name for fill colors
  #colorCol == character string of column name for border colors
  #colorCode is an option named vector of colors. If specifying color and fill manually, create a list with  length 2, with names c("fill","color")
  #single.col.outline is T/F for wether to susbet to dataframe for the border colors.
  #toHighlight is a character vector for if single.col.outline == TRUE. Gives the factor to highlight with borders from the colorColumn.
  #ellipse is T/F for an ellipse based on fillCol column
  expnData <- expnData[,intersect(rownames(CDE),colnames(expnData))]

  # print(dim(expnData))
  pca <- prcomp(t(expnData), scale=TRUE)
  summ <- summary(pca)

  scores <- as.data.frame(pca$x) %>%
    tibble::rownames_to_column("sample_id") %>%
    dplyr::inner_join(., dplyr::select(CDE,sample_id=matches("sample_id$"),
                                       everything()),
                      by="sample_id") %>%
    dplyr::select(sample_id, fillCol, colorCol, everything())

  #Plot function for  PC1 and either PC2 or anyother
  pca.plot_function <- function(scores,PC){

    idx <- as.numeric(gsub("[A-Za-z]{2}","", PC))


    pca.plot  <- ggplot2::ggplot(scores, ggplot2::aes_string(x="PC1", y=PC)) +
      labs(x=paste("PC1: ", round(summ$importance[2,1], digits=3)*100, "% variance"),
           y=paste(paste0(PC,": "), round(summ$importance[2,idx], digits=3)*100, "% variance")) +
      theme_numX  +
      theme(legend.text = element_text(size=14),
            legend.title = element_text(size=16))

    if(single.col.outline){
      pca.plot <- pca.plot +
        geom_point(size=5,stroke=0.2, alpha=1,shape=21,color="white",
                   ggplot2::aes_string(fill=fillCol)) +
        geom_point(data=subset(scores, scores[,colorCol] == toHighlight),
                   ggplot2::aes_string(x="PC1", y=PC, color=colorCol),
                   size=4, stroke=1, alpha=1,shape=21)

    }else{
      pca.plot <- pca.plot +
        geom_point(size=5, stroke=2, alpha=0.85,shape=21,
                   ggplot2::aes_string(fill=fillCol, color=colorCol))
    }

    if(!is.null(colorCode)){

        if(is.list(colorCode)){
          pca.plot <- pca.plot +
            scale_fill_manual(values=colorCode[["fill"]]) +
            scale_color_manual(values=colorCode[["color"]])
        }else{
          pca.plot <- pca.plot +
            scale_fill_manual(values=colorCode)
        }
    }

    if(ellipse){
        pca.plot <- pca.plot +
          stat_ellipse(data=scores, type="norm",
                       ggplot2::aes_string(x="PC1", y=PC, fill=fillCol),
                       geom="polygon", alpha=0.1)
    }

    return(pca.plot)
  }

  #PC1 and PC2 plot
  pc1.pc2 <- pca.plot_function(scores=scores,PC="PC2")

  #PC1 and PC3 plot
  if(PC3){
    pc1.pc3 <- pca.plot_function(scores=scores, PC="PC3")

  }


  res <- list("pca"=pca,"scores"=scores,"plot.1"=pc1.pc2)

  if(PC3){
    res[["plot.2"]] <- pc1.pc3
  }

  return(res)

}






