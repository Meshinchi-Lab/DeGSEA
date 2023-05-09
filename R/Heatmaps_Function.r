#Jenny Smith

#April 21, 2017

#purpose: to create a edgeR dge list given a vector of phenotypes, and a dataframe with counts.
#Also create colored dendrograms and heatmaps.

matchMatrix <- function(annodf, ExpnMatrix){
  #annodf is a cleaned CDE subset with patient IDs are rownames
  #expn matrix has patient IDs as colnames, genes as rownames

  #match the names in the expn matrix to names in a phenovector to ensure proper matching/order
  annoCol <- annodf[match(colnames(ExpnMatrix), rownames(annodf)),]

  # #find columns that are character class and update with "Not Available". keep numeric class as NA
  characterCols <- sapply(annoCol, class) == "character"
  idx <- is.na(annoCol[,characterCols])[,1]
  annoCol[idx,characterCols] <- "Unknown"

  #Set the rownames to columnames in expnmatrix to remove NAs from rownames
  rownames(annoCol) <- colnames(ExpnMatrix)

  return(annoCol)
}

colorVectors_asList <- function(df){
  library(RColorBrewer)
  #df of cleaned clinical characteristics. patient IDs as rownames

  cols <- colnames(df)

  list <- NULL
  for (col in cols){
    vector <- df[,col]
    groups <- unique(df[,col])
    colors <- c("turquoise3", "yellow", "blue", "firebrick1",
                "black","seagreen2", "maroon", "orchid", "cornflowerblue",
                "darkblue", "azure4", "chartreuse1", "darkmagenta","orange1",
                "deeppink", "darkslategray1", "green4", "navajowhite2",
                "brown3", "darkgoldenrod3", "deepskyblue1", "lightcoral",
                "mediumorchid", "saddlebrown")

    for (i in 1:length(groups)){
      vector[vector == groups[i]] <- colors[i]
    }
    list[[col]] <- vector
  }
  return(list)
}



#Create a color codes
colorCodes_aheatmap <- function(df,random=FALSE){
  #df of cleaned clinical characteristics, such as annoDF,   patient IDs as rownames
  library(RColorBrewer)
  cols <- colnames(df)

  list <- NULL
  for (col in cols){
    groups <- unique(df[,col])[order(unique(df[,col]))] #alphabetical order
    colors <- c("#E41A1C", "#377EB8" ,"#4DAF4A" ,
                 "#F781BF","blue1", "darkslategray3","burlywood3", "#984EA3", "#FF7F00",
                "seagreen2", "maroon", "orchid", "cornflowerblue", "yellow2",
                "darkblue", "azure4", "chartreuse1", "orange1",
                "deeppink", "darkslategray1", "green4", "navajowhite2",
                "brown3", "darkgoldenrod3", "deepskyblue1", "lightcoral",
                "mediumorchid", "darkmagenta")

    n <- length(groups)
    cc <- colors[1:n]
    names(cc) <- groups

    #update based on common group names used in my heatmaps
    if(any(grepl("Other$", groups))){
      cc["Other"] <- "lightgrey"}
    if(any(grepl("OtherAML", groups))){
      cc["OtherAML"] <- "lightgrey"}
    if(any(grepl("No$", groups))){
      cc["No"] <- "lightgrey"}
    if(any(grepl("^NBM$", groups))){
      cc["NBM"] <- "black"}

    #append color vector to list
    list[[col]] <- cc

  }

  return(list)
}



#' Create dendrogram from DGE object
#'
#' @param expnData normalized expression data or count data, patient IDs are column names and rows are genes.
#' @param pheno  a character vector with patient IDs as names, and the status for each in each group (eg pos,neg)
#' @param method cluster method
#' @param genelist a character vector with the genes of interest
#' @param add.count psuedocount value for log2 transformation
#' @param percent is the % of samples in the input expn matrix that must express a gene at 1 CPM. Filter to remove low count genes.
#' @param filterTopGenes boolean. If TRUE filter top 1000 most varied genes.
#' @param createDGE boolean. TRUE if raw counts are input.
#' @param log boolean. Yes/No for log2 transformation.
#'
#' @return list
#' @export
#'
#' @examples
#'
#' ex <- c('TBD')
dge_dendrograms <- function(expnData, pheno, method,
                            genelist=NULL,add.count=0.01, percent=0.01,
                            filterTopGenes=FALSE, createDGE=TRUE,log=FALSE){
  #df with count data, patient IDs are column names and rows are genes.
  #pheno is a character vector with patient IDs as names, and the status for each in each group (eg pos,neg)
  #genelist is a character vector with the genes of interest
  #percent is the % of samples in the input expn matrix that must express a gene at 1 CPM. Filter to remove low count genes.
  #set log=TRUE if providing a log2 expression dataframe.
  #filterTopGenes shuold be a logical. If TRUE filter top 1000 most varied genes.
  require(edgeR)
  suppressPackageStartupMessages(library(dendextend))
  library(matrixStats)

  expnData <- expnData[, intersect(names(pheno), colnames(expnData))] #ensure correct order, drop rows with nas just in case


  if(createDGE){
    #updated on 10/11/17 to add AML samples CPM cutoff before calc. norm factors.
    AML <- ! grepl("BM[0-9]|R[O0][0-9]", colnames(expnData))
    AMLsamples <- ncol(expnData[,AML])

    dge <- DGEList(counts = expnData)
    #updated on 7/26/18 - really only changes the heatmap coloring since the genes for each heatmap are already known to pass this threshold.
    keep.dge <- rowSums(cpm(dge)[,AML] >= 1) >= max(2,(percent*AMLsamples)) #X% of AML samples has cpm of at least 1 for a gene. This should help if looking for genes that are low or absent in NBMs vs AMLs

    dge <- dge[keep.dge,] #subset for those genes with cmp >= 1 per gene in AML samples
    dge <- calcNormFactors(dge)

    TMMCPM <- cpm(dge, normalized.lib.sizes = TRUE) #log = TRUE, prior.count = add.count

    #Use +1 to avoid outliers.... not necessarily true will revisit this. Using 0.01 actually cleans up the rare fusions clustering.
    TMMCPM <- as.data.frame(apply(TMMCPM, 2, function(x) log2(x + add.count)))


  #added else statement on 11/26/18
  }else{
    #yes, not accurate, bc it may not be TMMCPM norm. but for now, keep it.
    TMMCPM <- expnData

    #Use +0.5 to avoid outliers.... not necessarily true will revisit this.
    if(!log){
      TMMCPM <- as.data.frame(apply(TMMCPM, 2, function(x) log2(x + add.count))) #log2 transform counts
    }
  }

  if(is.null(genelist)){

    if(filterTopGenes){
      # calculate the variance for each gene
      rv <- rowVars(as.matrix(TMMCPM))

      # select the ntop genes by variance
      select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]

      #select the top 1000 most varied genes
      TMMCPM <- TMMCPM[select,]
      # print(dim(TMMCPM))

    }else{
      TMMCPM <- TMMCPM
    }

  }else{
    TMMCPM <- TMMCPM[which(rownames(TMMCPM) %in% genelist), ] #subset the matrix to genes of interest
  }


  d1 <- dist(t(TMMCPM), method = "euclidean", diag = FALSE,
             upper = FALSE) #sample distances WITHOUT SCALING
  d2 <- dist(TMMCPM, method = "euclidean", diag = FALSE,
             upper = TRUE) #gene distances WITHOUT SCaling

  samp.c1 <- hclust(d1, method = method, members = NULL) #sample clustering
  gene.c2 <- hclust(d2, method = method, members = NULL) #gene clustering


  list <- list(TMMCPM,samp.c1,gene.c2)
  names(list) <- c("TMMCPM","samp.c1", "gene.c2")

  return(list)
}


dge_dends_scale <- function(df, pheno, genelist, method){
  #df with count data, patient IDs are column names and rows are genes.
  #pheno is a character vector with patient IDs as names, and the status for each in each group (eg pos,neg)
  #genelist is a character vector with the genes of interest
  require(edgeR)
  suppressPackageStartupMessages(library(dendextend))

  df <- df[, intersect(names(pheno), colnames(df))] #ensure correct order, drop rows with nas just in case

  AML <- ! grepl("^BM|^RO", colnames(expData))
  AMLsamples <- ncol(expData[,AML])
  #updated 7/26/18
  keep.dge <- rowSums(cpm(dge)[,AML] >= 1) >= max(2, (0.01*AMLsamples)) #5% of AML samples has cpm of at least 1 for a gene

  dge <- DGEList(counts = df)
  dge <- dge[keep.dge,] #subset for those genes with cmp >= 1 per gene in AML samples
  dge <- calcNormFactors(dge)

  TMMCPM <- cpm(dge, normalized.lib.sizes = TRUE)
  TMMCPM <- TMMCPM[which(rownames(TMMCPM) %in% genelist), ] #subset the matrix to genes of interest

  names <- rownames(TMMCPM)
  #updated to +1 to avoid outliers
  TMMCPM <- as.data.frame(apply(TMMCPM, 2, function(x) log2(x + 1))) #log2 transform counts
  # TMMCPM10 <- as.data.frame(apply(TMMCPM, 2, function(x) log10(x + 0.01))) #log10 transform counts

  TMMCPM.n <- scale(t(TMMCPM)) #row-wise center scale (z-scores by genes)
  TMMCPM.tn <- t(TMMCPM.n) #put back in original orientation

  d1 <- dist(TMMCPM.n, method = "euclidean", diag = FALSE,
             upper = FALSE) #sample distances
  d2 <- dist(TMMCPM.tn, method = "euclidean", diag = FALSE,
             upper = TRUE) #gene distances

  c1 <- hclust(d1, method = method, members = NULL) #sample clustering
  c2 <- hclust(d2, method = method, members = NULL) #gene clustering


  list <- list(TMMCPM,names, d1,d2,c1,c2)
  names(list) <- c("TMMCPM","names", "d1", "d2", "c1", "c2")

  return(list)
}


colorDends <- function(hclustObject, colorCodes, group, textsize){
  #https://stackoverflow.com/questions/29265536/how-to-color-the-branches-and-tick-labels-in-the-heatmap-2
  #hclust object is the ward or other clustering matrix from hclust()
  #colorcodes is a vector with group status="color" format.
  #group is the character vector with groups status, one for each patient (patientIDs as names).  Usually called phenovector in my code.
  #textsize is a numeric vector of length 2 with size for 1) the labels of axes, and 2) patient IDs on leaves.

  suppressPackageStartupMessages(require(dendextend))


  N <- length(group)
  groups <- as.factor(group)

  dend <- dendextend::rotate(as.dendrogram(hclustObject),
                     order = c(1:N))

  branchCol <- colorCodes[as.numeric(groups)]
  branchCol <- branchCol[order.dendrogram(dend)]
  branchCol <- factor(branchCol, unique(branchCol))

  labels_colors(dend) <- colorCodes[group][order.dendrogram(dend)]
  dend <- color_branches(dend, clusters = as.numeric(branchCol), col = levels(branchCol))

  #example if want to add colored branches too
  # branch_col <- c("darkred", "forestgreen", "orange", "blue")[1:k]
  # dend <- as.dendrogram(dend)
  # dend <- color_branches(d, k=4, col=branch_col)


  par(mfrow=c(1,1), cex=textsize[1],
      mar=c(6.5, 7.5, 8.5, 2), pty="m",
      cex.main = 1, cex.lab = 0.85)
  plot(dend,
         xlab = " ", ylab=" ",  main=" ",
         axes=TRUE,
         type = "rectangle",
         cex.axis=textsize[2],
         horiz = FALSE)

}


#Function for Grouping and Splitting Trees
colorDends_Groups <- function(dendrogram, phenovector, k=NULL,h=NULL,
                              branchCol=NULL,colorcodes){
  library(RColorBrewer)
  suppressPackageStartupMessages(library(dendextend))
  #branchCol is a vector of length K, with the colors for the branches.
  #color codes is a named character vector for the known groups (eg c(pos="red, neg="blue)) for the leaves.

  if(class(dendrogram) == "list"){
    d <- dendrogram
    dend <- as.dendrogram(d$samp.c1)
  }else{
    dend <- dendrogram
  }

  #order_clusters_as_data allows you to order the labels from L to R bsed on the tree - NOT input dataframe for clustering.
  groups <- dendextend::cutree(dend, k=k, h=h, order_clusters_as_data = FALSE) #d$c1
  n <- unique(groups)

  #create vector of colors (as coded by 1=black, 2=red, etc)
  if(is.null(branchCol)){
    colors=n
    # colors=seq(1,k, by=1)
  }else{
    colors=branchCol
  }

  #select colors for the branches
  dend <- color_branches(dend, k=k, h=h, col=colors)
  labels_colors(dend) <- colorcodes[phenovector][order.dendrogram(dend)]
  labels_dend <- labels(dend)

  #Initialize variables for the sub-dendrograms
  dends <- list()
  group_labels <- list()

  for (i in 1:n[length(n)]){ #k
    group_labels[[i]] <- labels_dend[groups == i]
    labels_to_keep <- labels_dend[i != groups]
    dends[[i]] <- prune(dend, labels_to_keep)
  }

  res <- list(dend, groups, group_labels, dends)
  names(res) <- c("dend", "groups", "group_labels", "split_dends")

  return(res)
}



basicHeatmap <- function(ExpnMatrix, geneDend, sampleDend, colors,title){
  require(gplots)
  library(colorspace)
  suppressPackageStartupMessages(library(dendextend))

  #ExpnMatrix is the genes as rownames, patient IDs as colnames
  #genedend is from hclust object
  #sample dend is from hclust objest
  #rowlabels is the rownames of the initial expn matrix
  #colors is a character vector of colors of equal length of samples to illustrate the different groups


  # colors = c(seq(-3,-2,length=100),seq(-2,0.5,length=100),seq(0.5,6,length=100))
  # breaks = seq(-5,5,length.out = 300), #this MUST be checked between datasets or it will be inaccurate
  # colorPal <- colorRampPalette(c("blue", "white", "red"))(n=299)
  # colorPal <- colorRampPalette(c("darkgreen", "forestgreen", "green3", "green2", "black", "firebrick1", "red3", "red4", "darkred"))(n=299)
  colorPal <- colorRampPalette(c("deepskyblue4", "deepskyblue3", "deepskyblue2", "deepskyblue1","white","red1", "red2", "red3", "red4"))(n=299)
  N <- ncol(ExpnMatrix)
  N.genes <- nrow(ExpnMatrix)
  rowLabels <- rownames(ExpnMatrix)

  if ( N.genes < 50){
    cex=1.75
  }else if ( N.genes >= 50 & N.genes < 100){
    cex=1.25
  }else if(N.genes >= 100 & N.genes < 175){
    cex=0.75
  }else if (N.genes >= 175 & N.genes < 350){
    cex=0.6
  }else{
    cex=0.2
  }

  par(cex.main=1.5, cex=0.75, font=2, font.axis=1, lend=1)
  heatmap.2(as.matrix(ExpnMatrix),
              Colv=dendextend::rotate(as.dendrogram(sampleDend),
                                      order = c(N:1)),
              Rowv=as.dendrogram(geneDend),
              labRow=rowLabels,
              labCol = "",
              ColSideColors = colors,
              density.info="density", #density.info="density",
              trace = "none",
              scale="row",
              col = colorPal,
              cexRow=cex,
              margins=c(2,10),
              lwid=c(.8,3),
              lhei=c(.8,3),
              srtCol=75,
              adjCol=c(1,1),
              keysize=0.75,
              key.title="",
              key.ylab ="",
              key.par = list(cex=0.75),
              main=title)

}


annotationHeatmap <- function(ExpnMatrix, geneDend, sampleDend,annoDF, annoColors,main=NULL){
  require(NMF)
  suppressPackageStartupMessages(library(dendextend))
  library(colorspace)
  library(grid)

  # colorPal <- colorRampPalette(c("blue", "white", "red"))(n=299)
  #colorPal <- colorRampPalette(c("darkgreen", "forestgreen", "green3", "green2", "black", "firebrick1", "red3", "red4", "darkred"))(n=299)
  colorPal <- colorRampPalette(c("deepskyblue4", "deepskyblue3", "deepskyblue2", "deepskyblue1","white","red1", "red2", "red3", "red4"))(n=299)
  # colorPal <- colorRampPalette(c("chartreuse4", "chartreuse3", "chartreuse2", "chartreuse1","black", "magenta1", "magenta2", "magenta3", "magenta4"))(n=299)
  N <- ncol(ExpnMatrix)
  N.genes <- nrow(ExpnMatrix)
  rowLabels <- rownames(ExpnMatrix)

  if (N.genes <= 50 ){
    cex=1.5
  }else if ( N.genes >=  50 & N.genes <= 100){
    cex=0.8
  }else if(N.genes > 100 & N.genes < 175){
    cex=0.5
  }else{
    cex=0.2
  }

  par(cex.main=1.5, cex=0.75, font=2, font.axis=1, lend=1, mar=c(5,4,4,100)) #margins dont work...
  aheatmap(as.matrix(ExpnMatrix),
             Colv=dendextend::rotate(as.dendrogram(sampleDend),
                                     order = c(N:1)),
             Rowv=as.dendrogram(geneDend),
             annCol = annoDF,
             annColors = annoColors,
             scale="row",
             color = colorPal,
             cexRow=cex,
             cexCol=0.15,
             labCol = NA,
             breaks = 0,
             treeheight = 25,
             fontsize = 15,
             gp = gpar(cex=1.25,fontsize=20),
             main=main)
  # par(mar=c(5,4,4,100))

  return()
}



#' Create ComplexHeatmap annotation object
#'
#' @param expn the normalized expression values with genes as rownames
#' @param geneList a character vector
#' @param goi genes of interest to label on the Heatmap. Character vector of gene symbols
#' @param cc color codes as a named character vector.
#' @param CDE clinical data. has column called USI
#' @param cols character vector of column names in CDE
#' @param colorbar.height numeric heigh in cm
#'
#' @return list
#' @export
#'
#' @examples
#' ex <- c('TBD')
create_HA_Labs_Hmap <- function(expn,geneList, goi=NULL,cc=NULL, CDE, cols,
                                colorbar.height=5){
  #expn is the normalized expression values with genes as rownames
  #gene list a character vector
  #goi are genes of interest to highlight on the Heatmap. Character vector of gene symbols
  #cc are color codes in WHICH FORMAT?
  #CDE has column called USI
  #cols is a character vector of column names


  #Example on how to use the reuslts

  #hmap <- ComplexHmap(XXXX, hmap_anno_obj=res$annoColumn, XXX)
  # draw(hmap + res$geneLabels, heatmap_legend_side="right", annotation_legend_side="right")

    library(dplyr)
    suppressPackageStartupMessages(library(ComplexHeatmap))

    #subset the expression matix
    expn <- expn[geneList,]

    #select annation columns and create factor levels
    anno <- CDE %>%
      filter(USI %in% colnames(expn)) %>%
      dplyr::select(USI,all_of(cols)) %>%
      mutate(USI=factor(USI, levels = colnames(expn))) %>%
      arrange(USI)#ensure same order as the expn matrix


    #legend graphical parameters
    params <- list(show_legend=TRUE,
                  labels_gp= gpar(fontsize=12),
                  title_gp= gpar(fontsize=16),
                  nrow = 6, ncol=3,
                  by_row=TRUE)

    if(is.null(cc)){
      annoCol <- suppressWarnings(HeatmapAnnotation(df = anno[,-1],
                                 name="Main Groups",
                                 which="column",
                                 gap=unit(1,"mm"),
                                 border = T,
                                 annotation_height=unit(colorbar.height, "cm"),
                                 annotation_name_side="left",
                                 show_annotation_name = TRUE,
                                 annotation_legend_param = params)) #`annotation_height` is set with length of one while with multiple annotations, `annotation_height` is treated as `height`.")

    }else {
      #Dont know how to get a dummy vlaue for col(colors). NA and NULL produce an error.
      annoCol <- suppressWarnings(HeatmapAnnotation(df = dplyr::select(anno, all_of(cols)),
                                   name="Main Groups",
                                   col=cc,
                                   which="column",
                                   gap=unit(1,"mm"),
                                   border = T,
                                   show_annotation_name = TRUE,
                                   annotation_name_gp=gpar(fontsize=12),
                                   annotation_name_offset=unit(1,"mm"),
                                   annotation_name_side="left",
                                   annotation_height=unit(colorbar.height, "cm"),
                                   annotation_legend_param = params,
                                   simple_anno_size_adjust=TRUE))
    }

    res <- list("annoColumn"=annoCol)

    if(!is.null(goi)){
      regex <- paste0("^",goi,"$", collapse = "|")
      goi.idx <- grep(regex, rownames(expn))
      labels <- rownames(expn)[goi.idx] #needs to be numeric indexes for complexheatmap

      labs <- rowAnnotation(link = anno_mark(at=goi.idx,
                                             labels=labels,
                                             which="row",
                                             link_width=unit(1, "mm")),
                            width= unit(1, "mm") + max_text_width(labels),
                            gp=gpar(fontsize=4))

      res[["geneLabels"]] <-  labs
    }


    return(res)
}



#' Create ComplexHeatmap plot
#'
#' @param mat the normalized, log2 (usually) transformed counts
#' @param name is the title
#' @param scale boolean whether to scale by row
#' @param threshold whether to make all z-scores in a certain range.
#' @param hmap_anno_obj from HeatmapAnnotation() function
#' @param hmap_anno_obj_genes from HeatmapAnnotation() function
#' @param space.type value for color palette, like sRGB or LAB
#' @param color_palette colors vector from circlize::colorRamp2
#' @param split data.frame with the groups identifying which rows to split.
#' @param cluster.method method for dendrogram.
#' @param dge_dendrograms.res the object from DeGSEA::dge_dendrograms()
#' @param samp_dend_order  numeric vector or character vector of column names from the matrix (mat) or the dge_dengrograms.res$TMMCPM matix, in the desired order.
#'
#' @return complexHeatmap object
#' @export
#'
#' @examples
#' ex <- c('TBD')
ComplexHmap <- function(mat, name="z-scores",
                        scale=TRUE,threshold=FALSE,
                        hmap_anno_obj,
                        hmap_anno_obj_genes=NULL,
                        space.type="sRGB",
                        color_palette=NULL,
                        split=NULL, cluster.method="ward.D2",
                        dge_dendrograms.res=NULL,
                        samp_dend_order=NULL){
  #mat is the normalized, log2 (usually) transformed counts
  #name is the title
  #scale is whether to scale by row
  #color palette is a colorRamp2() object.
  #threshold is whether to make all z-scores in a certain range.
  #hmap_anno_obj is from HeatmapAnnotation() function
  #space.type is for the color/shades on the heatmap. See ?Heatmap for all the choices.
  #dge_dendrograms.res is the list object output from the dge_dendrograms() function.
  #samp_dend_order is the numeric vector or character vector of column names from the matrix (mat) or the dge_dengrograms.res$TMMCPM matix, in the desired order.

  suppressPackageStartupMessages(library(ComplexHeatmap))
  suppressPackageStartupMessages(require(circlize))
  library(RColorBrewer)
  suppressPackageStartupMessages(library(dendextend))
  ht_opt$message = FALSE

  if(threshold){
    mat[mat > 4] <- 4
    mat[mat < -4] <- -4
  }

  if(is.null(color_palette)){
    pal <- colorRamp2(c(-4,-2, 0, 2, 4),
                      c("deepskyblue3", "deepskyblue","white", "red", "red3"),
                      space=space.type)
  }else{
    pal <- color_palette
  }
  # col <- colorRampPalette(c("cyan1", "cyan2", "cyan3", "cyan4","azure4","magenta4", "magenta3", "magenta2", "magenta1"))(n=299)
  #colorRamp2(c(-2, 0, 4), c("deepskyblue","white", "red"), space="RGB") #use for breaks.
  # colorPal <- colorRampPalette(c("deepskyblue4", "deepskyblue3", "deepskyblue2", "deepskyblue1","white","red1", "red2", "red3", "red4"))(n=299)

  #legend graphical parameters
  params <-  list(color_bar="continuous",
       legend_direction="horizontal",
       title_position="leftcenter",
       legend_width=unit(5,"cm"),
       legend_height=unit(5,"cm"),
       title_gp=gpar(fontsize=10,
                     fontface="bold"))

  if(is.null(dge_dendrograms.res)){

    if(scale){
      mat <- t(scale(t(mat))) ##rowwise scaling
      print(range(mat))
    }

    hmap <- Heatmap(mat,
                    name=name,
                    col=pal, #use for breaks.,

                    heatmap_legend_param=params,
                    row_title="Genes",
                    row_title_side="left",
                    row_title_gp=gpar(fontsize=15, fontface="bold"),
                    show_row_names=FALSE,
                    row_names_gp=gpar(fontsize=3),

                    column_title="",
                    column_title_side="top",
                    column_title_gp=gpar(fontsize=15, fontface="bold"),
                    column_title_rot=0,
                    show_column_names=FALSE,

                    row_dend_width=unit(8,"mm"),
                    column_dend_height=unit(22.5,"mm"),

                    top_annotation=hmap_anno_obj,
                    right_annotation = hmap_anno_obj_genes,
                    split=split,

                    #NOTE: This will cluster samples based on the z-score transformed matrix if scale=TRUE here
                    clustering_distance_columns="euclidean",
                    clustering_method_columns=cluster.method,
                    column_dend_reorder=ncol(mat):1,
                    clustering_distance_rows="euclidean",
                    clustering_method_rows=cluster.method)

                    #NOTE: use code below for correlation as distance.
                    # clustering_distance_columns=function(x) as.dist(1-cor(t(x))))

  # Use this if defining your own clustering
  }else if(!is.null(dge_dendrograms.res)){

    if(scale){
      mat <- t(scale(t(dge_dendrograms.res[[1]]))) #rowwise scaling
      print(range(mat))
    }else{
      mat <- dge_dendrograms.res[[1]]
    }

    # min <- min(mat)+2
    # max <- max(mat)-2
    # print(c(min,max))

    if(is.null(samp_dend_order)){
      clust <- dendextend::rotate(as.dendrogram(dge_dendrograms.res$samp.c1),
                         order=c(ncol(mat):1))
    }else{
      clust <- dendextend::rotate(as.dendrogram(dge_dendrograms.res$samp.c1),
                                order=samp_dend_order)

    }

    hmap <- Heatmap(mat,
                    name=name,
                    col=pal,

                    heatmap_legend_param=params,
                    row_title="Genes",
                    row_title_side="left",
                    row_title_gp=gpar(fontsize=15,
                                      fontface="bold"),
                    show_row_names=FALSE,
                    row_names_gp=gpar(fontsize=3),

                    column_title="",
                    column_title_side="top",
                    column_title_gp=gpar(fontsize=15,
                                         fontface="bold"),
                    column_title_rot=0,
                    show_column_names=FALSE,

                    row_dend_width=unit(8,"mm"),
                    column_dend_height=unit(22.5,"mm"),

                    top_annotation=hmap_anno_obj,
                    right_annotation = hmap_anno_obj_genes,
                    split=split,

                    #cluster the rows (genes) using the scaled matrix rather than pre-defined order in the dge.dendrograms object.
                    #the clustering on non-scaled then visualizing the scaled rows really changes the expression pattern on the rows.
                    clustering_distance_rows="euclidean",
                    clustering_method_rows=cluster.method,
                    # cluster_rows = as.dendrogram(dge_dendrograms.res$gene.c2),
                    # row_dend_reorder = FALSE,
                    cluster_columns = clust,
                    column_dend_reorder=FALSE)
  }

  return(hmap)
}



#Pheatmap
quickPheatmap <- function(expn,geneList, clinData,
                          cols=c("Age.Category", #default columns
                                 "Cytogenetic.Category.1", "SNVs",
                                 "Cytogenetic.Category.2","Rare.Fusions"),
                          Add.Anno.Cols=NULL,
                          annots_col_colors=NULL){
  #expn has rows with genes as rownames
  #clinData has rowns with patient USIs as rownames - matches the expn data column names

  library(pheatmap)

  len <- 299
  col <- colorRampPalette(c("deepskyblue4", "deepskyblue3", "deepskyblue2", "deepskyblue1","white","red1", "red2", "red3", "red4"))(n=len)

  e <- expn[geneList, intersect(colnames(expn), rownames(clinData))]
  clinData <- clinData[intersect(colnames(expn), rownames(clinData)), ]

  #if scaling, this doesnt work. would need to scale the expression matrix myself
  # Breaks <- c(seq(min(e), 0, length.out=ceiling(len/2) + 1),
  #             seq(max(e)/len, max(e), length.out=floor(len/2)))

  cols.colorbar <- cols
  cols.colorbar <- c(Add.Anno.Cols, cols.colorbar)

  anno.col <- clinData %>%
    dplyr::select(all_of(cols.colorbar))

  if(is.null(annots_col_colors)){
    annots_col_colors <- lapply(anno.col, function(column){
        grps <- unique(column)
        n <- length(grps)

        if(n < 3 | n > 12){
          pal <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n)
        }else{
          pal <- RColorBrewer::brewer.pal(n, "Paired")
        }

        names(pal) <- grps
        return(pal)
      })
  }


  p <- pheatmap::pheatmap(mat= e ,
                     color=col,
                     border_color="black",
                     scale = "row",
                     # breaks = Breaks,
                     cluster_rows=FALSE,
                     cluster_cols=TRUE,
                     clustering_method="complete",
                     clustering_distance_cols="euclidean",
                     annotation_names_col=TRUE,
                     annotation_col = anno.col,
                     annotation_colors = annots_col_colors,
                     # gaps_row=11,
                     # cutree_cols = 2,
                     # labels_row = labs,
                     show_rownames=TRUE,
                     fontsize = 5,
                     show_colnames=FALSE)


  return(p)

}
