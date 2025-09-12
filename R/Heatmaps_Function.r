#Jenny Smith

#April 21, 2017

#purpose: to create a edgeR dge list given a vector of phenotypes, and a dataframe with counts.
#Also create colored dendrograms and heatmaps.


#' Create a list of color codes for all columns in a datafram
#'
#' @param df df of cleaned clinical characteristics. patient IDs as rownames
#'
#' @returns list
#' @export
#'
#' @examples
#' ex <- c("TBD")
colorVectors_asList <- function(df){
  # library(RColorBrewer)
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



#
#' Create a color codes
#'
#' @param df cleaned clinical characteristics, with patient IDs as rownames
#' @param random boolean. should colors be selected by random sampling.
#'
#' @return list
#' @export
#'
#' @examples
#' ex <- c('TBD')
colorCodes_aheatmap <- function(df,random=FALSE,AML_groups=FALSE){
  #df of cleaned clinical characteristics, such as annoDF,   patient IDs as rownames
  # library(RColorBrewer)
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

    if (AML_groups){
      #update based on common group names used in my heatmaps
      if(any(grepl("Other$", groups))){
        cc["Other"] <- "lightgrey"}

      if(any(grepl("OtherAML", groups))){
        cc["OtherAML"] <- "lightgrey"}

      if(any(grepl("No$", groups))){
        cc["No"] <- "lightgrey"}

      if(any(grepl("^NBM$", groups))){
        cc["NBM"] <- "black"}
    }

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
                            filterTopGenes=FALSE, createDGE=TRUE,log=FALSE,
                            exclude_samples = NULL, scale_data = FALSE){
  #df with count data, patient IDs are column names and rows are genes.
  #pheno is a character vector with patient IDs as names, and the status for each in each group (eg pos,neg)
  #genelist is a character vector with the genes of interest
  #percent is the % of samples in the input expn matrix that must express a gene at 1 CPM. Filter to remove low count genes.
  #set log=TRUE if providing a log2 expression dataframe.
  #filterTopGenes shuold be a logical. If TRUE filter top 1000 most varied genes.
  # require(edgeR)
  # suppressPackageStartupMessages(library(dendextend))
  # library(matrixStats)

  expnData <- expnData[, intersect(names(pheno), colnames(expnData))] #ensure correct order, drop rows with nas just in case

  if( ! is.null(exclude_samples)){
    #updated on 10/11/17 to add AML samples CPM cutoff before calc. norm factors.
    selected_samps <- ! grepl(exclude_samples, colnames(expnData))
  } else {
    selected_samps <- colnames(expnData)
  }

  if(createDGE){
    dge <- edgeR::DGEList(counts = expnData)
    #updated on 7/26/18 - really only changes the heatmap coloring since the genes for each heatmap are already known to pass this threshold.
    N <- ncol(expnData[,selected_samps])
    keep.dge <- rowSums(edgeR::cpm(dge)[,selected_samps] >= 1) >= max(2,(percent*N)) #X% of samples has cpm of at least 1 for a gene. This should help if looking for genes that are low or absent in NBMs vs AMLs

    dge <- dge[keep.dge,] #subset for those genes with cmp >= 1 per gene in samples
    dge <- edgeR::calcNormFactors(dge)

    TMMCPM <- edgeR::cpm(dge, normalized.lib.sizes = TRUE) #log = TRUE, prior.count = add.count

    #Use +1 to avoid outliers.... not necessarily true will revisit this. sample_idng 0.01 actually cleans up the rare fsample_idons clustering.
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

  if( scale_data ){
    TMMCPM <-  t(scale(t(TMMCPM))) #row-wise center scale (z-scores by genes)
  }

  d1 <- dist(t(TMMCPM), method = "euclidean", diag = FALSE,
             upper = FALSE) #sample distances WITHOUT SCALING
  d2 <- dist(TMMCPM, method = "euclidean", diag = FALSE,
             upper = TRUE) #gene distances WITHOUT SCaling

  samp.c1 <- stats::hclust(d1, method = method, members = NULL) #sample clustering
  gene.c2 <- stats::hclust(d2, method = method, members = NULL) #gene clustering


  list <- list(TMMCPM,samp.c1,gene.c2)
  names(list) <- c("TMMCPM","samp.c1", "gene.c2")

  return(list)
}

#' Generate a dendrogram with branches colored by groups
#'
#' @param hclustObject hclust object is the ward or other clustering matrix from hclust()
#' @param colorCodes colorcodes is a vector with group status="color" format
#' @param group group is the character vector with groups status, one for each patient (sample_id as names)
#' @param textsize textsize is a numeric vector of length 2 with size for 1) the labels of axes, and 2) sample_id on leaves
#'
#' @returns plot
#' @export
#'
#' @examples
#' ex <- c("TBD")
colorDends <- function(hclustObject, colorCodes, group, textsize){
  #https://stackoverflow.com/questions/29265536/how-to-color-the-branches-and-tick-labels-in-the-heatmap-2
  #hclust object is the ward or other clustering matrix from hclust()
  #colorcodes is a vector with group status="color" format.
  #group is the character vector with groups status, one for each patient (patientIDs as names).  Usually called phenovector in my code.
  #textsize is a numeric vector of length 2 with size for 1) the labels of axes, and 2) patient IDs on leaves.


  #example if want to add colored branches too
  # branch_col <- c("darkred", "forestgreen", "orange", "blue")[1:k]
  # dend <- as.dendrogram(dend)
  # dend <- color_branches(d, k=4, col=branch_col)


  N <- length(group)
  groups <- as.factor(group)

  dend <- dendextend::rotate(dendextend::as.dendrogram(hclustObject),
                     order = c(1:N))

  branchCol <- colorCodes[as.numeric(groups)]
  branchCol <- branchCol[dendextend::order.dendrogram(dend)]
  branchCol <- factor(branchCol, unique(branchCol))

  labels_colors(dend) <- colorCodes[group][dendextend::order.dendrogram(dend)]
  dend <- dendextend::color_branches(dend,
                         clusters = as.numeric(branchCol),
                         col = levels(branchCol))

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

#' Function for Grouping and Splitting Trees
#'
#' @param dendrogram a dendrogram object
#' @param phenovector character vector with groups status, one for each patient (sample_id as names)
#' @param k number of clusters
#' @param h height to cut the dendrogram
#' @param branchCol branchCol is a vector of length K, with the colors for the branches.
#' @param colorcodes color codes is a named character vector for the known groups (eg c(pos="red, neg="blue)) for the leaves.
#'
#' @returns list
#' @export
#'
#' @examples
#' ex <- c("TBD")
colorDends_Groups <- function(dendrogram, phenovector, k=NULL,h=NULL,
                              branchCol=NULL,colorcodes){

  # branchCol is a vector of length K, with the colors for the branches.
  # color codes is a named character vector for the known groups (eg c(pos="red, neg="blue)) for the leaves.

  if(inherits(dendrogram, "list")){
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
  dend <- dendextend::color_branches(dend, k=k, h=h, col=colors)
  labels_colors(dend) <- colorcodes[phenovector][dendextend::order.dendrogram(dend)]
  labels_dend <- labels(dend)

  #Initialize variables for the sub-dendrograms
  dends <- list()
  group_labels <- list()

  for (i in 1:n[length(n)]){ #k
    group_labels[[i]] <- labels_dend[groups == i]
    labels_to_keep <- labels_dend[i != groups]
    dends[[i]] <- dendextend::prune(dend, labels_to_keep)
  }

  res <- list(dend, groups, group_labels, dends)
  names(res) <- c("dend", "groups", "group_labels", "split_dends")

  return(res)
}


#' Create ComplexHeatmap annotation object
#'
#' @param expn the normalized expression values with genes as rownames
#' @param geneList a character vector
#' @param goi genes of interest to label on the Heatmap. Character vector of gene symbols
#' @param cc color codes as a named character vector.
#' @param CDE clinical data. has column called sample_id
#' @param cols character vector of column names in CDE
#' @param colorbar.height numeric heigh in cm
#'
#' @return list
#' @export
#'
#' @examples
#' ex <- c('TBD')
create_HA_Labs_Hmap <- function(expn,geneList,
                                CDE, cols,
                                goi=NULL,cc=NULL,
                                colorbar.height=5){
  #expn is the normalized expression values with genes as rownames
  #gene list a character vector
  #goi are genes of interest to highlight on the Heatmap. Character vector of gene symbols
  #cc are color codes in WHICH FORMAT?
  #CDE has column called sample_id
  #cols is a character vector of column names


    #Example on how to use the reuslts:
    #hmap <- ComplexHmap(XXXX, hmap_anno_obj=res$annoColumn, XXX)
    # draw(hmap + res$geneLabels, heatmap_legend_side="right", annotation_legend_side="right")


    #subset the expression matix
    expn <- expn[geneList,]

    #select annation columns and create factor levels
    anno <- CDE %>%
      dplyr::filter(sample_id %in% colnames(expn)) %>%
      dplyr::select(sample_id, dplyr::all_of(cols)) %>%
      dplyr::mutate(sample_id=factor(sample_id, levels = colnames(expn))) %>%
      dplyr::arrange(sample_id)#ensure same order as the expn matrix


    #legend graphical parameters
    params <- list(show_legend=TRUE,
                  labels_gp = grid::gpar(fontsize=12),
                  title_gp = grid::gpar(fontsize=16),
                  nrow = 6, ncol=3,
                  by_row=TRUE)

    if(is.null(cc)){
      annoCol <- suppressWarnings(
          ComplexHeatmap::HeatmapAnnotation(df = anno[,-1],
                                   name="Main Groups",
                                   which="column",
                                   gap=unit(1,"mm"),
                                   border = T,
                                   annotation_height=unit(colorbar.height, "cm"),
                                   annotation_name_side="left",
                                   show_annotation_name = TRUE,
                                   annotation_legend_param = params)
        ) #`annotation_height` is set with length of one while with multiple annotations, `annotation_height` is treated as `height`.")

    }else {
      #Dont know how to get a dummy vlaue for col(colors). NA and NULL produce an error.
      annoCol <- suppressWarnings(
        ComplexHeatmap::HeatmapAnnotation(df = dplyr::select(anno, all_of(cols)),
                                   name="Main Groups",
                                   col=cc,
                                   which="column",
                                   gap=unit(1,"mm"),
                                   border = T,
                                   show_annotation_name = TRUE,
                                   annotation_name_gp=grid::gpar(fontsize=12),
                                   annotation_name_offset=unit(1,"mm"),
                                   annotation_name_side="left",
                                   annotation_height=unit(colorbar.height, "cm"),
                                   annotation_legend_param = params,
                                   simple_anno_size_adjust=TRUE)
        )
    }

    res <- list("annoColumn"=annoCol)

    if(!is.null(goi)){
      regex <- paste0("^",goi,"$", collapse = "|")
      goi.idx <- grep(regex, rownames(expn))
      labels <- rownames(expn)[goi.idx] #needs to be numeric indexes for complexheatmap

      labs <- ComplexHeatmap::rowAnnotation(link = ComplexHeatmap::anno_mark(at=goi.idx,
                                             labels=labels,
                                             which="row",
                                             link_width=unit(1, "mm")),
                              width = unit(1, "mm") +
                              max_text_width(labels),
                              gp = grid::gpar(fontsize=4))

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

  if(threshold){
    mat[mat > 4] <- 4
    mat[mat < -4] <- -4
  }

  if(is.null(color_palette)){
    pal <- circlize::colorRamp2(c(-4,-2, 0, 2, 4),
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
       title_gp=grid::gpar(fontsize=10,
                     fontface="bold"))

  if(is.null(dge_dendrograms.res)){

    if(scale){
      mat <- t(scale(t(mat))) ##rowwise scaling
      print(range(mat))
    }

    hmap <- ComplexHeatmap::Heatmap(mat,
                    name=name,
                    col=pal, #use for breaks.,

                    heatmap_legend_param=params,
                    row_title="Genes",
                    row_title_side="left",
                    row_title_gp=grid::gpar(fontsize=15, fontface="bold"),
                    show_row_names=FALSE,
                    row_names_gp=grid::gpar(fontsize=3),

                    column_title="",
                    column_title_side="top",
                    column_title_gp=grid::gpar(fontsize=15, fontface="bold"),
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

    if(is.null(samp_dend_order)){
      clust <- dendextend::rotate(as.dendrogram(dge_dendrograms.res$samp.c1),
                         order=c(ncol(mat):1))
    }else{
      clust <- dendextend::rotate(as.dendrogram(dge_dendrograms.res$samp.c1),
                                order=samp_dend_order)

    }

    hmap <- ComplexHeatmap::Heatmap(mat,
                    name=name,
                    col=pal,

                    heatmap_legend_param=params,
                    row_title="Genes",
                    row_title_side="left",
                    row_title_gp=grid::gpar(fontsize=15,
                                      fontface="bold"),
                    show_row_names=FALSE,
                    row_names_gp=grid::gpar(fontsize=3),

                    column_title="",
                    column_title_side="top",
                    column_title_gp=grid::gpar(fontsize=15,
                                         fontface="bold"),
                    column_title_rot=0,
                    show_column_names=FALSE,

                    row_dend_width=unit(8,"mm"),
                    column_dend_height=unit(22.5,"mm"),

                    top_annotation = hmap_anno_obj,
                    right_annotation = hmap_anno_obj_genes,
                    split=split,

                    #cluster the rows (genes) sample_idng the scaled matrix rather than pre-defined order in the dge.dendrograms object.
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



# Pheatmap
quickPheatmap <- function(expn,geneList, clinData,
                          cols=NULL,
                          Add.Anno.Cols=NULL,
                          annots_col_colors=NULL){
  #expn has rows with genes as rownames
  #clinData has rowns with patient sample_ids as rownames - matches the expn data column names

  # library(pheatmap)

  len <- 299
  col <- grDevices::colorRampPalette(c("deepskyblue4", "deepskyblue3", "deepskyblue2", "deepskyblue1","white","red1", "red2", "red3", "red4"))(n=len)

  e <- expn[geneList, intersect(colnames(expn), rownames(clinData))]
  clinData <- clinData[intersect(colnames(expn), rownames(clinData)), ]

  #if scaling, this doesnt work. would need to scale the expression matrix myself
  # Breaks <- c(seq(min(e), 0, length.out=ceiling(len/2) + 1),
  #             seq(max(e)/len, max(e), length.out=floor(len/2)))

  if(is.null(cols)){
    #default columns
    cols <- c("Age.Category",
           "Cytogenetic.Category.1", "SNVs",
           "Cytogenetic.Category.2","Rare.Fsample_idons")
  }else{
    cols.colorbar <- cols
  }

  cols.colorbar <- c(Add.Anno.Cols, cols.colorbar)

  anno.col <- clinData %>%
    dplyr::select(dplyr::all_of(cols.colorbar))

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
