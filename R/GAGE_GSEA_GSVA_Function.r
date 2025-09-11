#Jenny Smith

#June 21, 2017

#Purpose: Perform gene set enrichement analysis with A) Inpur normalized counts (from limma)
# and B) input log2 fold changes for significantly differentially expressed genes.




#
#' DE analysis with the gsva matrix.
#'
#' @param gsva_matrix output scores from gsva
#' @param clinData dataframe of clinical data with patients as rownames
#' @param col covariate column name
#' @param p.value p-value threshold
#'
#' @returns list
#' @export
#'
#' @examples
#' ex <- c("TBD")
gsva_DE <- function(gsva_matrix,clinData,col,p.value=0.001){
  #clinData must be factor leveled!

  #Example:
  # X <-  sample_info[samps_all, ] %>%
  #   mutate(USI=Sample) %>%
  #   mutate_at(vars(NUP98.Rearranged.Groups),
  #             ~factor(., levels=c("OtherAML","NUP98-X"))) %>%
  #   set_rownames(.$Sample)
  # gsva_DE(gsva_matrix = gsva.res.all[,X$Sample],
  #         clinData=X,
  #         col="NUP98.Rearranged.Groups",
  #         p.value = 0.25)


  # library(edgeR)
  # library(limma)

  gsva <- gsva_matrix[,clinData$Sample]
  design <- model.matrix(formula(paste("~0 +", col)),
                         data=clinData)
  colnames(design) <- c("Ref","Comparitor")



  cm <- makeContrasts("Comparitor-Ref", levels = design)


  fit <- lmFit(gsva, design = design)
  fit2 <- contrasts.fit(fit, contrasts = cm)
  fit2 <- eBayes(fit2)

  allGeneSets <- topTable(fit2, p.value = p.value,number=Inf) %>%
    rownames_to_column("GeneSet") %>%
    dplyr::arrange(desc(logFC), desc(adj.P.Val))

  res <- list("contrast"=cm,"fit"=fit2, "gene_sets"=allGeneSets)

  return(res)
}


#' GAGE gene-set enrichment analysis
#'
#' @param df df has genes symbols as rownames and patients as columns with normalized expression data
#' @param type type is either limma differential expression results (FC) or expression matrix of normalized counts (expn)
#' @param geneset optional list gene-sets
#' @param ref column numbers for reference samples
#' @param samp column numbers for samples in condition of interest
#' @param min_size min size of gene-sets
#'
#' @returns list
#' @export
#'
#' @examples
#' ex <- c("TBD")
GAGE_GSA <- function(df, type, geneset=NULL, ref=NULL,samp=NULL, min_size=15){
  #df has genes symbols as rownames and patients as columns with normalized expression data
  #type is either limma differential expression results (FC) or expression matrix of normalized counts (expn)


  if(is.null(geneset)){
    message(
      'Load the kegg human datasets with data("kegg.gs") and data("go.sets.hs") and data("egSymb")'
      )

    #create objects to hold the information from the datasets loaded above
    kg.hsa=kegg.gsets("hsa")

    #Examine all Kegg genesets, signalling, metabolic, and disease.
    kegg.gs <- kg.hsa$kg.sets

    #convert to gene symbols
    geneset <- lapply(kegg.gs, eg2sym)

  }else{
    geneset = geneset
  }

  if (type == "FC"){
    FCs <- df$logFC #Note: the FCs should include ALL genes tested - so approx. 17-20,000 regarless of p-value.
    names(FCs) <- rownames(df)

    gsa <- gage(FCs, gsets=geneset, same.dir = TRUE,
                saaTest = gs.tTest, #gs.KSTest
                ref=ref, samp=samp, set.size = c(20, 450)) #examine all differntially expressed genes.

    sigRes <- subset(gsa$greater, gsa$greater[,"q.val"] <  0.1 &
                       !is.na(gsa$greater[,"q.val"]))

    #Final results list
    res <- list(gsa, sigRes)
    names(res) <- c("gsa", "SigPaths")

    if (nrow(sigRes) > 1){
      essSets.Up <- esset.grp(gsa$greater, FCs, gsets=geneset,
                           use.q = TRUE,cutoff = 0.1,
                           ref=NULL, sampl=NULL,
                           output = FALSE,
                           test4up = FALSE,
                           same.dir = FALSE)

      res[["essSets"]] <-  essSets.Up
    }

  }else if(type=="expn"){
    gsa <- gage(df, gsets=geneset,
                same.dir = TRUE, ref= ref,
                samp=samp,compare ="unpaired",
                set.size = c(20, 450))

    sigRes.Up <- subset(gsa$greater, gsa$greater[,"q.val"] <  0.1 &
                       !is.na(gsa$greater[,"q.val"]))

    sigRes.Dn <- subset(gsa$less, gsa$less[,"q.val"] <  0.1 &
                          !is.na(gsa$less[,"q.val"]))

    #Final results list
    res <- list(gsa, sigRes.Up, sigRes.Dn)
    names(res) <- c("gsa", "SigPaths.Up","SigPaths.Dn")

    if (nrow(sigRes.Up) > 1){
      essSets.Up <- esset.grp(gsa$greater, df, gsets=geneset,
                              use.q = TRUE,
                              ref=ref, sampl=samp,
                              compare = "unpaired",
                              output = FALSE,
                              cutoff = 0.1,
                              test4up = TRUE)

      res[["essSets.Up"]] <-  essSets.Up

      if (nrow(sigRes.Dn) > 1){
      essSets.Dn <- esset.grp(gsa$less, df, gsets=geneset,
                              use.q = TRUE,
                              ref=ref, sampl=samp,
                              compare = "unpaired",
                              output = FALSE,
                              cutoff = 0.1,
                              test4up = FALSE)

      res[["essSets.Dn"]] <-  essSets.Dn

      }
    }
  }

  return(res)
}

#


#' Wrapper Function to run GAGE GSEA on the output for the DEGs pipeline
#'
#' @param twoGroups_DEGs.res output from twoGroups_DEGs()
#' @param type type is either limma differential expression results (FC) or expression matrix of normalized counts (expn)
#' @param geneset optional list gene-sets
#' @param method voom or trend
#' @param min_size min size of gene-sets
#' @param ... others
#'
#' @returns list
#' @export
#'
#' @examples
#' ex <- c("TBD")
gage_from_pipeline <- function(twoGroups_DEGs.res,type,geneset=NULL, method="voom",min_size=15,... ){

  if (type == "FC"){
    #change this to use ALL genes tested (~20K genes regardless of significance)
    df <- twoGroups_DEGs.res$DE$DE
    ref <- NULL
    cond <- NULL

  }else if (type == "expn" ) {

    if(method=="voom"){
      df <- twoGroups_DEGs.res$DE$Voom$E
    }else if(method=="trend"){
      df <- twoGroups_DEGs.res$DE$dge
    }

    ref <- which(colnames(df) %in%
                   names(twoGroups_DEGs.res$phenovector[twoGroups_DEGs.res$phenovector == "GroupB"]))

    cond <- which(colnames(df) %in%
                    names(twoGroups_DEGs.res$phenovector[twoGroups_DEGs.res$phenovector == "GroupA"]))
  }

  #Run GSA
  print(min_size)
  GSA <- GAGE_GSA(df=df, type=type, geneset = geneset, ref=ref, samp=cond, min_size=min_size)

  return(GSA)

}


##Simple function to read in a .gmt file and return a list of pathways
#https://github.com/arcolombo/junk/blob/5d39a5893bb7dd6328d587791b9735de14b2ff45/R/qusageArm.R
#Directly from github. Cannot install the library.
#' Read in a .gmt file and return a list of pathways
#'
#' @param file file path to gmt
#'
#' @returns list
#' @export
#'
#' @examples
#' ex <- c("TBD")
read.gmt = function(file){
  if(!grepl("\\.gmt$",file)[1]){stop("Pathway information must be a .gmt file")}
  geneSetDB = readLines(file)                                ##read in the gmt file as a vector of lines
  geneSetDB = strsplit(geneSetDB,"\t")                       ##convert from vector of strings to a list
  names(geneSetDB) = sapply(geneSetDB,"[",1)                 ##move the names column as the names of the list
  geneSetDB = lapply(geneSetDB, "[",-1:-2)                   ##remove name and description columns
  geneSetDB = lapply(geneSetDB, function(x){x[which(x!="")]})##remove empty strings
  return(geneSetDB)
}

#
#' Simply function to write custom list of gene-sets to gmt file format.
#'
#' @param list_of_genesets a list of gene-sets
#' @param Description_list a paired list of gene-set descriptions
#' @param filename output filenames
#'
#' @returns nothing
#' @export
#'
#'
#' @examples
#' ex <- c("TBD")
write.gmt <- function(list_of_genesets,Description_list=NA, filename="gene_sets.gmt"){

  if(is.na(Description_list)){
    Description_list <- as.list(rep(NA,length(list_of_genesets)))
  }

  for (i in 1:length(list_of_genesets)){
    append <- ifelse(i==1,FALSE,TRUE)
    list_of_genesets[[i]] <- c(Description_list[[i]],list_of_genesets[[i]])

    line <- paste0(c(names(list_of_genesets)[i], list_of_genesets[[i]]),
                   collapse = "\t")
    line <- paste0(line,"\n")

    # cat(line,
    #     file = filename,
    #     sep = " ",
    #     fill = FALSE,
    #     labels = NULL,
    #     append = append)
  }

  print(paste("Finshed writing gene-sets to file", filename))
}


#Needs to be tested
extract_core_genes <- function(gage_from_pipeline.res,direction="up", ret.genesets=FALSE){

  gage.res.sig <- gage_from_pipeline.res[sapply(gage_from_pipeline.res,length) == 5]

  if(direction=="up"){
    idx <- "essSets.Up"
    idx.num <- 4
  }else if (direction=="down"){
    idx <- "essSets.Dn"
    idx.num <- 5
  }

  #select the correct index for the pathway core genesets
  core.genes <- lapply(gage.res.sig, `[[`, idx.num)
  core.genes <- lapply(core.genes, `[[`, 7)

  if(ret.genesets){
    return(core.genes)
  }

  #Define an empty list to append to
  paths.core.genes <- list()
  #index position for the new list
  n <- 1

  #for loop to un-nest the core.genes, so that they are not in a depth == 2 list
  for (i in 1:length(core.genes)){

    gs <- core.genes[[i]]

    if(i == 1){
      stop=length(gs)
    }else{
      stop <- n + c(length(gs)-1)
    }

    items <- n:stop
    paths.core.genes[items] <- gs
    names(paths.core.genes)[items] <- names(gs)

    n <- stop+1
  }

  return(paths.core.genes)

}


#MUST FINISH
# GAGE.heatmap <- function(GAGE_GSA.res){
# library(pheatmap)
#   #result object from GAGE_GSA() function above
#
#
#   essSets.Up <- GAGE_GSA.res$essSets$essentialSets
#   up <- kegg$SigPaths[Kegg.essSets.Up[1:15],1:5]
#   up[,"q.val"]  <- round(up[,"q.val"], digits = 3)
#   up <- up[up[,"q.val"] <= 0.001,]
#
#   # dim(up) #11 with FDR < 0.001
#
#   dn <- kegg$gsa$less[Kegg.essSets.Dn$essentialSets[1:10],1:5]
#   dn[,"q.val"]  <- round(dn[,"q.val"], digits = 3)
#   dn <- dn[dn[,"q.val"] <= 0.001,]
#
#   head(dn)
#   dim(dn) #8 with FDR < 0.001
#
#
#   # col <-  colorRampPalette(rev(brewer.pal(n = 5, name ="RdYlBu")))(299)
#   len <- 299
#   col <- colorRampPalette(c("deepskyblue4", "deepskyblue3", "deepskyblue2", "deepskyblue1","white","red1", "red2", "red3", "red4"))(n=len)
#
#   paths <- c(rownames(up), rownames(dn))
#   stat <- kegg$gsa$stats[paths, 2:ncol(kegg$gsa$stats)] #ranges from -5.4 to 4.5
#
#
#   Breaks <- c(seq(min(stat), 0, length.out=ceiling(len/2) + 1),
#               seq(max(stat)/len, max(stat), length.out=floor(len/2)))
#
#   q.val <- as.character(c(up[,4],dn[,4]))
#   q.val <- ifelse(q.val!="0", paste("FDR =", q.val), "FDR < 0.001")
#   labs <- gsub("KEGG_", "", paths ) %>% gsub("_"," ", .) %>%
#     paste0(.,", ", q.val )
#
#   anno.col <- merged %>%
#     filter(TARGET.USI.1 %in% names(DEGs$cts.hd.1031$phenovector)) %>%
#     filter(CBFA2T3.GLIS2=="Yes") %>%
#     dplyr::select(TARGET.USI.1,M7_AML) %>%
#     column_to_rownames("TARGET.USI.1")

# pheatmap::pheatmap(mat=stat,
#                    color=col,
#                    scale = "none",
#                    breaks = Breaks,
#                    cluster_rows=FALSE,
#                    clustering_method="complete",
#                    clustering_distance_cols="euclidean",
#                    kmeans_k = NA,
#                    # annotation_col = anno.col,
#                    # annotation_colors = list(M7_AML=c("Yes"="firebrick4","No"="red")),
#                    gaps_row=11,
#                    # cutree_cols = 2,
#                    labels_row = labs,
#                    fontsize = 20,
#                    fontsize_number = 14,
#                    gp=gpar(col="blue"),
#                    show_colnames=FALSE,
#                    display_numbers=FALSE,
#                    number_format="%.1f",
#                    number_color="black")
#
#
#
# }

