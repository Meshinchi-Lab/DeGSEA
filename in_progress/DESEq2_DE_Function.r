#Jenny Smith

#February 13,  2017

#Purpose: Perform Differential Expression analysis with two group comparisons. 


# setwd("/home/jlsmith3/RNA_seq_Analysis/2017.04.21_TARGET_AML_IF")
# counts <- get(load("RNAseq_Counts_HD_LD_batchCorrected.RData"))
# subsetCounts <- get(load("PatientGroups_HD_LD_batchCorrected.RData"))


GroupIDs <- function(clinData, col){
  #clindata has patient IDs as rownames. 
  #col is a chracter string of the column name in the clinical data with the factor/variable information.
  list <- list()
  grps <- unique(clinData[,col])
  N <- length(grps)
  
  for (i in 1:N){
    IDs <- rownames(clinData[grepl(grps[i], clinData[,col]), ])
    list[[i]] <- IDs
  }
  names(list) <- grps
  return(list)
}


#phenovectors is my term for what is really "coldata" in DESeq
phenoVectors <- function(groupA, groupB){
  library(magrittr)
  #groupA and GroupB are character vectors with the patients IDs in each group
  g1 <- as.character(substitute(groupA))
  g2 <- as.character(substitute(groupB)) 
  
  vector <- c(rep(g1, length(groupA)), rep(g2, length(groupB)))
  names(vector) <- c(groupA, groupB)
  
  return(vector)
}

# 
# DEGs_DESeq <- function(expnData, colData, colNames, reference){
#   #expndata is a matrix of counts, with patient IDs are colnames, genes as rows
#   #coldata is a dataframe with patient IDs as rownames, and variables as colnames (eg, mutation, flt3)
#   #reference is a charater vector with the group that is the comparison reference eg "neg". 
#   
#   library(DESeq2)
#   library(magrittr)
#   
#   expnData <- as.matrix(round(expnData, digits=0)) #if input is fractional counts
#   expnData <- expnData[,names(coldata)] #correct order and subset
#   # groups <- as.data.frame(groups) %>% setNames(. , "Status") #make into dataframe with status column
#   
#   dds <- DESeqDataSetFromMatrix(countData = expnData, 
#                                 colData = colData, 
#                                 design = ~ Status)
#   
#   dds$Status <- relevel(dds$Status, ref=reference) #FC will be based on reference=factor1/ factor2
#   dds <- dds[ rowSums(counts(dds)) > 1, ]
#   dds <- DESeq(dds)
#   
#   vst <- vst(dds)
#   p <- plotPCA(vst, intgroup=colnames(colData))
#   
#   res <- results(dds, pAdjustMethod = "BH")
#   res2 <- res[abs(res$log2FoldChange) > 1, ]
#   res2 <- subset(res2, padj < 0.05)
#   
#   list <- list(dds,vst, p, res, res2)
#   names(list) <- c("dds","vst","PCA","results", "Filt_Res")
#   
#   return(list)
# }


DEGs_DESeq <- function(expnData, groups, reference){
  #expndata is a matrix of counts, with patient IDs are colnames, genes as rows
  #groups is a character vector with Patient IDs as rownames, and status, eg. pos,neg, etc. for each patient.
  #reference is a charater vector with the group that is the comparison reference eg "neg".

  library(DESeq2)
  library(magrittr)

  expnData <- as.matrix(round(expnData, digits=0)) #if input is fractional counts
  expnData <- expnData[,names(groups)] #correct order and subset
  groups <- as.data.frame(groups) %>% setNames(. , "Status") #make into dataframe with status column

  dds <- DESeqDataSetFromMatrix(countData = expnData,
                                colData = groups,
                                design = ~ Status)

  dds$Status <- relevel(dds$Status, ref=reference) #FC will be based on reference=factor1/ factor2
  dds <- dds[ rowSums(counts(dds)) > 10, ]
  dds <- DESeq(dds)

  vst <- vst(dds)
  p <- plotPCA(vst, intgroup= "Status")

  res <- results(dds, pAdjustMethod = "BH")
  res2 <- res[abs(res$log2FoldChange) > 1, ]
  res2 <- subset(res2, padj < 0.05)

  list <- list(dds,vst, p, res, res2)
  names(list) <- c("dds","vst","PCA", "results", "Filt_Res")

  return(list)
}






