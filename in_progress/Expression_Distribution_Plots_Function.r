#Jenny Smith

#May 5, 2017

#purpose: A function to output density plots, boxplots, and waterfall plots for a specific set of genes. 



expnDist <- function(expnMatrix,phenoVector, geneList,unit, BM=FALSE,logCPM=NULL,plot=FALSE){
  #expnMatrix is a dataframe with genes are rownames, patient IDs as col names
  #phenoVector is character vector with group membership (eg pos, neg) with patient IDs as names
  #genelist is a character vector
  #BM is for whether to include BM samples into phenovector
  #logCPM is whether to log2 CPM normalize the expn data
  library(magrittr)
  library(ggplot2)
  library(plyr)
  
  if(BM == TRUE){
    BM <- rep("BM", length(grep("^BM", colnames(expnMatrix)))) %>% setNames(grep("^BM", colnames(expnMatrix), value = TRUE))
    phenoVector <- c(phenoVector, BM)
  }else if (BM==FALSE){
    phenoVector = phenoVector
  }

  expnMatrix <- expnMatrix[rownames(expnMatrix) %in% geneList, ] #subset for genes of interest
  expnMatrix <- expnMatrix[, intersect(names(phenoVector), colnames(expnMatrix))] #match the column names 

  
  if (is.null(logCPM)){
    expnMatrix = expnMatrix
    x <- unit
    y <- unit
  } else if (logCPM==TRUE){
    expnMatrix <- cpm(expnMatrix, log=TRUE, prior.count = 1)  #convert to log2 CPM
    x <- "Log2 Counts per Million (CPM)"
    y <- "Log2 Counts per Million (CPM)"
  } else if (logCPM == FALSE){
    expnMatrix <- apply(expnMatrix, 2, function(x) log2(x + 1))
    x <- paste("Log2", unit, sep=" ")
    y <- paste("Log2", unit, sep=" ")
  }
  
  tmp <- data.frame(t(expnMatrix),
                    Status=phenoVector) 
  
  if (plot==TRUE) {
    for (i in 1:nrow(expnMatrix)){
      gene <- rownames(expnMatrix)[i]
      means <- as.data.frame(tapply(tmp[,i], INDEX = tmp$Status, FUN=mean)) %>% cbind(., Status=rownames(.))
      
      dtitle <- paste("Density  Plot of ", gene, "Expression in TARGET AML")
      btitle <- paste("Distribution of ", gene, "Expression in TARGET AML")
      
      densityPlot <- ggplot(tmp, aes(x=tmp[,i], fill=Status)) +
        geom_histogram(aes(y=..density..), alpha=0.65, position="identity") +
        geom_density(alpha=0.5) +
        labs(title=dtitle, y="Density", x=x) +
        geom_vline(data=means, aes(xintercept = means[,1], color=Status),linetype="dashed") +
        theme_bw()
      
      boxPlot <- ggplot(tmp, aes(x=Status, y=tmp[,i], fill=Status)) +
        geom_boxplot() +
        labs(title=btitle, x=" ", y=y) +
        theme_bw()
      
      print(densityPlot)
      print(boxPlot)
    }
  }else if (plot==FALSE){
    return(tmp)
  }
}





# GroupIDs <- function(clinData, col){
#   #clindata has patient IDs as rownames. 
#   #col is a chracter string of the column name in the clinical data with the factor/variable information. 
#   list <- list()
#   grps <- unique(clinData[,col])
#   N <- length(grps)
#   
#   for (i in 1:length(grps)){
#     if (grepl("[^a-zA-Z0-9 :]", grps[i])){
#       grps[i] <- gsub("[^a-zA-Z0-9 :]", "\\.", grps[i]) #remove special characters and replace with a "."
#     }
#     IDs <- rownames(clinData[grepl(grps[i], clinData[,col]), ])
#     list[[i]] <- IDs
#   }
#   
#   names(list) <- grps
#   return(list)
# }
# 
# phenoVectors_MultipleGroups <- function(listofgoupsIDs){
#   library(magrittr)
#   #listofgoupsIDs contains a list with each item containing the IDs for each factor level.  
#   #See GroupIDs function - this produced the input called  "listofgroupIDs"
#   group <- names(listofgoupsIDs)
#   
#   vector <- NULL
#   names <- NULL
#   for (i in 1:length(listofgoupsIDs)){
#     g <- group[i]
#     vector <- c(vector, rep(g, length(listofgoupsIDs[[i]])))
#     names <- c(names, listofgoupsIDs[[i]])
#   }
#   
#   names(vector) <- names
#   return(vector)
# }
# 
# phenoVectors <- function(groupA, groupB){
#   library(magrittr)
#   #groupA and GroupB are character vectors with the patients IDs in each group
#   g1 <- as.character(substitute(groupA))
#   g2 <- as.character(substitute(groupB)) 
#   
#   vector <- c(rep(g1, length(groupA)), rep(g2, length(groupB)))
#   names(vector) <- c(groupA, groupB)
#   
#   return(vector)
# }



# expnDist <- function(expnMatrix,phenoVector, geneList,BM=FALSE,logCPM=FALSE){
#   #expnMatrix is a dataframe with genes are rownames, patient IDs as col names
#   #phenoVector is character vector with group membership (eg pos, neg) with patient IDs as names
#   #genelist is a character vector
#   #BM is for whether to include BM samples into phenovector
#   #logCPM is whether to log2 CPM normalize the expn data
#   
#   library(ggplot2)
#   library(plyr)
#   
#   if(BM == TRUE){
#     BM <- rep("BM", length(grep("^BM", colnames(expnMatrix)))) %>% setNames(grep("^BM", colnames(expnMatrix), value = TRUE))
#     phenoVector <- c(phenoVector, BM)
#   }else if (BM==FALSE){
#     phenoVector = phenoVector
#   }
#   
#   expnMatrix <- expnMatrix[rownames(expnMatrix) %in% geneList, ] #subset for genes of interest
#   expnMatrix <- expnMatrix[, intersect(names(phenoVector), colnames(expnMatrix))] #match the column names 
#   
#   if (logCPM==TRUE){
#     expnMatrix <- cpm(expnMatrix, log=TRUE, prior.count = 1)  #convert to log2 CPM
#     x <- "Log2 Counts per Million (CPM)"
#     y <- "Log2 Counts per Million (CPM)"
#   } else if (logCPM == FALSE){
#     expnMatrix <- apply(expnMatrix, 2, function(x) log2(x + 1))
#     x <- "Log2 Read Counts"
#     y <- "Log2 Read Counts"
#   }
#   
#   
#   tmp <- data.frame(t(expnMatrix),
#                     Status=phenoVector) 
# 
#   for (i in 1:nrow(expnMatrix)){
#     gene <- rownames(expnMatrix)[i]
#     means <- as.data.frame(tapply(tmp[,i], INDEX = tmp$Status, FUN=mean)) %>% cbind(., Status=rownames(.))
#     
#     dtitle <- paste("Density  Plot of ", gene, "Expression in TARGET AML")
#     btitle <- paste("Distribution of ", gene, "Expression in TARGET AML")
#     
#     densityPlot <- ggplot(tmp, aes(x=tmp[,i], fill=Status)) +
#       geom_histogram(aes(y=..density..), alpha=0.65, position="identity") +
#       geom_density(alpha=0.5) +
#       labs(title=dtitle, y="Density", x=x) +
#       geom_vline(data=means, aes(xintercept = means[,1], color=Status),linetype="dashed") +
#       theme_bw()
#     
#     boxPlot <- ggplot(tmp, aes(x=Status, y=tmp[,i], fill=Status)) +
#       geom_boxplot() +
#       labs(title=btitle, x="Group", y=y ) +
#       theme_bw()
#     
#     print(densityPlot)
#     print(boxPlot)
#     
#     name1 <- gsub("-", "\\.", gene) %>% paste(., "densityPlot.pdf", sep="_")
#     name2 <- gsub("-", "\\.", gene) %>% paste(., "boxPlot.pdf", sep="_")
#     
#   }
# }
