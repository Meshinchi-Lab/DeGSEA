#Jenny Smith

#Jan 5, 2018

#Purpose: To remove duplicate genes from a expression matrix. allowing for genes to be set as row names. 


getIDmap <- function(GTF,type="transcript"){
  #GTF is a dataframe from read.delim(gtf_file)
  #type is either transcript or gene
  library(dplyr)
  library(tibble)
  options(stringsAsFactors = FALSE)
  
  #standard ensembl GTF format and gencode GTF.
  df <- GTF %>%
    filter(grepl(type, V3)) %>% 
    dplyr::pull(V9) %>% #use pull() to create vector from a single column 
    str_split(., pattern = "; ") %>% 
    lapply(., function(x) t(str_split(x, pattern = " ", simplify = TRUE))) %>% 
    sapply(.,  function(x) set_colnames(x, value = x[1,])[-1,]) %>% #bapply ?
    sapply(., function(x) data.frame(as.list(x))) %>% 
    bind_rows(.) %>% 
    mutate(across(everything(), ~gsub("\"","",.x)))
  
  return(df)
}


#Function for the TPM conversion. 
# Based on https://groups.google.com/forum/#!topic/rsem-users/W9RQrZIOzA4
RPKM_to_TPM <- function(RPKM){
  #This is for use with one patient column at a time, so use apply() or sapply()
  conversionFactor <- sum(RPKM) / 1E6
  TPM <- RPKM / conversionFactor
  return(TPM)
}

rmDups <- function(count.matrix, ID.map, matrix.class=TRUE,rowname.GeneSym=TRUE){
  #count.matrix is a matrix class or dataframe with expression values. Genes as rownames, patients as columns 
  #ID.map is a ensemble ID to gene symbol map with colnames gene_id and gene_name.
  #rowname.GeneSym= TRUE return gene symbols as the rownames. If false, will have ensembl IDs as rownames
  
  df <- count.matrix %>%
    as.data.frame() %>%
    rownames_to_column("gene_id") %>%
    left_join(., ID.map, by=c("gene_id")) %>% 
    select(gene_id, gene_name, everything()) %>%
    filter(!grepl("_PAR_", gene_id)) #specific to gencode refs
  
  
  dup <- df$gene_name[which(duplicated(df$gene_name))]
  ddf <- df %>%
    filter(gene_name %in% dup) %>% 
    arrange(gene_name)%>% 
    mutate(Variance=genefilter::rowVars(select(., -gene_id, -gene_name))) %>%
    
    group_by(gene_name) %>% 
    mutate(High.Variance=Variance==max(Variance)) %>% 
    ungroup() %>% 
    
    filter(High.Variance)  %>%
    filter(!duplicated(gene_name))#if variance tied, pick first one
  
  
  if(rowname.GeneSym){
    rmDups.df <- df %>% 
      filter(! gene_name %in% dup) %>% 
      bind_rows(.,ddf) %>%
      select(everything(),-gene_id, -Variance,-High.Variance) %>%
      column_to_rownames("gene_name")
  }else{
    rmDups.df <- df %>% 
      filter(! gene_name %in% dup) %>% 
      bind_rows(.,ddf) %>%
      select(everything(),-gene_name, -Variance,-High.Variance) %>%
      mutate(gene_id=gsub("\\.[0-9]{1,2}", "", gene_id)) %>% 
      column_to_rownames("gene_id")
  }
  
  if(matrix.class){
    rmDups.df <- data.matrix(rmDups.df)
  }
  
  return(rmDups.df)
}  


#########Older Code from July 17, 2017. 
#NOTE:#this code is too memory intensive.  DOES NOT WORK FOR MANY NAs!!!  
rmDupGenes <-  function(expnData, geneCol,set.rownames=TRUE){
  library(genefilter)
  library(magrittr)
  #expnData is the matrix with a column for gene names.
  #geneCol a character vector with the name of the column for gene names.
  #DO NOT USE read_csv(). Tibble class causes the function, which used 100% base R, to crashhh. really dumb.
  #if using set.rowames=TRUE then it return gene symbols as rownames and all numeric (counts) columns. Else, return the filtered matrix with geneSymbol columns identical to input (allows for having more character columns in the ouput)
  idx <- which(duplicated(expnData[, geneCol]))
  dups <- unique(expnData[idx,geneCol])  
  
  #add regex anchors to duplicate gene names
  dup.regex <- paste0("^", dups, "$") 
  genes <- expnData[,geneCol]
  numeric <- sapply(expnData, is.numeric)
  vars <- genefilter::rowVars(expnData[,numeric])
  # vars <- genefilter::rowVars(expnData[,-(which(colnames(expnData) == geneCol))])
  names(vars) <- genes

  tokeep <- function(dup,vars){
    geneVar <-  vars[grepl(dup, names(vars))]
    dup <- gsub("\\^|\\$", "", dup) #strip the anchors
    keep <- which(names(vars) == dup & vars == max(geneVar))
    #if there are ties with same max variation.
    if (length(keep) > 1){
      keep <- keep[1]
    }
    return(keep)
  }

  keep <- sapply(dup.regex,tokeep,vars=vars)
  sel.Dups <- expnData[keep, ]
  noDups <- expnData[which(! expnData[,geneCol] %in% dups), ]
  # print(sapply(list(sel.dup, noDups), dim))

  remDups <- bind_rows(noDups, sel.Dups) 
  
  if(set.rownames){
    remDups <- remDups %>% 
      set_rownames(., .[,geneCol])
    remDups <- remDups[,numeric]
  }

  list <- list(dups,dup.regex, vars, keep, remDups)
  names(list) <- c("dups", "dup.regex", "vars", "keep", "remDups")
  return(list)
}

#Function to filter out genes with low expression. 
filterGs <- function(expnMatrix, cutoff,log2=FALSE){
  #genes are  rownames in expnMatrix
  #cutoff is  numeric value
  filtered <- expnMatrix %>%
    rownames_to_column() %>%
    filter(rowSums(select(.,-rowname)) > cutoff) 
  
  if(log2){
    filtered <- filtered %>% 
      mutate_if(is.numeric, function(x) log2(x+1)) 
  }
  
  filtered <- column_to_rownames(filtered, "rowname")
}


kallisto_rmDups <- function(kallisto.df, geneIDmap){
  
  #Fix the PAR regions. Need to sum them since the fasta sequence is identical between PAR_Y and the X chromosome gene 
  PAR <- filter(geneIDmap, grepl("_PAR_Y", gene_id))
  PAR <- PAR %>% 
    pull(gene_id) %>% 
    gsub("_PAR_Y", "", .) %>% 
    paste(., collapse="|")
  
  #For each PAR region, use ColSums to combine them 
  idx <- grep(PAR, rownames(kallisto.df), value=TRUE) 
  length(idx) #90 PAR_Y and their X-chromosome pair
  
  pairs <- lapply(seq(2, length(idx), by=2), 
                  function(x){ 
                    rows <- idx[c((x-1):x)]
                    colSums(kallisto.df[rows,])
                  }) %>% 
    do.call(rbind, .)
  
  #Remove the PAR_Y rows and add in the summed gene expression
  non_PAR <- idx[seq(2, length(idx), by=2)-1]
  kallisto.df <- kallisto.df[!grepl("_PAR_Y$", rownames(kallisto.df)), ]
  kallisto.df[non_PAR,] <- pairs
  
  
  # dim(kallisto.df) #59808  2116
  # head(kallisto.df[,1:5])
  
  #convert Ensembl IDs to Gene Symbols
  kallisto.df.sym <- as.data.frame(kallisto.df) %>%
    rownames_to_column("gene_id") %>%
    left_join(., select(geneIDmap,gene_id,gene_name),
              by="gene_id") %>%
    select(gene_id,gene_name, everything())
  
  #seperate out those genes without duplicate symbols
  kallisto.df.sym_noDups <- kallisto.df.sym %>%
    filter(!duplicated(gene_name), !duplicated(gene_name, fromLast = TRUE))
  
  #identify which columns are numeric counts columns
  numeric_cols <- sapply(kallisto.df.sym, function(x) is.numeric(x))
  numeric_cols <- names(numeric_cols)[which(numeric_cols)]
  
  #calculate the averate for all duplicate gene symbols.
  kallisto.df.sym_Dups <- kallisto.df.sym %>%
    filter(duplicated(gene_name) | duplicated(gene_name, fromLast = TRUE)) %>%
    arrange(gene_name) %>%
    
    #Calculate IQR to measure variance
    rowwise() %>%
    mutate(IQR=IQR(c_across(all_of(numeric_cols)))) %>%
    ungroup() %>%
    
    #use the max IQR expression to address duplicate gene symbols
    group_by(gene_name) %>%
    mutate(Rank=rank(IQR,ties.method= "first")) %>%
    mutate(Keep=ifelse(Rank==max(Rank), TRUE, FALSE)) %>%
    ungroup() %>%
    select(IQR:Keep, everything())
  
  # head(kallisto.df.sym_Dups[,1:10])
  # dim(kallisto.df.sym_Dups) #1719 2118
  
  #remove duplicated gene_names
  kallisto.df.sym_rmDups <- kallisto.df.sym_Dups %>%
    filter(Keep) %>%
    select(-IQR,-Keep, -Rank)
  
 
  # dim(kallisto.df.sym_rmDups) #174 2118
  # head(kallisto.df.sym_rmDups)
  
  kallisto.df.sym_final <- bind_rows(kallisto.df.sym_noDups, kallisto.df.sym_rmDups)
  
  # dim(kallisto.df.sym_final) #58263  2118
  # head(kallisto.df.sym_final[,1:5])
}


# Compatible with 60,000 by 20,000 size dataframes. and compatible with dplyr
#Not faster, despite having fewer intermediate dataframes in memory and instead is a vector of T/F
#with one subset of the entire df, the o

# filterDups <- function(df, geneCol){
#   #Function to input a expression matrix and remove duplicate gene entries
#   #df is the expression matrix with a column for gene names as a character class.
#   #geneCol is the character vector with the column name.
#   #best usage: Arrange the dataframe by gene name BEFORE using this function.
#   library(magrittr)
#   library(genefilter)
#   library(dplyr)
# 
#   #Genes and duplicate genes in the dataframe
#   genes <-unlist(df[,geneCol])
#   dup.genes <- genes[duplicated(genes) | duplicated(genes, fromLast=TRUE)]
# 
#   #vector of all TRUEs to be the base vector to be updated based on rowVars.
#   keep <- rep(TRUE, nrow(df)) %>%
#     set_names(genes)
# 
#   #subset the data frame for only duplicate genes
#   df.dup <- df %>%
#     filter_(paste(geneCol,"%in%","dup.genes"))
# 
#   #remove the larger dataframe from RAM
#   rm(df)
# 
#   #For loop to determine which duplicate gene has the largest variance.
#   for (dup in unique(dup.genes)){
# 
#     # print(dup)
#     temp <- df.dup %>%
#       filter_(paste(geneCol, "==", paste0('"', dup,'"'))) %>%
#       # filter_(paste(geneCol, "==", paste0("\'",dup,"'"))) %>% #noticed possibly typo. fixed above
#       select(-which(colnames(.) == geneCol))
# 
#     geneVars <- rowVars(temp)
# 
#     #if the variance is the same for all duplicate genes, pick the first one.
#     if (max(geneVars) - min(geneVars) == 0){
#       # print("equal var")
#       update <- c(TRUE, rep(FALSE, nrow(temp)-1))
#     }else{
#       #else, select the gene with highest variation.
#       update <- apply(temp, 1, function(x) var(x) >= max(geneVars))
# 
#       #more than one duplicate with the same max variance, select the first one.
#       if (sum(update) > 1){
#         update[which(update)]  <- c(TRUE,rep(FALSE, length(which(update))-1))
#       }
#     }
# 
#     keep[names(keep) == dup] <- update
#   }
# 
# 
#   return(keep)
# }



