#Jenny Smith

#Jan 5, 2018

#Purpose: To remove duplicate genes from a expression matrix. allowing for genes to be set as row names.

#
# Based on https://groups.google.com/forum/#!topic/rsem-users/W9RQrZIOzA4

#' Function for the TPM conversion.
#'
#' @param RPKM expression matrix in RPKM
#'
#' @returns matrix
#' @export
#'
#' @examples
#' ex <- c("TBD")
RPKM_to_TPM <- function(RPKM){
  #This is for use with one patient column at a time, so use apply() or sapply()
  conversionFactor <- sum(RPKM) / 1E6
  TPM <- RPKM / conversionFactor
  return(TPM)
}


#' Remove duplicate genes from an expression matrix
#'
#' @param count.matrix  matrix or dataframe with expression values. Genes as rownames, samples as columns
#' @param ID.map  ID.map is a ensemble ID to gene symbol dataframe with colnames gene_id and gene_name.
#' @param matrix.class boolean to return a matrix or dataframe object
#' @param rowname.GeneSym return gene symbols as the rownames. If false, will have ensembl IDs as rownames
#'
#' @returns an expression matrix or dataframe
#' @export
#'
#' @examples
#' ex <- c("TBD")
rmDups <- function(count.matrix, ID.map, matrix.class=TRUE,rowname.GeneSym=TRUE){
  #count.matrix is a matrix class or dataframe with expression values. Genes as rownames, patients as columns
  #ID.map is a ensemble ID to gene symbol map with colnames gene_id and gene_name.
  #rowname.GeneSym= TRUE return gene symbols as the rownames. If false, will have ensembl IDs as rownames

  df <- count.matrix %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene_id") %>%
    dplyr::left_join(., ID.map, by=c("gene_id")) %>%
    dplyr::select(gene_id, gene_name, everything()) %>%
    dplyr::filter(!grepl("_PAR_", gene_id)) #specific to gencode refs


  dup <- df$gene_name[which(duplicated(df$gene_name))]
  ddf <- df %>%
    dplyr::filter(gene_name %in% dup) %>%
    dplyr::arrange(gene_name)%>%
    dplyr::mutate(Variance=genefilter::rowVars(select(., -gene_id, -gene_name))) %>%

    dplyr::group_by(gene_name) %>%
    dplyr::mutate(High.Variance=Variance==max(Variance)) %>%
    dplyr::ungroup() %>%

    dplyr::filter(High.Variance)  %>%
    dplyr::filter(!duplicated(gene_name))#if variance tied, pick first one

  rmDups.df <- df %>%
    dplyr::filter(! gene_name %in% dup) %>%
    dplyr::bind_rows(.,ddf) %>%
    dplyr::select(everything(),-gene_id, -Variance,-High.Variance)

  if(rowname.GeneSym){
    rmDups.df <- rmDups.df %>%
      tibble::column_to_rownames("gene_name")
  }else{
    rmDups.df <- rmDups.df %>%
      dplyr::mutate(gene_id=gsub("\\.[0-9]{1,2}", "", gene_id)) %>%
      tibble::column_to_rownames("gene_id")
  }

  if(matrix.class){
    rmDups.df <- data.matrix(rmDups.df)
  }

  return(rmDups.df)
}


### To Do

# kallisto_rmDups <- function(kallisto.df, geneIDmap){
#
#   #Fix the PAR regions. Need to sum them since the fasta sequence is identical between PAR_Y and the X chromosome gene
#   PAR <- filter(geneIDmap, grepl("_PAR_Y", gene_id))
#   PAR <- PAR %>%
#     pull(gene_id) %>%
#     gsub("_PAR_Y", "", .) %>%
#     paste(., collapse="|")
#
#   #For each PAR region, use ColSums to combine them
#   idx <- grep(PAR, rownames(kallisto.df), value=TRUE)
#   length(idx) #90 PAR_Y and their X-chromosome pair
#
#   pairs <- lapply(seq(2, length(idx), by=2),
#                   function(x){
#                     rows <- idx[c((x-1):x)]
#                     colSums(kallisto.df[rows,])
#                   }) %>%
#     do.call(rbind, .)
#
#   #Remove the PAR_Y rows and add in the summed gene expression
#   non_PAR <- idx[seq(2, length(idx), by=2)-1]
#   kallisto.df <- kallisto.df[!grepl("_PAR_Y$", rownames(kallisto.df)), ]
#   kallisto.df[non_PAR,] <- pairs
#
#
#   # dim(kallisto.df) #59808  2116
#   # head(kallisto.df[,1:5])
#
#   #convert Ensembl IDs to Gene Symbols
#   kallisto.df.sym <- as.data.frame(kallisto.df) %>%
#     rownames_to_column("gene_id") %>%
#     left_join(., select(geneIDmap,gene_id,gene_name),
#               by="gene_id") %>%
#     select(gene_id,gene_name, everything())
#
#   #seperate out those genes without duplicate symbols
#   kallisto.df.sym_noDups <- kallisto.df.sym %>%
#     filter(!duplicated(gene_name), !duplicated(gene_name, fromLast = TRUE))
#
#   #identify which columns are numeric counts columns
#   numeric_cols <- sapply(kallisto.df.sym, function(x) is.numeric(x))
#   numeric_cols <- names(numeric_cols)[which(numeric_cols)]
#
#   #calculate the averate for all duplicate gene symbols.
#   kallisto.df.sym_Dups <- kallisto.df.sym %>%
#     filter(duplicated(gene_name) | duplicated(gene_name, fromLast = TRUE)) %>%
#     arrange(gene_name) %>%
#
#     #Calculate IQR to measure variance
#     rowwise() %>%
#     mutate(IQR=IQR(c_across(all_of(numeric_cols)))) %>%
#     ungroup() %>%
#
#     #use the max IQR expression to address duplicate gene symbols
#     group_by(gene_name) %>%
#     mutate(Rank=rank(IQR,ties.method= "first")) %>%
#     mutate(Keep=ifelse(Rank==max(Rank), TRUE, FALSE)) %>%
#     ungroup() %>%
#     select(IQR:Keep, everything())
#
#   # head(kallisto.df.sym_Dups[,1:10])
#   # dim(kallisto.df.sym_Dups) #1719 2118
#
#   #remove duplicated gene_names
#   kallisto.df.sym_rmDups <- kallisto.df.sym_Dups %>%
#     filter(Keep) %>%
#     select(-IQR,-Keep, -Rank)
#
#
#   # dim(kallisto.df.sym_rmDups) #174 2118
#   # head(kallisto.df.sym_rmDups)
#
#   kallisto.df.sym_final <- bind_rows(kallisto.df.sym_noDups, kallisto.df.sym_rmDups)
#
# }
