#Jenny Smith

#Differential Expression Analysis pipeline

#Purpose: Given matrix of RNA-seq raw counts and clinical annoations,
#create a series of differential expression Analyses for two group comparisions.
#Annotate DEGs with TM information, provide Heatmaps and PCA/MDS clustering.
#Heatmaps have option for additional color bars


#Function to create a gene ID map for ensembl gene_id, and transcript_ids
#' Create a gene ID map for ensembl gene_id, and transcript_ids
#'
#' @param GTF a gtf file
#'
#' @returns a dataframe
#' @export
#'
#' @examples
#' gtf <- data.frame(V3 = "transcript", V9 = "attr1 red; attr2 blue; attr green")
#' getIDmap(gtf)
getIDmap <- function(GTF){
  #standard ensembl GTF format and gencode GTF.
  tx <- GTF %>%
    dplyr::filter(grepl("transcript", V3)) %>%
    dplyr::select(V9) %>%
    unlist(.) %>%
    stringr::str_split(., pattern = "; ") %>%
    lapply(., function(x) t(stringr::str_split(x, pattern = " ", simplify = TRUE))) %>%
    sapply(.,  function(x) set_colnames(x, value = x[1,])[-1,]) %>%
    sapply(., function(x) data.frame(as.list(x))) %>%
    dplyr::bind_rows(.)

  return(tx)
}


#' Collapse duplicate rows so that all information is retain,optionally unique, and separated by a semi-colon.
#'
#' @param col colname
#' @param uniq boolean
#' @param split boolean
#' @param sep delimiter
#'
#' @returns string vector
#' @export
#'
#' @examples
#' metadata <- data.frame(mutation = rep("yes","no", length.out = 10))
#' collapseRows(col = "mutation")
collapseRows <- function(col, uniq=FALSE, split=FALSE,sep=""){
  #designed for dplyr so that "col" paramter is a vector of that column.
  #Similar to collapseDuplicates(), but for fewer columns, plus preselection of columns. Where are collapseDuplicates() you don't need to know the exact column names before hand.

  if (uniq){
    if(split){
      col <- stringr::str_split(col, pattern = sep) %>% unlist() %>% unique()
    }else{
      col <- unique(col)
    }
  }
  collapsed <- ifelse(all(is.na(col)), NA, paste(col, collapse = "; "))
  return(collapsed)
}



#' Voom DE with Batch effect in the Model
#'
#' @param expnData gene counts matrix
#' @param clinData dataframe - has patients as rownames
#' @param col a vector of 2 column names for co-variate/batch effect variable
#' @param percent percent of samples for expression threshold
#' @param trend boolean limma trend
#' @param logCPM boolean - true if input counts are already log2 CPM
#' @param normalization boolean to run upper quantile norm
#' @param GOI vector of gene names - genes of interest (GOI)
#'
#' @returns list
#' @export
#'
#' @examples
#' mat <- matrix(rnorm(10), nrow = 10, ncol = 10)
#' metadata <- data.frame(mutation = rep("yes","no", length.out = 10), group=rep("A","B", each = 2))
#' voom_DE_BE(mat, metadata, col = c("mutation","group"))
voom_DE_BE <- function(expnData,clinData,col,percent=0.05,
                       trend=FALSE, logCPM=FALSE,
                       normalization=FALSE,
                       GOI=NULL, exclude_samples = NULL, scale_data = FALSE){
  # library(edgeR)
  # library(limma)

  ##ensure correct order
  expnData <- expnData[,match(rownames(clinData), colnames(expnData))]

  if (!all(complete.cases(expnData))){
    print("Names DO NOT match in between phenovector and colnames of expression matrix")
    return(list(expnData=expnData,pheno=clinData))
  }

  #NOTE ClinData MUST BE already factor  leveled!
  dge <- edgeR::DGEList(counts = expnData, samples = clinData[,col])

  if( ! is.null(exclude_samples)){
    #updated on 10/11/17 to add AML samples CPM cutoff before calc. norm factors.
    # example: selected_samps = "BM[0-9]|R[O0][0-9]"
    # selected_samps <- ! grepl(exclude_samples, colnames(expnData))
    # AML <- ! grepl("BM[0-9]|R[O0][0-9]", colnames(expnData))
    # AMLsamples <- sum(AML)
    selected_samps <- ! grepl(exclude_samples, colnames(expnData))
  } else {
    selected_samps <- colnames(expnData)
  }


  #remove low count genes. This is to allow for 2/3 samples must have 1 cpm. Since 3 is the minimum # of samples I will allow for DE analysis.
  N <- ncol(expnData[,selected_samps])
  keep.dge <- rowSums(edgeR::cpm(dge)[,selected_samps] >= 1) >= max(2,(percent*N))
  dge <- dge[keep.dge,] #subset for those genes with cmp >= 1 per gene in AML samples
  dge <- edgeR::calcNormFactors(dge) #Do TMM normalization

  #Create a design and contrasts with  the groups to compare and the co-variate/batch effect variable
  #since the columns in dge$samples dataframe are already factor leveled - the ref is in the first column in the design
  design <- stats::model.matrix(formula(paste(c("~0", col), collapse="+")),
                         data=dge$samples) #~0 means no intercept.
  colnames(design) <- c("Ref","Comparitor","BatchEffect")


  #contrast is the comparison-referenceGroup, or Pos-Neg, or Mut-WT
  cont.matrix <- limma::makeContrasts(contrasts = paste(c("Comparitor","Ref"),collapse = "-"),
                               levels = design)


  if (is.null(GOI)){
    GOI <- 1:nrow(dge)
  }else{
    GOI <- intersect(rownames(dge), GOI)
    print(paste0("Length of GOI: ", length(GOI)))
  }


  if (logCPM){
    dge.norm <- edgeR::cpm(dge, log=TRUE, prior.count = 1) #log2 + 1 CPM transformatio for normalization for sample to sample comparisons.
    norm_method <- "TMMCPM"

  }else if (all(!logCPM & !normalization)){
    dge.norm <- limma::voom(dge, design, plot = TRUE) #can I use voom transformed values like CPM? yes
    norm_method <- "Voom"  #voom transformed counts for sample to sample comparisons.

  }else if(all(!logCPM & normalization=="qt")){
    dge.norm <- limma::voom(dge, design, plot = FALSE, normalize.method = "quantile")
    norm_method <- "Voom.quantile" #voom and quantilte normalized counts for sample to sample comparisons.
  }

  print(norm_method) #to confirm which type of DE is performed, trend or voom

  #fit the linear model.
  fit <- limma::lmFit(dge.norm, design)
  fit <- limma::contrasts.fit(fit, contrasts = cont.matrix)

  #compute moderated t-statistics using empirical bayes moderation.
  if(all(trend & logCPM)){ #only use limma trend method with CPM values, as per manual.
    fit2 <- limma::eBayes(fit,trend = trend)[GOI,]
  }else{
    fit2 <- limma::eBayes(fit)[GOI,]
  }

  # select differentially expressed genes.
  DE <- limma::topTable(fit2,adjust.method="BH",sort.by="P",
                  number=nrow(dge),
                  p.value=0.05, lfc=1) #abs(logFC) >= 1 for all genes

  list <- list(dge.norm,fit2, DE)
  names(list) <- c(norm_method,"eBayesFit", "DE")

  return(list)
}


#
#' Function for Limma Voom differential expression
#'
#' @param expnData a matrix or data frame with the raw counts. Patient IDs as colnames, genes as rownames
#' @param pheno a character vector with patient IDs as names, and the status for each in each group(eg pos,neg)
#' @param ref  a character vector of the reference level for DE. for example ref="No".
#' @param percent percent is the fraction (0-1 numberic) of AML samples to include when setting an expression threshold. eg 5% of AMLs, percent=0.05
#' @param logCPM boolean. use logCPM values or not.
#' @param trend trend is for using limma trend method with log2 CPMs
#' @param normalization method of normalization such as quantile if necessary. should be either FALSE or "qt" so far
#' @param GOI a character vector of genes (or numeric vector of row indices) of interest to subset at the end. keeps BH adjuted p-values more accurate.
#' @param eBayesRobust boolean
#' @param lmMethod least squares (ls)
#'
#' @return list
#' @export
#'
#' @examples
#' ex <- c('TBD')
#'
voom_DE <- function(expnData, pheno,ref,percent,
                    logCPM=FALSE,trend=FALSE,
                    normalization=FALSE,GOI=NULL,
                    eBayesRobust=FALSE,
                    lmMethod="ls", exclude_samples = NULL, scale_data = FALSE) {
  # expnData is a matrix or data frame with the raw counts. Patient IDs as colnames, genes as rownames
  # pheno is a character vector with patient IDs as names, and the status for each in each group(eg pos,neg)
  #ref is a chacter vector of the reference level for DE. for example ref="No".
  # percent is the fraction (0-1 numberic) of AML samples to include when setting an expression threshold. eg 5% of AMLs, percent=0.05.
  #trend is for using limma trend method with log2 CPMs
  #normalization is for an extra method of normalization such as quantile if necessary. should be either FALSE or "qt" so far
  #GOI is a character vector of genes (or numeric vector of row indices) of interest to subset at the end. keeps BH adjuted p-values more accurate.

  expnData <- expnData[,match(names(pheno), colnames(expnData))] ##ensure correct order

  if (!all(complete.cases(expnData))){
    print("Names DO NOT match in between phenovector and colnames of expression matrix")
    return(list(expnData=expnData,pheno=pheno))
  }

  groups <- unique(pheno)
  groups <- c(groups[groups != ref], ref) #order so that reference is second
  pheno.f <- factor(pheno, levels=groups)
  dge <- edgeR::DGEList(counts = expnData, group = pheno.f)

  if( ! is.null(exclude_samples)){
    #updated on 10/11/17 to add AML samples CPM cutoff before calc. norm factors.
    # example: selected_samps = "BM[0-9]|R[O0][0-9]"
    # selected_samps <- ! grepl(exclude_samples, colnames(expnData))
    # AML <- ! grepl("BM[0-9]|R[O0][0-9]", colnames(expnData))
    # N <- sum(AML)
    selected_samps <- ! grepl(exclude_samples, colnames(expnData))
  } else {
    selected_samps <- colnames(expnData)
  }

  #This is to allow for 2/3 samples must have 1 cpm. Since 3 is the minimum # of samples I will allow for DE analysis.
  N <- ncol(expnData[,selected_samps])
  keep.dge <- rowSums(edgeR::cpm(dge)[,selected_samps] >= 1) >= max(2,(percent*N)) #X% of AML samples has cpm of at least 1 for a gene
  dge <- dge[keep.dge,] #subset for those genes with cmp >= 1 per gene in AML samples
  dge <- edgeR::calcNormFactors(dge) #Do TMM normalization

  design <- stats::model.matrix(~0 + pheno.f, data=dge$samples)#~0 means no intercept.
  colnames(design) <- levels(pheno.f)

  cont.matrix <- limma::makeContrasts(contrasts = paste(groups, collapse = "-"), levels = design) #contrast is approx. log2(mean(Pos)) - log2(mean(Neg)) per gene.

  if (is.null(GOI)){
    GOI <- 1:nrow(dge)
  }else{
    GOI <- intersect(rownames(dge), GOI)
    print(paste0("Length of GOI: ", length(GOI)))
  }


  if (logCPM){
    dge.norm <- edgeR::cpm(dge, log=TRUE, prior.count = 1) #log2 + 1 CPM transformatio for normalization for sample to sample comparisons.
    norm_method <- "TMMCPM"

  }else if (all(!logCPM & !normalization)){
    dge.norm <- limma::voom(dge, design, plot = FALSE) #can I use voom transformed values like CPM? yes
    norm_method <- "Voom"  #voom transformed counts for sample to sample comparisons.

  }else if(all(!logCPM & normalization=="qt")){
    dge.norm <- limma::voom(dge, design, plot = FALSE, normalize.method = "quantile")
    norm_method <- "Voom.quantile" #voom and quantilte normalized counts for sample to sample comparisons.
  }

  print(norm_method) #to confirm which type of DE is performed, trend or voom

  #fit the linear model.
  print(lmMethod)
  fit <- limma::lmFit(dge.norm, design,method=lmMethod)
  fit <- limma::contrasts.fit(fit, contrasts = cont.matrix)

  #compute moderated t-statistics using empirical bayes moderation.
  if(all(trend & logCPM)){ #only use limma trend method with CPM values, as per manual.
    fit2 <- limma::eBayes(fit,trend = trend, robust=eBayesRobust)[GOI,]
  }else{
    fit2 <- limma::eBayes(fit, robust=eBayesRobust)[GOI,]
  }

  # select differentially expressed genes.
  DE <- limma::topTable(fit2,adjust.method="BH",sort.by="P",
                 number=nrow(dge),p.value=0.05, lfc=1) #abs(logFC) >= 1 for all genes

  list <- list(dge.norm,fit2, DE)
  names(list) <- c(norm_method,"eBayesFit", "DE")


  return(list)
}


################### Pipeline DE Analysis Function #########################

#' Pipeline DE Analysis Function
#'
#' @param expnData a matrix or data frame with the raw counts. Patient IDs as colnames, genes as rownames
#' @param clinData data.frame of clindata has sample_id as rownames. must have column 'sample_id'
#' @param col a character string of the factor column of interest
#' @param ref the character strign of the reference group level (eg NBM, Neg, or control)
#' @param percent.cutoff the fraction (0-1 numberic) of AML samples to include when setting an expression threshold. eg 5% of AMLs, percent=0.05.
#' @param logCPM boolean whether use logCPM values
#' @param BM boolean whether to compare group to normal bone marrows
#' @param GOI a character vector of genes (or numeric vector of row indices) of interest to subset at the end. keeps BH adjuted p-values more accurate.
#' @param anno boolean whether to add the annotation of the DEGs with TM and cell localization patterns.
#' @param ids2symbols ids2symbols can be NULL, NA, or a file path. NULL is for BCCA id mapping + backwards compatibility.
#' @param gene.name.col a character vector for column name containing gene_names if a dataframe has an alternative column name for gene symbols (eg BCL2,TP53). This CANNOT be gene_id.
#' @param method hierarchical clustering method
#' @param Add.Anno.Col character vector of additional column names from clinData to add one or two more columns for heatmap annotation color bars
#' @param Custom.Cols use entirely custom annotation cols
#' @param SkipPlots boolean whether to generate PCA, MDS, and Heatmaps
#'
#' @return list
#' @export
#'
#' @examples
#'
#' ex <- c('TBD')
#'
twoGroups_DEGs <- function(expnData, clinData, col, ref,
                           percent.cutoff=0.05,logCPM=FALSE,
                           GOI=NULL,
                           anno=TRUE, ids2symbols=NULL,
                           gene.name.col="gene",
                           method="ward.D2",
                           Add.Anno.Col=NULL,
                           Custom.Cols=NULL,
                           SkipPlots=FALSE,
                           BM=FALSE){
  # expnData is a matrix or data frame with the raw counts. Patient IDs as colnames, genes as rownames
  #clindata has patient IDs as rownames.
  #col is a character string of the factor column of interest
  #ref is the character strign of the reference group level (eg BM, Neg, or control)
  #anno is for the annotation of the DEGs with TM and cell localization patterns.
  #Add.Anno.Col is to add one or two more columns for heatmap annotation color bars
  #Anno.Cols is to use entirely custom annotation cols.

  #remove unknown categories from the datasets since only want yes/no or 0/1 groups
  rmUnknowns <- function(clinData, cols){
    removeUnknowns <- clinData

    for (i in 1:length(cols)){
      idx <- ! grepl("Unknown",removeUnknowns[, cols[i]], ignore.case=TRUE)
      removeUnknowns <- removeUnknowns[idx, ]
    }
    return(removeUnknowns)
  }

  # For plot file names
  variantName <- col

  # Remove unknowns from clindata
  clinData <- rmUnknowns(clinData, col)

  if (BM){
    ref <- "NBM"
    normalBMs <- grep("BM[0-9]|R[O0][0-9]", colnames(expnData), value=TRUE)
    N_normalBMs <- length(normalBMs)
    clinData <- clinData %>%
      tibble::add_row( {{ col }} := rep("NBM", N_normalBMs),
                       sample_id = normalBMs)
  }
  # Define Groups to compare based on group IDs from clinical data.
  groups <- clinData %>%
    dplyr::pull(dplyr::all_of(col), name = sample_id)
  # Intersect with expression matrix to subset.
  groups <- groups[base::intersect(names(groups), colnames(expnData))]

  # Only analyze at least 3x3 comparisons at minimum
  if ( any(table(groups) < 3) ){
    message("Fewer than 3 samples in one of the conditions")
    list <- list(expnData, clinData, groups)
    names(list) <- c("InputExpnMatrix","InputClinData", "Groups")
    return(list)
  }

  #Define the phenotype vector
  phenoVector <- groups %>%
    base::as.factor(.) %>%
    stats::relevel(.,ref = ref)
  phenoVector <- phenoVector[order(phenoVector)]

  #subset and order the dataframe.
  expnData <- expnData[ ,names(phenoVector)] #mutant (groupA), then WT (groupB)
  clinData <- clinData[names(phenoVector), ]

  if (any(is.na(expnData))){
    message("NAs Introduced. Check Rownames and colnames of inputs")
    list <- list(expnData, clinData, groups)
    names(list) <- c("InputExpnMatrix","InputClinData", "Groups")
    return(list)
  }

  # Calculate Differential Expression
  DE <- voom_DE(expnData = expnData,
                pheno = phenoVector, #mutant - wild type. logCPM and voom transform the counts
                ref=ref,
                percent=percent.cutoff,
                logCPM=logCPM,
                normalization=FALSE,
                GOI=GOI)

  #add column for linear scale fold-changes
  DE$DE$FoldChange <- gtools::logratio2foldchange(DE$DE$logFC)

  #begin list of results to return
  res <- list("phenovector"=phenoVector, "DE"=DE)

  #Add Protein Annotations to the DEGs, if there are any DEGs
  if(nrow(DE$DE) == 0){
    message("No DEGs identified")
    return(res)
  }
  # will add in later
  # DE[["DE.Anno"]] <- gene_protein_anno(df=rownames_to_column(DE$DE, "gene"),
  #                                        gene.name.col = gene.name.col,
  #                                        ids2symbols = ids2symbols)

  #to avoid PCA/Heatmaps when not needed.
  if(SkipPlots){
    return(res)
  }

  # if too few genes DEGs then return without more data viz
  if (nrow(DE$DE) <= 9){
    cc <- levels(phenoVector) %>%
      magrittr::set_names(., value = c("black","firebrick"))

    eigen <- PCA(expnData, phenoVector,colorCode=cc, title=variantName)
    res[c("InputClinData", "InputExpnMatrix","PCA")] <-  list(clinData, expnData, eigen)
    return(res)
  }

    # Unsupervised Heirachrach clustering

    cols.colorbar <- c(col)
    if(!is.null(Add.Anno.Col)){
      # cols.colorbar <- c("Cytogenetic.Category.1","Cytogenetic.Category.2", "SNVs","Rare.Fusions", col) #"Cytogenetic.Category.2
      cols.colorbar <- c(cols.colorbar, Add.Anno.Col)
    }else if (!is.null(Custom.Cols)){
      cols.colorbar <- c(col, Custom.Cols)
    }

    if(!all(cols.colorbar %in% colnames(clinData))){
      message("Not all requested columns are in the provided clinData")
      return(res)
    }
      cc <- colorCodes_aheatmap(df=clinData[,cols.colorbar])

      #90th percentile absolute logFC
      mostDE <- quantile(abs(DE$DE$logFC), probs = seq(0,1,length.out = 11))[10]
      genes.to.plot <- DE$DE[order(abs(DE$DE$logFC),decreasing=TRUE), ]
      genes.to.plot <- rownames(genes.to.plot)[abs(genes.to.plot$logFC) >= mostDE]
      genes.to.plot <- genes.to.plot[seq_len(min(30, length(genes.to.plot)))]


      # dendrogram (for easier pulling out "like" groups that cluster together and clustering is done on log2 valuies that are not scaled before hand)
      dends_DE <- dge_dendrograms(expnData = expnData , #expnData can be counts or TMM normalized (but set createDGE=FALSE)
                                  pheno = phenoVector,
                                  genelist = rownames(DE$DE), #subset for only DEGs
                                  createDGE=TRUE,
                                  add.count=0.01,
                                  method=method)

      # annotation color bars
      HA <- create_HA_Labs_Hmap(expn=dends_DE$TMMCPM,
                              geneList=rownames(DE$DE),
                              CDE=clinData,
                              cols=cols.colorbar,
                              goi=genes.to.plot,
                              cc=cc) #warning: `annotation_height` is set with length of one while with multiple annotations, `annotation_height` is treated as `height`."
      # heatmap
      heatmap <- ComplexHmap(mat=dends_DE$TMMCPM,
                             hmap_anno_obj=HA$annoColumn,
                             dge_dendrograms.res=dends_DE,
                             hmap_anno_obj_genes=HA$geneLabels)



    # PCA Clustering
    cc <- c("GroupB"="grey50", "GroupA"="firebrick")
    eigen <- PCA(expnData, phenoVector, PC3=TRUE, colorCodes=cc,
                 title=variantName, GOI=GOI, ntop = 1000) # PCA on top 1000 varied genes in dataset.

    # Unconstrained Cluster Analysis/ PCoA
    genes <- rownames(DE$DE)
    MDS <- plotPCoA(assay(eigen$vst),phenoVector,geneList=genes,
                    colorCode=cc, title=variantName)

    # return the objects
    res[c( "dendrogram", "Heatmap", "MDS", "PCA")] <- list( dends_DE, heatmap, MDS, eigen)
    return(res)

}





####### Extraction methods to get items of interest. ############


#' Title
#'
#' @param twoGroups_DEGs.res output of twoGroups_DEGs()
#' @param filter boolean remove duplicate genes
#' @param goi genes of interest
#' @param anno boolean for gene annots
#' @param geneLevel boolean for gene annots or transcript annots
#'
#' @returns dataframe
#' @export
#'
#' @examples
#' ex <- c("TBD")
#'
extract_DEGs <- function(twoGroups_DEGs.res, filter=FALSE,goi=NULL,anno=FALSE, geneLevel=FALSE){
  # library(dplyr)

  #function to extract gene level annotation infromation. Added 3/22/19
  extract_anno_subset <- function(extract_DEGs.res){
    #for use after extracting all the DEGs with Annotations.

    if(is.null(dim(extract_DEGs.res))){
      return("No DEGs")
    }else{
      df <- extract_DEGs.res %>%
        dplyr::select(-Number_of_Transcripts,-Transcript_ID, -TM_Protein_Regions) %>%

        #any gene with at least 1 TMhelix detected for  any of its transcripts, has a TMhelix.
        group_by(gene) %>%
        mutate(Predicted_Transmembrane_Structure = case_when(
          any(grepl("TM", Predicted_Transmembrane_Structure)) ~ "TMhelix",
          TRUE ~ Predicted_Transmembrane_Structure)) %>%
        arrange(Cellular.Compartment_Membrane) %>%
        filter( grepl("^[A-Z].+", Cellular.Compartment_Membrane) | (!duplicated(gene, fromLast = TRUE) )) %>%
        ungroup()  %>%

        group_by(gene) %>%
        mutate_at(vars(Ensembl_ProteinID:Cellular.Compartment_Receptors),
                  ~collapseRows(., uniq = TRUE, split = TRUE, sep="; ")) %>%
        filter(!duplicated(gene)) %>%
        ungroup() %>%
        arrange(desc(logFC))

      return(df)
    }
  }

  #If statement to avoid errors whenDEGs is empty
  if(length(nrow(twoGroups_DEGs.res$DE$DE)) < 1 | ! grepl("phenovector|DE", names(twoGroups_DEGs.res))){
    return("No DEGs or object is not from twoGroups_DEGs()")
  }else{

    if(anno){
      DE <- twoGroups_DEGs.res$DE$DE.Anno

      if(geneLevel){
        DE <- extract_anno_subset(DE)
      }

    }else{
      DE <- twoGroups_DEGs.res$DE$DE %>%
        mutate(gene=rownames(.))
    }

    DE <- DE %>%
      arrange(dplyr::desc(logFC)) %>%
      dplyr::select(gene, everything())

    if (filter){
      DE <- DE %>%
        filter(gene %in% goi)
    }

    return(DE)
  }
}


# extract_MDS <- function(twoGroups_DEGs.res){
#   twoGroups_DEGs.res$MDS$plot
# }
#
#
# extract_PCA <- function(twoGroups_DEGs.res){
#   twoGroups_DEGs.res$PCA$pca_plot
# }
#
# extract_N.DE_NormFact <- function(twoGroups_DEGs.res){
#   # library(magrittr)
#   # cytogenetics <- names(twoGroups_DEGs.res)
#
#   N.DE <- nrow(twoGroups_DEGs.res$DE$DE)
#   NormFactors <- range(twoGroups_DEGs.res$DE$NormFactors$norm.factors) %>% paste(., collapse="-")
#   N.DE_NormFactor <- cbind(N.DE, NormFactors)
#   return(N.DE_NormFactor)
# }



# Annote transmembrane, subcellular location, and signal peptide genes, and ADC drug information.
#
# @param df dataframe with genes as rows. May be expression dataset or DE genes list.
# @param gene.name.col character vector of column name if any dataframe has an alternative column name for gene symbols (eg BCL2,TP53)
# @param ids2symbols an be NULL, NA, or a file path. NULL is for BCCA id mapping + backwards compatibility.
# @param mart.37 NULL or provide biomaRt object for GRCh37
# @param mart.38 NULL or provide biomaRt object for GRCh38
# @param makeQuery boolean
# @param attempts numeric
#
# @return data.frame
#
# @examples
# ex <- c('TBD')
#
# gene_protein_anno <- function(df,gene.name.col="gene",
#                               ids2symbols=NULL,
#                               mart.37=NULL,
#                               mart.38=NULL,
#                               makeQuery=TRUE,
#                               attempts=5){
#   #Modified by J.Smith to include more on the Compartments.
#
#   #df is a dataframe with genes as rows. May be expression dataset or DE genes list.
#   #ids2symbols can be NULL, NA, or a file path. NULL is for BCCA id mapping + backwards compatibility.
#   #while NA is for forwards compatibility if you don't need a ID map like with Kallisto counts.
#   #NOTE: ids2symbols must have the gene_names column be the first column!!
#   #gene.name.col is if any dataframe has an alternative column name for gene symbols (eg BCL2,TP53)
#   #gene.name.col CANNOT be gene_id.
#
#   # library(dplyr)
#   # library(stringr)
#   # library(tidyr)
#   # suppressPackageStartupMessages(library(rDGIdb))
#   # library(biomaRt)
#   # options(stringsAsFactors = FALSE)
#
#   #Function for mining certain keywords from the compartemnts data
#   matchCompartment <- function(protein, gene, geneStableID, ref.df, keywords.regex){
#
#
#     if (any(protein %in% ref.df$V1)){ #if a ENSP ID is matched, then use that information preferentially
#       comp <- ref.df %>%
#         filter(V1 %in% protein ) %>%
#         filter(grepl(keywords.regex, V4, ignore.case = TRUE)) %>%
#         #confidence score must be greater than or equal to 3 to be considered for annotation. add 3.14.19
#         filter(V7 >= 3) %>%
#         dplyr::select(V4)
#
#       res <- rep("",length(protein))
#       res[protein %in% ref.df$V1] <- paste(unlist(unique(comp)),collapse="; ")
#
#     }else{ #if not, search for the gene symbol and ENSG ID
#       comp <- ref.df %>%
#         filter( V2 %in% toupper(gene) |  V2 %in% geneStableID ) %>%
#         filter(grepl(keywords.regex, V4, ignore.case = TRUE)) %>%
#         #confidence score must be greater than or equal to 3 to be considered for annotation. add 3.14.19
#         filter(V7 >= 3) %>%
#         dplyr::select(V4)
#
#       res <-  paste(unlist(unique(comp)),collapse="; ")
#     }
#     return(res)
#   }
#
#   #Read in the external database information
#   compartment_knowledge_data <- read.delim(file.path(PROJHOME,"0000.00.02_Reference_GeneInfo/Compartments_database","human_compartment_knowledge_full_3.5.21.tsv"),
#                                            sep = "\t",as.is = TRUE, header = FALSE)
#
#   ADCs <- read.csv(file.path(PROJHOME,"0000.00.02_Reference_GeneInfo/ADC_and_CARTcell_Targets_Database_ADCReview_rmDups_clinicaltrialsGov.csv"),
#                    as.is=TRUE)
#
#   if(is.null(ids2symbols)){
#     # First, converts gene symbols present in data to gene stable IDs
#     ids2symbols <- read.csv(file.path(PROJHOME,"0000.00.02_Reference_GeneInfo","GeneSymbol_Ensembl_ID_Conversion_GRCh37.69_FromBCCA.csv"),
#                             header = TRUE)
#
#   }else if (is.na(ids2symbols)){
#     ids2symbols <- NA
#
#   }else{
#     #if want a custom gene name to ensemble gene ID file. Can contain additional annotation like ADCs or drug targets.
#     ids2symbols <- read.csv(ids2symbols,header = TRUE)
#
#   }
#
#
#   # #rename the column containing the gene symbols/names
#   # #this avoid conficlts with column names that are commonly used like gene_id, gene_name, etc
#   df <- df %>%
#     dplyr::select(geneSymbol=all_of(gene.name.col),everything())
#
#   #if df already has ENSG ids for rownames, just rename the ENSG ID column
#   if(all(is.na(ids2symbols))){
#     print(paste("No ID mapping is performed, ids2symbols is ", ids2symbols))
#
#   }else{
#     #merge in the ensembl gene IDs to the input DEGs dataframe
#     col <- colnames(ids2symbols)[1]
#     df <- df %>%
#       left_join(., ids2symbols,
#                 by=c(geneSymbol=col))
#   }
#
#   #find which column has the ensemble IDs
#   ensCol <- sapply(df, function(x) any(grep(pattern="^ENSG", x)))
#   ensCol <- names(which(ensCol))
#
#   #Create a new column called geneStableID which contains the Ensembl gene IDs
#   #this allows one to not re-name any columns in the input dataframe
#   df <- df %>%
#     mutate( geneStableID := !! as.name(ensCol))
#
#
#   #Check that the ENSG ID doesnt have NAs if using a provided gene ID map.
#   if(any(is.na(df$geneStableID) | is.null(df$geneStableID))){
#     print("NAs introduced in ID mapping. Check Reference is correct")
#     return(list("reference_used"=ids2symbols))
#   }
#
#
#   if(is.null(mart.37) & makeQuery){
#     # Uses gene stable IDs to query Ensembl for further info about the gene & associated transcripts, proteins etc.
#     # mart.92 <- mart2 <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",host = "http://apr2018.archive.ensembl.org")
#     mart.37 <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL",
#                           dataset = "hsapiens_gene_ensembl",
#                           GRCh = 37)   #GRCh37 does not have TSL (transcript support level ). Can't use mirror with this GRCh parameter
#     #if biomartr still cannot load
#     if(!exists("mart.37")){
#       mart.37 <- readRDS(file.path(PROJHOME,"0000.00.02_Reference_GeneInfo/biomaRt.GRCh37.RDS"))
#       message("Cannot load GRCh37 from BiomaRt currently. Loaded older local version.")
#     }
#   }
#
#
#   if(is.null(mart.38) & makeQuery){
#     #try to load the database. Having a ton of SSL errors, but oddly sometimes it works
#     #the errors are new the Rstudio server and R v4.0.4
#     #I need to find a way to avoid this query step. its just too buggy and unreliable...
#     # https://github.com/grimbough/biomaRt/issues/31
#     #ugh honestly fuck biomaRt. I'm going to need to just use the rest API apparently. Which will take time to figure out.
#     #I cant even get it to reliably load the goddamn mart object....
#     httr::set_config(httr::config(ssl_verifypeer = FALSE))
#     for(i in 1:attempts){
#       try(mart.38 <- useEnsembl("ensembl",
#                                 mirror = "uswest",
#                                 dataset = "hsapiens_gene_ensembl"), silent = T)
#     }
#
#     #if biomartr still cannot loaded load an older saved version
#     if(is.null(mart.38)){
#       mart.38 <- readRDS(file.path(PROJHOME,"0000.00.02_Reference_GeneInfo/biomaRt.GRCh38.RDS"))
#       message("Cannot load GRCh38 from BiomaRt currently. Loaded older local version.")
#     }
#   }
#
#   #Query Biomart for protein information.
#   attr.mart <-  c("ensembl_gene_id",
#                   "external_gene_name",
#                   "transcript_count",
#                   "ensembl_transcript_id",
#                   "ensembl_peptide_id",
#                   "tmhmm", "tmhmm_start", "tmhmm_end") #"transcript_tsl",
#
#   if(makeQuery){
#     #GRCh38 results
#     res.anno1 <- getBM(attributes = attr.mart,
#                        filters = "ensembl_gene_id",
#                        values = df$geneStableID,
#                        mart = mart.38)
#     #GRCh37 results
#     res.anno2 <- getBM(attributes = attr.mart,
#                        filters = "ensembl_gene_id",
#                        values = df$geneStableID,
#                        mart = mart.37)
#   }else{
#     #due to getBM() time-outs and useEnsembl() memory errors, the above query in real-time has become too burdensom
#     #instead use these files as the default for now on and will need to periodically update them.
#
#     #This was not finished clearly.. Need to save a stable copy for use on Rhino/Gizmo due to SSL cert issues for everything...
#     # res.anno1 <- readRDS(file.path(PROJHOME,"0000.00.02_Reference_GeneInfo/")) %>%
#     #   filter(ensemble_gene_id %in% df$geneStableID)
#     #
#     #
#     # res.anno2 <- readRDS(file.path(PROJHOME,"0000.00.02_Reference_GeneInfo/")) %>%
#     #   filter(ensemble_gene_id %in% df$geneStableID)
#   }
#
#   #gene_id's which are in GRCh37 but NOT GRCh38
#   g.idx <- which(! res.anno2$ensembl_gene_id %in% res.anno1$ensembl_gene_id)
#
#   #update results by adding GRCh37 results to the GRCh38 dataframe
#   final.res <- res.anno1 %>%
#     bind_rows(., res.anno2[g.idx,])
#
#   #protein_id's which are in GRCh37 but NOT GRCh38
#   p.idx <- which(! res.anno2$ensembl_peptide_id %in% final.res$ensembl_peptide_id)
#
#   #final update for protien ID by adding GRCh37 results to the GRCh38 dataframe (since CompartmentsDB uses ENSP IDs from GRCh37)
#   final.res <- final.res %>%
#     bind_rows(.,res.anno2[p.idx,]) %>%
#     arrange(ensembl_gene_id)
#
#   #Rename columns for clarity
#   colnames(final.res) <- c("geneStableID",
#                            "external_gene_name",
#                            "Number_of_Transcripts",
#                            "Transcript_ID",
#                            "Ensembl_ProteinID",
#                            "Predicted_Transmembrane_Structure",
#                            "Start_TM_Region", "End_TM_Region") #"Ensembl_TSL",
#
#   # groups by gene_id, then concatenates the start and stop positions of each transmembrane protein into 1 column
#   results_by_gene <- final.res %>%
#     mutate_at(vars(Start_TM_Region:End_TM_Region), ~as.character(.)) %>%
#     unite(TM_Protein_Regions, Start_TM_Region, End_TM_Region, sep = "-") %>%
#
#     #combine TM regions by protien
#     group_by(geneStableID,Ensembl_ProteinID) %>%
#     mutate_at(vars(TM_Protein_Regions), ~collapseRows(., uniq = FALSE))  %>%
#     ungroup() %>%
#
#     mutate_at(vars(TM_Protein_Regions), funs(gsub("NA-NA",NA, .))) %>%
#     mutate_at(vars(Number_of_Transcripts:TM_Protein_Regions), ~gsub("^$", NA, .)) %>%
#     unique(.)
#
#   # Adds expression data back onto the newly queried data and Annotate cellular compartments.
#   results <- df %>%
#     left_join(., results_by_gene, by = "geneStableID") %>%
#     group_by(geneSymbol) %>%
#     mutate(Cellular.Compartment_Membrane=matchCompartment(protein = Ensembl_ProteinID,
#                                                           gene = geneSymbol,
#                                                           geneStableID = geneStableID,
#                                                           ref.df=compartment_knowledge_data,
#                                                           keywords.regex = c("extracellular|plasma membrane|transmembrane|Cell periphery")),
#
#            Cellular.Compartment_Receptors=matchCompartment(protein = Ensembl_ProteinID,
#                                                            gene = geneSymbol,
#                                                            geneStableID = geneStableID,
#                                                            ref.df=compartment_knowledge_data,
#                                                            keywords.regex = c("receptor|EGFR"))) %>%
#     ungroup() %>%
#     #change the original column name back
#     dplyr::select(!! as.name(gene.name.col) := geneSymbol, everything()) #revert to original column name
#
#   #Identify small molecule inhibitors if available
#   DGI_Filter <- queryDGIdb(pull(results,gene.name.col),
#                            geneCategories = c("CLINICALLY ACTIONABLE")) #for now, only this filter since no filters has soo many drugs that are not relevant
#   DGI_Final <- detailedResults(DGI_Filter)
#
#
#   #Append ADCs and Small Molecular Inhibitors to the results
#   #https://stackoverflow.com/questions/28399065/dplyr-join-on-by-a-b-where-a-and-b-are-variables-containing-strings
#   if(nrow(DGI_Final) < 1 ){
#     print("No Interactions in Drug Gene Database")
#
#     #Merge in the ADC drugs
#     results <- results %>%
#       left_join(., ADCs,
#                 by=setNames("Gene.symbol.of.target..Final.",gene.name.col)) %>%
#       dplyr::select(gene.name.col,
#                     everything())
#
#   }else{
#
#     #Collapse multiple drugs for 1 gene-target into a single row per gene-target.
#     DGI_Final <- DGI_Final %>%
#       group_by(Gene) %>%
#       #collapse genes with multiple drugs into a single row
#       mutate_at(vars(Drug:PMIDs),
#                 ~collapseRows(col = ., uniq = FALSE, sep="; ")) %>%
#       ungroup()  %>%
#       dplyr::select(-SearchTerm) %>%
#       unique()
#
#     #Merge in the ADC and Drug-Gene interactions
#     results <- results %>%
#       left_join(., ADCs,
#                 by=setNames("Gene.symbol.of.ADC.target..Final.",gene.name.col)) %>%
#       left_join(.,DGI_Final,
#                 by=setNames("Gene",gene.name.col)) %>%
#       dplyr::select(gene.name.col,
#                     everything())
#   }
#   return(results)
# }
