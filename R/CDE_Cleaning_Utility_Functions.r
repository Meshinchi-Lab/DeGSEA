#Jenny Smith
#12/12/18
#helper functions for cleaning the CDEs and using the fusion columns


#' For each duplicate USI, collapse the columns' information so that it will be seperated by a semi-colon.
#'
#' @param df the data.frame
#' @param ID.column the name of the column with sample IDs
#' @param duplicate boolean if sample ID is duplicated or not
#'
#' @return data.frame
#' @export
#'
#' @examples
#' my_data <- data.frame(Sample=paste0("s",1:20), mutation=rep(c("Yes","No"), length.out=20))
#' collapseDuplicates(my_data, ID.column = "Sample", duplicate=c("s2","s5"))
#'
collapseDuplicates <- function(df,ID.column,duplicate){
  #Purpose: for each duplicate USI, collapse the columns' information so that it will be seperated by a semi-colon.
  #This retains all information for each sample, but can remove duplicate USIs.

  #df is the datframe with multiple patient entries to collapse
  #ID column is the column to match the USI
  #duplicate is the USI of the dups.

  #Use example:
  # dups <- unique(df$Patient.ID[duplicated(df$Patient.ID)])
  #collapsed.subset <-  bind_rows(lapply(dups, function(x) collapseDuplicates(df=df, ID.column="Patient.ID", duplicate=x))
  # rmDups <- df[!(df$Patient.ID %in% dups), ] %>%
  # bind_rows(. , collapsed.subset)


  idx <- which(df[, ID.column] == duplicate)
  cde <- df[idx,]

  if (length(unique(cde)) == 1){
    cde <- unique(cde)
  }else{
    #Examine eac column seperately
    for (i in 1:ncol(cde)){
      #if all identical, just unique it
      if (length(unique(cde[,i])) == 1){
        cde[1,i] <- unique(cde[,i])
      }else{
        #otherwise, collap
        cde[1,i] <- paste(unique(cde[,i]), collapse = ";")
      }
    }
  }

  #update the clinical annotations with only the merged cde.
  cde <- cde[1,]


  return(cde)
}


collapseRows <- function(col, uniq=FALSE){
  #designed for dplyr so that "col" paramter is a vector of that column.
  #Similar to collapseDuplicates(), but for fewer columns, plus preselection of columns. Where are collapseDuplicates() you don't need to know the exact column names before hand.
  if (uniq){col <- unique(col)}

  collapsed <- ifelse(all(is.na(col)), NA, paste(col, collapse = "; "))
  return(collapsed)
}


GroupIDs <- function(clinData, col){
  #clindata has patient IDs as rownames.
  #col is a chracter string of the column name in the clinical data with the factor/variable information.
  list <- list()
  grps <- unique(clinData[,col])
  N <- length(grps)

  for (i in 1:length(grps)){
    # if (grepl("[^a-zA-Z0-9 :]", grps[i])){
    #   grps[i] <- gsub("[^a-zA-Z0-9 :]", "\\.", grps[i]) #remove special characters and replace with a "."
    # }
    group=grps[i]
    Column=clinData[,col]
    IDs <- rownames(clinData[grepl(group,Column), ])
    list[[i]] <- IDs
  }

  names(list) <- grps
  return(list)
}


phenoVectors_MultipleGroups <- function(listofgoupsIDs){
  # library(magrittr)
  #listofgoupsIDs contains a list with each item containing the IDs for each factor level.
  #See GroupIDs function - this produced the input called  "listofgroupIDs"
  group <- names(listofgoupsIDs)

  vector <- NULL
  names <- NULL
  for (i in 1:length(listofgoupsIDs)){
    g <- group[i]
    vector <- c(vector, rep(g, length(listofgoupsIDs[[i]])))
    names <- c(names, listofgoupsIDs[[i]])
  }

  names(vector) <- names
  return(vector)
}

phenoVectors <- function(groupA, groupB){
  # library(magrittr)
  #groupA and GroupB are character vectors with the patients IDs in each group
  g1 <- as.character(substitute(groupA))
  g2 <- as.character(substitute(groupB))

  vector <- c(rep(g1, length(groupA)), rep(g2, length(groupB)))
  names(vector) <- c(groupA, groupB)

  return(vector)
}



pheno_bars <- function(CDE,IDCol,cols){
  #CDE is the clinical data frame with patietns as rows.
  #IDcol is the name of the column with patient USIs or COG#s
  #cols are the colnames to be combined.

  replace_yes <- function(col,name){
    name <-gsub(".RNASeqCalls|.positive.", "", name)
    col <- ifelse(grepl("Yes", col, ignore.case = TRUE), name, col)
    return(col)
  }

  phenobar.df <- CDE %>%
    select(IDCol,cols)

  if(length(cols) > 1){
    phenobar.df <- bind_cols(phenobar.df, mapply(replace_yes, CDE[,cols], cols, SIMPLIFY = FALSE))
  }else{
    new <- data.frame(replace_yes(CDE[[cols]],cols)) %>% set_colnames(cols)
    phenobar.df <- bind_cols(phenobar.df, new) #dplyr bind_cols throws error Error in cbind_all(x) : Argument 2 must have names??
  }


  p <- NULL
  for (col in cols){p <- paste(p,phenobar.df[[paste0(col,1)]], sep="_")}


  phenobar <- p %>%
    gsub("No|Unknown|_", "", .) %>%
    gsub("^$", "OtherAML",.) %>%
    set_names(CDE[[IDCol]])


  return(phenobar)

}

#Create function to remove duplicate fusions
removeDups <- function(col.rowwise){
  #use with rowwise() in dplyr
  str_split(col.rowwise, "; ") %>%
    unlist() %>%
    gsub(" ","",.) %>%
    unique(.) %>%
    .[. != ""] %>%
    paste(.,collapse = "; ")
}


#Create function to remove duplicate fusions
removeDups <- function(col.rowwise){
  #use with rowwise() in dplyr
  str_split(col.rowwise, "; ") %>%
    unlist() %>%
    gsub(" ","",.) %>%
    unique(.) %>%
    .[. != ""] %>%
    paste(.,collapse = "; ")
}


#Function to split the fusion columns into multiple single columns, coded as Yes/No/Unknown
createMultiCols <- function(col,split=FALSE,suffix){
  #col is the column with many factors (eg all fusions). character vector lenght 1.
  #suffix is a character string for the suffix on all column names
  #desinged for dplyr bind_cols() function

  #example
  # df <- df %>%
  #   bind_cols(.,createMultiCols(.$Fusion, suffix=".sensefusion"))

  if(split){
    groups <- unique(gsub(" ","",unlist(str_split(col, "; "))))
    groups <- groups[groups != ""]
  }else{
    groups <- unique(col)
  }
  list <- lapply(groups, function(x) ifelse(grepl(paste0(x, ";"), col) | grepl(paste0(x, "$"), col) , "Yes",
                                            ifelse(grepl("^$|Unknown",col) | is.na(col), "Unknown","No")))
  list.names <- gsub("-", ".",  groups)
  names(list) <- paste0(list.names,suffix)
  return(list)
}


splitFusions <- function(FusionCol,regex){
  #FusionCol is the column name containing the semi-colon seperated values of fusions, patients as rows
  #regex is the fusion of interest
  #use with rowwise in dplyr

  # library(stringr)

  fusions <- unlist(str_split(FusionCol, "; |, "))
  fus <- unique(grep(regex, fusions, value = TRUE, ignore.case = TRUE))



  if (length(fus) < 1){
    fus <- ""
  }else if (length(fus) == 1){
    fus <- fus
  }else if (length(fus) > 1){

    #Order by alphabetical.
    fus <- str_split(fus, "-",simplify = FALSE) %>%
      lapply(., function(x) x[order(x)])

    #Check for reciprocal fusions
    test.reciprocal <- sapply(1:length(fus), function(x) identical(x=x,y=fus[1]))

    if(all(test.reciprocal)){
      #if identical keep the first one. Reassemble to GeneA-GeneB notation.
      fus <- paste(fus[[1]],collapse = "-")

    }else{
      #if there are more than 2 fusions, will check on later to make function more generalizable. Now, just collapse them.
      fus <- sapply(fus, function(y) paste(y, collapse = "-")) %>%
        paste(., collapse = "; ")
    }
  }

  # fus <- paste(fus[grep(regex,fus)], fus[-grep(regex,fus)], sep="-")

  # #second check
  # if (length(fus) > 1){
  #   fus <- paste(fus, collapse=";")

  return(fus)
}



#Function to combine reciprocol fusions.
combineFusions <- function(df, regex, identifier){
  cols <- grep(regex, colnames(df), value=TRUE)
  new <- NULL

  for (col in cols){
    new <- paste(new,unlist(df[,col]), sep="; ")
  }

  new <- unlist(lapply(new, removeDups))
  new <- gsub("Yes; No|No; Yes", "Yes", new) #change any calls to "Yes" if either reciprocol fusion is defined as "Yes"

  if (length(new) > 0){
    names(new) <- identifier
  }

  return(new)
}


# pheno_bars <- function(CDE,IDCol,cols){
#   #CDE is the clinical data frame with patietns as rows.
#   #IDcol is the name of the column with patient USIs or COG#s
#   #cols are the colnames to be combined.
#
#   #NOTE: all the unlist() functions are because tibbles do not play nice with base R...
#
#   replace_yes <- function(col,name){
#     name <-gsub(".RNASeqCalls|.positive.", "", name)
#     col <- ifelse(grepl("Yes", col, ignore.case = TRUE), name, col)
#     return(col)
#   }
#
#   phenobar.df <- CDE %>%
#     select(IDCol,cols)
#
#   if(length(cols) > 1){
#     phenobar.df <- bind_cols(phenobar.df, mapply(replace_yes, CDE[,cols], cols, SIMPLIFY = FALSE))
#   }else{
    # new <- data.frame(replace_yes(unlist(CDE[,cols]),cols)) %>% set_colnames(cols)
#     phenobar.df <- bind_cols(phenobar.df, new) #dplyr bind_cols throws error Error in cbind_all(x) : Argument 2 must have names??
#   }
#
#
#   p <- NULL
#   for (col in cols){p <- paste(p,unlist(phenobar.df[,paste0(col,1)]), sep="_")}
#
#
#
#   phenobar <- p %>%
#     gsub("No|Unknown|_", "", .) %>%
#     gsub("^$", "OtherAML",.) %>%
#     set_names(unlist(CDE[,IDCol]))
#
#
#   return(phenobar)
#
# }
