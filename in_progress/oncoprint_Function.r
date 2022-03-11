#Jenny Smith 
#10/26/18


#Purpose: Create an Oncoprint with geom_tile(). Of note, this is less of a true function bc so much is hard coded. The sample ordering has been the absolute most difficult. 



oncoprint_plot <- function(CDE, Colunms=NULL,Regex=NULL, subset=FALSE,subset.col=NULL, subset.group="Yes"){
  library(dplyr)
  library(stringr)
  library()
  #CDE is the clinical data frame
  #Columns are the charatcter vector of column names to include. 
  #Regex would be the pattern for matches() function. 
  
  #If you want to filter for a particular subgroup, like normal karyotype or only NUP98-KDM5A positives. 
  if(subset){
    reg <- paste0("^", subset.group, "$")
    
    CDE <- 
      filter(grepl(reg, subset.col))
    
  }
  
  #Select columns for Oncoprint
  if(all(!is.null(Columns) &  !is.null(Regex))){
    tile.df <- CDE %>% 
      select(contains("USI",ignore.case = F),Columns, matches(Regex)) 
      
  }else if(!is.null(Columns)){
    tile.df <- CDE %>% 
      select(contains("USI",ignore.case = F), Columns) 
    
  }else if(!is.null(Regex)){
    tile.df <- CDE %>% 
      select(contains("USI",ignore.case = F), matches(Regex)) 
    
  }
  
  temp <- tile.df %>%
    mutate_all(funs(as.character(.))) %>% #remove any factor levels that might already be in the original data
    
    #if the primary cyto code is used. 
    mutate(Primary.Cytogenetic.Code= if(exists("Primary.Cytogenetic.Code", where = .)) gsub("Unknown", "Other",
                                         gsub("\\(|\\)", "\\.", Primary.Cytogenetic.Code)) else NA) %>%
    
    #change Unknowns to No/Negative for plotting.   
    mutate_all(funs(ifelse(grepl("Unknown",.), "No", .))) %>%
    
    mutate_if(.predicate = all(grepl("Yes|No", ., ignore.case = TRUE)), .funs = funs(factor(. , levels=c("No","Yes")))) %>%
    mutate_at(vars(matches("CBFA|tri|NP|CE|FLT|M7|M6")),
              funs(Num.Status=as.numeric(.))) %>%
    
    mutate(Num.CytoCode=as.numeric(factor(Primary.Cytogenetic.Code,levels=c("Normal","inv.16.","Other")))+ 4) %>%
    mutate(Num.AgeCat=as.numeric(factor(Age.Category, levels=c("Between 3 and 5 years","Less than 3 years" ))) + 7) %>%
    gather(var,val,matches("^Num|_Num"))
  
  
  
  tile.w.perc <- tile.df %>% 
    select(TARGET.USI.1,var,val) %>%
    
    group_by(var) %>%
    mutate(Total=sum(val)) %>%
    ungroup() %>%
    
    arrange(desc(Total))
  
  CytoOrder <- unique(tile.w.perc$var)
  
  
  CytoOrder
  
  Sample.Order <- tile.w.perc %>%
    filter(!grepl("Primary|Age|CBFA2T3",var)) %>% #these variable don't matter for thier order
    mutate(var=factor(var,levels=CytoOrder)) %>%
    arrange(var) %>%
    
    group_by(TARGET.USI.1) %>%
    mutate(Total_PerPatient=sum(val)) %>%
    ungroup() %>%
    
    filter(!duplicated(TARGET.USI.1)) %>% 
    inner_join(., unique(select(tile.df, TARGET.USI.1, Primary.Cytogenetic.Code, Age.Category,M7_AML)), by="TARGET.USI.1") %>%
    mutate(Primary.Cytogenetic.Code=factor(Primary.Cytogenetic.Code, levels=c("Other","Normal","inv.16."))) %>%
    arrange(Primary.Cytogenetic.Code,desc(Total_PerPatient),desc(M7_AML), Age.Category)
  
  
  Sample.Order
  dim(Sample.Order)
  
  tile.order <- tile.df %>%
    mutate(var=factor(var, levels = rev(CytoOrder))) %>%
    mutate(TARGET.USI.1=factor(TARGET.USI.1, levels = Sample.Order$TARGET.USI.1))
  
  
  factor.labs <- c("Negative","Positive","Normal Karyotype","Inv.16", "Other", "Between 3 and 5 years", "Less than 3 years")
  factor.cols <- c("1"="navy",
                   "2"="red",
                   "5"="#FFFF33",
                   "6"="#984EA3",
                   "7"="#A65628",
                   
                   "8"="#F781BF",
                   "9"="#999999")
  
  
  labs <- gsub("Num\\.|_Num.+|\\.positive.","",levels(tile.order$var)) %>%
    gsub("_|\\.", " ", .) %>%
    ifelse(grepl("AgeCat", .), "Age Category", .) %>%
    ifelse(grepl("Cyto", .), "Primary Cytogenetic Code",.) %>%
    ifelse(grepl("CBFA2T3", .), "CBFA2T3.GLIS2", .)
  
  oncotile <- ggplot(tile.order, aes(x=TARGET.USI.1, y=var, fill = factor(val))) + 
    geom_tile(color="white", size=1.2) +
    scale_y_discrete(labels=labs) +
    scale_fill_manual(values = factor.cols,
                      labels =  factor.labs) +
    labs(title="Co-occuring Cytogenetic Abnormalities in CBFA2T3-GLIS2", y="", x="") +
    theme(axis.text.x = element_text(angle=45,hjust = 1, vjust=1, color="black"),
          axis.ticks.y = element_blank(),
          # axis.ticks.x = element_blank(),
          axis.text.y = element_text(size=16, color="black"),
          plot.title = element_text(size=18, hjust=0.5), 
          legend.title = element_blank(),
          legend.text = element_text(size=16))
  
  
}



