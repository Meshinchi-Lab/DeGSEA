#Jenny Smith 

#6/20/18 

#Purpose: Function to save multiple ggplot or base plot objects produced using dplyr do() function. 



saveMultiPlots.dplyr <- function(dplyr.do,w=8,h=5){
  #This is for the relapse results from plots with do() command in 
  
  N <- nrow(dplyr.do)
  chars <- sapply(dplyr.do,class) == "character"
  
  if (sum(chars) == 1){
    name <- function(i){paste(names(dplyr.do)[1], dplyr.do[[1]][i],col,".tiff", sep="_")}
  }else{
    name <- function(i){paste( paste(names(dplyr.do)[chars], collapse = "_"),  paste(dplyr.do[chars][i,], collapse = "_"), col,".tiff", sep="_")}
  }
  
  #For comp risk regression 
  cols <- colnames(dplyr.do)[!chars] %>% .[. != "compRiskReg"]
  
  for (col in cols){
    # print(col)
    if(grepl("ggplot", class(dplyr.do[[col]][[1]]))){
        lapply(1:N, function(x) ggsave(filename = name(x),
                                       plot = dplyr.do[[col]][[x]],
                                       device = "tiff",
                                       width = w,
                                       height = h,
                                       dpi=600))
    }else{
      for (i in 1:N){
      # print(name(i))
      tiff(name(i), height = 5, width = 6, units="in", res=600)
      plot(dplyr.do[[col]][[i]], color = c(1:8), lwd=3, xlab = "Days")
      dev.off()
      }
    }
  }
}