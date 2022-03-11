#Jenny Smith

#May 15, 2017

#Themes for ggplot that are are customized. 

library(pryr)


#Function to save multple plots when created with dplyrs do() function. 
saveMultiPlots.dplyr <- function(dplyr.do,w=8,h=5){
  #This is for the relapse results from plots with do() command in 
  
  N <- nrow(dplyr.do)
  
  name <- function(i){paste(dplyr.do[[1]][i],col,".tiff", sep="_")}
  cols <- colnames(dplyr.do)[-1]
  
  for (col in cols){
    lapply(1:N, function(x) ggsave(filename = name(x),
                                   plot = dplyr.do[[col]][[x]],
                                   device = "tiff",
                                   width = w,
                                   height = h,
                                   dpi=600))
    
  }
}



theme_JS %<a-% { theme(plot.title = element_text(hjust = 0.5, size = 18),
                       panel.background = element_rect(fill="white"),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.border = element_rect(color = "black", fill=NA),
                       axis.text = element_text(color = "black"),
                       axis.text.x = element_text(angle = 20,hjust=1,vjust = 1, size = 18),
                       axis.text.y = element_text(size = 22),
                       axis.title = element_text(size = 20), 
                       legend.key=element_blank()) 
}

theme_rotateLabs %<a-% { theme(plot.title = element_text(hjust = 0.5, size = 18),
                       panel.background = element_rect(fill="white"),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.border = element_rect(color = "black", fill=NA),
                       axis.text = element_text(color = "black"),
                       axis.text.x = element_text(angle = 45,hjust=1,vjust = 1, size = 15),
                       axis.text.y = element_text(size = 18),
                       axis.title = element_text(size = 18), 
                       legend.text = element_text(size=15),
                       legend.key=element_blank())
}

theme_numX %<a-% { theme(plot.title = element_text(hjust = 0.5, size = 18),
                       panel.background = element_rect(fill="white"),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.border = element_rect(color = "black", fill=NA),
                       axis.text = element_text(color = "black"),
                       axis.text.x = element_text(angle = 0,hjust=0.5,vjust = 0.5, size = 16),
                       axis.text.y = element_text(size = 16),
                       axis.title = element_text(size = 18), 
                       legend.key=element_blank())
}




