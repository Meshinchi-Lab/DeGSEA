#Jenny Smith
#3/28/18


color_bars <- function(list.labels,colorDends_Groups.res){
  
  cb1 <- list.labels[[1]][names(colorDends_Groups.res$groups)] #subset & order
  cb2 <- list.labels[[2]][names(colorDends_Groups.res$groups)]
  
  cb.all <- data.frame(Cytogenetics=cb1,
                       FLT3.ITD=cb2,
                       Cytogenetics.Num=as.numeric(as.factor(cb1)),
                       FLT3.ITD.Num=as.numeric(as.factor(cb2)),
                       Hier.Cluster.Group=colorDends_Groups.res$groups)
  
  c <- c("red","deepskyblue","darkorchid1","blue3","gold2")
  c2 <- c("darkolivegreen2","darkslategray")
  c3 <- c("firebrick1", "dodgerblue3")
  
  colors <- data.frame(Fusions=c[cb.all[,"Cytogenetics.Num"]],
                       FLT3=c2[cb.all[,"FLT3.ITD.Num"]],
                       Group=c3[cb.all[,"Hier.Cluster.Group"]],
                       stringsAsFactors = TRUE)
  return(colors)
}


#CDE

merged <- read.csv("~/reference_mapping-files/TARGET_AML_1031_0531_Merged_CDE_3.30.18.csv",
                   stringsAsFactors = FALSE)

#Color Bars

cols <- c("CBFA2T3.GLIS2","NUP98.KDM5A","RBM15.MKL1","DEK.NUP214")
cols2 <- c("CEBPA.mutation","NPM.mutation")
labels1 <- pheno_bars(CDE=merged, IDCol = "TARGET.USI.1", cols=cols)
labels2 <- pheno_bars(CDE=merged, IDCol = "TARGET.USI.1", cols="FLT3.ITD.positive.")
labels3 <- pheno_bars(CDE=merged, IDCol="TARGET.USI.1",cols=cols2)
color.df <- color_bars(list(labels1,labels2), CBF.GLIS.like)


#Genes of Interest (GOI)
topDEGs <- CBFvsOtherAML.1031 %>%
  filter(logFC < -3.488684  | logFC >  7.541174 ) %>%
  select(gene) %>%
  unlist()


#Cluster on GOI
d.top <- dge_dendrograms(expnData=CBFvsAML$cts.hd.1031$InputExpnMatrix, 
                         pheno=CBFvsAML$cts.hd.1031$phenovector,
                         genelist = topDEGs,
                         method = "ward.D2")


#Color dengrogram labels and branches

col.top <- ifelse(CBFvsAML$cts.hd.1031$phenovector == "GroupA", "Red", "dark blue")
cc.top <- c(GroupA="red",GroupB="dark blue")


CBF.GLIS.like <- colorDends_Groups(dge_dendrograms.res = d.top, 
                                   phenovector = CBFvsAML$cts.hd.1031$phenovector,
                                   k=2,
                                   branchCol = c("firebrick","navy"),
                                   colorcodes = cc.top)

plot(CBF.GLIS.like$split_dends[[1]])
plot(CBF.GLIS.like$split_dends[[2]])


#Add Color bar below color dendrogram
par(mfrow=c(1,1), cex=0.125, mar=c(35, 7.5, 8.5, 2), pty="m")
plot(CBF.GLIS.like$dend, axes=TRUE,cex.axis=9, horiz=FALSE)
par(cex=0.8, cex.main = 1, cex.lab = 0.85)
colored_bars(colors = color.df, y_scale=80, rowLabels=c("", ""))





