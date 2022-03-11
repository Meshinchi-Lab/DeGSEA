#Jenny Smith 

#December 2017 

#Purpose: Create a barcode plot using either SNR (or is Kolmogorov-Smirnov statistic?) statistics from Broad GSEA, or limma moderated t-stats. 
# Will check on that. I know genes are ranked by Signal-to-noise (SNR) and KS stat is used to 
#determine the significance of genes "inside" the geneset vs genes "outside" the geneset. So that is just the nominal P-value. 
#or right, so using the rank stat is SNR. 


#Create Barcode Plots of Enrichement


GSEA.res <- read.csv("GSEA/only_1031/GSEA_50TPM.Filter_1000xPhenoPerm_C2.Kegg.Gsea.1510966589373/gsea_report_for_pos.1031_1510966589373.csv", stringsAsFactors = FALSE)

HH.res <- read.csv("GSEA/only_1031/GSEA_50TPM.Filter_1000xPhenoPerm_C2.Kegg.Gsea.1510966589373/KEGG_HEDGEHOG_SIGNALING_PATHWAY.csv", stringsAsFactors = FALSE)

Rank.Stat <- read.csv("GSEA/only_1031/GSEA_50TPM.Filter_1000xPhenoPerm_C2.Kegg.Gsea.1510966589373/ranked_gene_list_pos.1031_versus_neg.1031_1510966589373.csv", stringsAsFactors = FALSE)

head(Rank.Stat)




# library(limma)
idx <- Rank.Stat$NAME %in% HH.res$PROBE #49 genes

# SigtoNoise.stat <- HH.res$RANK.METRIC.SCORE %>%
#   set_names(HH.res$PROBE)

SigtoNoise.stat <- Rank.Stat$SCORE %>%
  set_names(Rank.Stat$NAME)


# idx <- HH.res$CORE.ENRICHMENT == "Yes"
# HH.res$RANK.METRIC.SCORE[idx] <- HH.res$RANK.METRIC.SCORE[idx] + 1

# tiff("GBCGLIS_vs_OtherAML_C2.HH.Kegg_Barcodeplot.tiff", height= 5, width = 7, units="in", res=600)
barcodeplot(SigtoNoise.stat, index=idx)
# dev.off()




t.stat <- CBFGLISvsOtherAML$DE$eBayesFit$t[,1]
idx <- names(t.stat) %in% HH.res$PROBE #36
# tiff("GBCGLIS_vs_OtherAML_C2.HH.Kegg_Limma_tstat_Barcodeplot.tiff", height= 5, width = 7, units="in", res=600)
barcodeplot(t.stat, index = idx)
# dev.off()



head(GSEA.res)



sel.gsea <- GSEA.res %>%
  slice(1:5) %>%
  select(NAME,FDR.q.val) %>%  
  mutate(FDR.q.val.fix = ifelse(FDR.q.val == 0, 1e-3, FDR.q.val)) %>%
  mutate(negLog2FDR = -log2(FDR.q.val.fix)) %>%
  mutate(path=gsub("_", "\n", NAME))

# tiff("GBCGLIS_vs_OtherAML_C2.Kegg_FDR_barplot_addExtraRoom.tiff", height= 7, width = 7, units="in", res=600)
ggplot(sel.gsea, aes(y=negLog2FDR,x=reorder(sel.gsea$path, -sel.gsea$FDR.q.val), fill=NAME)) +
  geom_bar(stat="identity") +
  scale_y_continuous(limits = c(0,14),breaks = seq(0,14,by=2)) +
  scale_fill_brewer(palette="Dark2") +
  geom_hline(aes(yintercept = -log2(0.05)), linetype="dashed") +
  theme_numX + labs(x="", y="-log2(FDR)") + guides(fill=FALSE) +
  coord_flip()
# dev.off()


