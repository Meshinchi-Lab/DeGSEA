#Jenny Smith
#Jan 18, 2019





# tiff("DnGenes_UpMirs_DotPlot.tiff", height = 18.5, width = 22, units="in", res=300)
ggplot(apop.ints.p, aes(x=miR.labs, y= Gene_symbol, fill=Correlation, size=Sum)) +
    geom_point(shape=21, stroke=0.1) + #aes(stroke=Validate)
     geom_point(data = filter(apop.ints.p, Validate==TRUE),
               aes(x=miR.labs, y=Gene_symbol, color=Validate),
               shape=21,stroke=1.65) +

  labs(x="", y="", title="miRNA-mRNA Interactions: Up-Regulated miRNAs vs Down-Regulated Genes", color="Validated", size="Number of Predictions") +
  theme(plot.margin = margin(t = 3, unit = "mm"),
        plot.title = element_text(hjust = 0.5, size = 30),
                       panel.background = element_rect(fill="white"),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.border = element_rect(color = "black", fill=NA),
                       axis.text = element_text(color = "black"),
                       axis.text.x = element_text(angle = 35,hjust=1,vjust = 1, size = 22, face="bold"),
                       axis.text.y = element_text(size = 22, face="bold"),
                       axis.title = element_text(size = 10),
                       legend.key=element_blank(),
                        legend.text = element_text(size=18)) +
  scale_fill_gradientn(colors=c("red4","red","orange"), values = c(0,0.75,1)) +
  scale_size_area(max_size = 16.5) +
  scale_color_manual(values = c("TRUE"="blue"))

# dev.off()
