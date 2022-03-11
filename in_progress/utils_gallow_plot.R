#Ivan Petrov
#https://gitlab.com/oncobox/watermelon_multisection/blob/master/utils/gallow_plot.R
## A combined dendrogram + barplot representing hierarchical structure in the data plus a single numeric value
## for each of the objects; can also display classes or clusters by colorizing the bar plot
## supports all hclust methods and custom resolution / file name (doesn't show a plot in the window though)
gallow_plot <- function(
  data_or_hc,
  num_feature,
  color_labs,
  n_clusters,
  y_bar_breaks,
  y_bar_palette,
  label_y_hist = "y",
  method = "ward.D2",
  dendro_line_size = 1.4,
  font_size_labs = 1.5,
  font_size_y_text = 10,
  font_size_y_name = 17,
  font_size_legend = 10,
  legend_key_size = 0.33,
  tick_line_size = 1,
  legend_position = "bottom",
  legend_nrow = 1,
  plot_name = "./gallow_plot.png",
  resolution = c(1600, 900)
  ) {

  # Performing hierarchical clustering
  require(magrittr)
  require(stringr)
  
  if (missing(y_bar_breaks)) {
    y_bar_breaks <- c(0:ceiling(max(num_feature)), log10(3500000))
    names(y_bar_breaks) <- c(0:ceiling(max(num_feature)), "QC")
  }

if (is.null(names(y_bar_breaks))) {
  names(y_bar_breaks) <- round(y_bar_breaks, 3)
}
  # if the data is large, we can cluster it beforehand
  if (!inherits(data_or_hc, "hclust")) {
  	hclust.data <-
  		data_or_hc %>%
    	t() %>% dist %>%
    	hclust(method = method)
  } else {
  	hclust.data <- data_or_hc
    method <- hclust.data[["method"]]
  }
  
  obj_names <- hclust.data$labels
  rearrange <- hclust.data$order
  n <- length(obj_names)
  
  require(dendroextras)
  # slicing data into clusters
  if (missing(color_labs)) {
    if (missing(n_clusters)) n_clusters <- 4
    if (is.null(n_clusters)) {
      color_labs <- rep("darkgrey", n)
    } else {
      color_labs <- dendroextras::slice(hclust.data, k = n_clusters)
    }
  }
  suppressPackageStartupMessages(require(dplyr))
  
  ### Plotting the dendrogram
  require(ggdendro)
  dendro4plot <- 
    hclust.data %>%
    dendro_data(type = "rectangle")
  
  require(ggplot2)
  dendro.ylab <- 
    ifelse(str_detect(tolower(method), "ward"),
           "Total within-cluster variance", "Between-cluster distance")
  
  data.ggdendro <-
    ggplot(segment(dendro4plot)) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend), size = dendro_line_size) +
    scale_x_continuous(limits = c(0, n + 1)) +
    scale_y_continuous(name = dendro.ylab) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = font_size_y_name),
    axis.text.y = element_text(size = font_size_y_text),
    panel.background = element_blank()
    )
  

  ### Plotting the barplot (=deviations)
  data.devbar <-
    ggplot(
      data_frame(
      	x = 1:n, 
        y = num_feature[rearrange], 
        cluster = color_labs[rearrange]
        ), 
      aes(x = x, y = y, fill = cluster)
      ) +
    geom_col(alpha = 0.75) +
    scale_x_continuous(limits = c(0, n + 1)) +
    scale_y_reverse(
      name = label_y_hist,
      breaks = y_bar_breaks,
      labels = names(y_bar_breaks)
      ) +
    geom_text(
      inherit.aes = FALSE,
      aes(
      	label = obj_names[rearrange],
        angle = 90,
        x = 1:n,
        y = -min(c(abs(num_feature), 0)),
        hjust = 1
        ),
      size = font_size_labs
      ) +
    theme(
      legend.position = "bottom",
      panel.background = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size = font_size_y_name),
      axis.text.y = element_text(size = font_size_y_text),
      panel.grid.major.y = element_line(size = tick_line_size, color = "black"),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.text = element_text(size = font_size_legend),
      legend.key.size = unit(legend_key_size, "cm")
      ) +
    guides(fill = guide_legend(title = NULL, nrow = legend_nrow))

    if(!missing(y_bar_palette)) {
      if(length(y_bar_palette) == 1) data.devbar <- data.devbar + scale_fill_brewer(palette = y_bar_palette)
      else data.devbar <- data.devbar + scale_fill_manual(values = y_bar_palette)
      }

  require(gridExtra)
  require(grid)
  gTree <- ggplotGrob(data.ggdendro)
  gBarplot <- ggplotGrob(data.devbar)
  # grid::grid.newpage()
  png(plot_name, width = resolution[1], height = resolution[2])
  # ggsave(str_c(plot_name, ".png"), arrangeGrob(gA, gB, ncol = 1, heights = c(1,1)))
  # this "rbind" variant is aligning the plots right underneath each other; it"s not exactly easy to do that by other means
  grid::grid.draw(gridExtra::gtable_rbind(gTree, gBarplot))
  dev.off()
}