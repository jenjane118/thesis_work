#function to make pca biplot with pca object from pca tools with selected pcas using ggplot

ggbiplot <- function(pca_obj, pca_list, label_var, col_var, shape_var, leg="left"){
  require(ggplot2)
  th <- theme_bw(base_size = 24) +
    theme(
      legend.background = element_rect(),
      plot.title = element_text(angle = 0, size = 16,
                                face = 'bold', vjust = 1),
      plot.subtitle = element_text(angle = 0, size = 12,
                                   face = 'plain', vjust = 1),
      plot.caption = element_text(angle = 0, size = 12,
                                  face = 'plain', vjust = 1),
      axis.text.x = element_text(angle = 0, size = 16,
                                 hjust = 0.5, vjust = 0.5),
      axis.text.y = element_text(angle = 0, size = 16,
                                 hjust = 0.5, vjust = 0.5),
      axis.title = element_text(size=16),
      legend.position = leg,
      legend.key = element_blank(),
      legend.key.size = unit(0.5, 'cm'),
      legend.text = element_text(size = 12),
      title = element_text(size = 12),
      legend.title = element_text(size = 14))
  
  plotobj <- NULL
  plotobj$x <- pca_obj$rotated[,pca_list[1]]
  plotobj$y <- pca_obj$rotated[,pca_list[2]]
  plotobj$lab <- as.character(pca_obj$metadata[[label_var]])
  plotobj <- as.data.frame(plotobj)
  plotobj$col <- pca_obj$metadata[[col_var]]
  plotobj$shape <- pca_obj$metadata[[shape_var]]
  var1 <- round(p$variance[[pca_list[1]]])
  var2 <- round(p$variance[[pca_list[2]]])
  plot <- ggplot(plotobj, aes(x = x, y = y, color = col, shape=shape, label=lab)) +
    geom_point(size = 4) +
    geom_text(hjust= 2, colour="black") +
    th +
    scale_shape_manual(name="sgRNA", labels=c("+ve", "-ve"), values=c(16, 17)) +
    xlab(paste0(pca_list[1],": ",var1,"% variance")) +
    ylab(paste0(pca_list[2], ": ", var2, "% variance")) +
    labs(color = "aTC treatment")
  plot
}

#usage: ggbiplot(p, c("PC1", "PC2"), "experiment", "treatment", "plasmid")
