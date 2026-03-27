library(Seurat)
library(ggplot2)

#draws violin plots
make_v_plot <- function (obj, feature, nPerRow, x_text = element_blank()){
  return (VlnPlot(obj, features = feature,
                  pt.size = 0.1, ncol = nPerRow) +
            theme(axis.title.x = element_blank(),
                  axis.text.x = x_text,
                  axis.ticks.x = element_blank()))
}

#draws spatial feature plot
make_sf_plot <- function (SeruratObj, feature, nPerRow){
  return(SpatialFeaturePlot(SeruratObj, features = feature, ncol = nPerRow) &
           theme(legend.position = "right"))
}

#make&format title labels
title_plot <- function(txt, padding = 22, top_margin = 0, bottom_margin = 5){
  ggplot() + 
    theme_void() + 
    ggtitle(txt) +
    theme(
      plot.title = element_text(
        hjust = 0.5,             # center
        vjust = -padding,               
        size = 20,
        face = "bold",
        margin = margin(t = top_margin, b = bottom_margin)
      )
    )
}