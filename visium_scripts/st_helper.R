library(Seurat)
library(ggplot2)
library(hdf5r)
library(Matrix)

# Load a Visium sample from an h5ad file + spatial/ folder.
# Returns a SeuratObject with a "Spatial" assay and spatial image attached,
# equivalent to what Load10X_Spatial() produces from an .h5 file.
load_h5ad_visium <- function(data_dir) {
  h5ad_files <- list.files(data_dir, pattern = "\\.h5ad$", full.names = TRUE)
  if (length(h5ad_files) == 0) stop("No .h5ad file found in: ", data_dir)
  if (length(h5ad_files) > 1) warning("Multiple .h5ad files found; using the first: ", h5ad_files[1])
  h5ad_path <- h5ad_files[1]
  spatial_dir <- file.path(data_dir, "spatial")
  
  h5 <- H5File$new(h5ad_path, mode = "r")
  on.exit(h5$close_all())
  
  # Gene names and spot barcodes
  gene_names <- h5[["var"]][["_index"]][]
  spot_barcodes <- h5[["obs"]][["_index"]][]
  
  # Count matrix stored as CSC sparse matrix in h5ad
  x_data    <- h5[["X"]][["data"]][]
  x_indices <- h5[["X"]][["indices"]][]   # 0-based row indices
  x_indptr  <- h5[["X"]][["indptr"]][]    # column pointers
  
  counts <- sparseMatrix(
    i    = x_indices + 1L,          # convert 0-based → 1-based
    p    = x_indptr,
    x    = x_data,
    dims = c(length(gene_names), length(spot_barcodes)),
    dimnames = list(gene_names, spot_barcodes)
  )
  
  seu <- CreateSeuratObject(counts = counts, assay = "Spatial")
  
  # Attach the Visium spatial image from the spatial/ folder
  image <- Read10X_Image(image.dir = spatial_dir, filter.matrix = TRUE)
  image <- image[Cells(seu)]
  DefaultAssay(object = image) <- "Spatial"
  seu[["slice1"]] <- image
  
  return(seu)
}

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

# Spatial image feature plot for Xenium (ImageFeaturePlot)
make_if_plot <- function(xenium_obj, features, ncol = 3) {
  ImageFeaturePlot(xenium_obj, features = features, ncol = ncol) &
    theme(legend.position = "right")
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
        margin = ggplot2::margin(t = top_margin, b = bottom_margin)
      )
    )
}
