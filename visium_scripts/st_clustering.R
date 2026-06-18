source("visium_scripts/st_helper.R")
library(Seurat)
library(patchwork)

library(SingleR)
library(celldex)

# INPUT: 
#   (Large Seurat) Seurat Obj after QC
#   (str) name of the data set 
#   (str) output path for graphical output (same as data folder by default)
# OUTPUT: SeuratObject of .h5 after qc
# make sure corresponding libraries are already downloaded
clustering <- function(data_obj, data_name, out_path, single_cell = TRUE, verbose_output = FALSE) {
  cat("===============Starting clustering===============\n")
  ref <- HumanPrimaryCellAtlasData()

  if (single_cell) {
    cat("-----running single cell labeling-----\n")
    cat("clustering\n")
    data_obj <- RunPCA(data_obj, assay = "SCT", verbose = verbose_output)
    data_obj <- RunUMAP(data_obj, reduction = "pca", dims = 1:30, verbose = verbose_output)
  } else {
    cat("-----running cluster level labeling-----\n")
    #Playing around with the dimensions of PCA/UMAP can result in different clustering algorithms.
    data_obj <- RunPCA(data_obj, assay = "SCT", verbose = verbose_output)
    data_obj <- FindNeighbors(data_obj, reduction = "pca", dims = 1:30)
    data_obj <- FindClusters(data_obj, resolution = 0.5, verbose = verbose_output) #changing the resolution changes how broadly or specifically the algorithm clusters the cell populations. A low resolution (close to 0) results in very broad clustering (fewer clusters), while high resolution results in very specific clustering (more clusters)
    data_obj <- RunUMAP(data_obj, reduction = "pca", dims = 1:30, verbose = verbose_output)
    
    #Visualization prep
    combined <- DimPlot(data_obj, reduction = "umap", label = TRUE) 
      #| SpatialDimPlot(data_obj, label = TRUE, pt.size.factor = 2.5, label.size = 3)
    
    #Table representing the number of cells in each cluster number
    table(data_obj@active.ident)
    
    # Get all cluster IDs
    clusters <- levels(data_obj)
    
    # Find markers for each cluster, store in a list
    data_obj <- PrepSCTFindMarkers(data_obj)
    
    all_markers <- lapply(clusters, function(clust) {
      FindMarkers(data_obj, ident.1 = clust, min.pct = 0.25, assay = "SCT")
    })
    
    # Name each element by cluster number
    names(all_markers) <- clusters
    cat("-----Done with clustering, starting Single R cluster identification-----\n")
  }
  
  test_data <- GetAssayData(data_obj, layer  = "data", assay = "SCT")
  
  if (single_cell){
    cat("running singleR\n")
    pred <- SingleR(
      test = test_data,
      ref = ref,
      labels = ref$label.main,
    )
    cat("adding singleR label to data_obj\n")
    data_obj$SingleR_label <- pred$labels
    
    print(sort(table(data_obj@meta.data[["SingleR_label"]]), decreasing = TRUE))
    
    #visualization
    cat("plotting\n")
    #plot combined annotated umap
    combined <- (DimPlot(data_obj, reduction = "umap", group.by = "SingleR_label", label = TRUE)) 
    png(paste0(out_path, "/", data_name, "_annotated_merged_cluster_plot.png"), width = 1500, height = 1000, res = 150)
    print(combined)
    dev.off()
    
    
    #plot annotated spatial plot for each sample
    cols <- 3
    combined <- (SpatialDimPlot(data_obj, group.by = "SingleR_label", pt.size.factor = 2.5, ncol = cols))
    n <- length(data_obj@images)
    rows <- ceiling(n / cols)
    
    png(paste0(out_path, "/", data_name, "_annotated_sample_spatial_plots.png"), width = cols*750, height = 500*rows, res = 150)
    print(combined)
    dev.off()
    
    
  } else{ 
    pred <- SingleR(
      test = test_data,
      ref = ref,
      labels = ref$label.main,
      clusters = Idents(data_obj)
    )
    cluster_labels <- pred$labels
    names(cluster_labels) <- levels(data_obj)
    data_obj <- RenameIdents(data_obj, cluster_labels)
    combined <- combined / (DimPlot(data_obj, reduction = "umap", label = TRUE) )
                            # | SpatialDimPlot(data_obj, label = TRUE, pt.size.factor = 2.5, label.size = 3))
    png(paste0(out_path, "/", data_name, "_annotated_group_cluster_plots.png"), width = 750, height = 1000, res = 150)
    print(combined)
    dev.off()
  }
  
  cat("==========Done with singleR, hope the results are good!==========\n")
  
  return (data_obj)
  
}

# ----------FUNSIES THAT MAKE VIOLIN GRAPH FOR TOP 5 GENES IN EACH CLUSTER----------
# # keep only the top 5 markers per cluster by log fold change
# top_markers <- lapply(all_markers, function(df) {
#   df <- df[order(-df$avg_log2FC), ]   # sort descending by log fold change
#   head(rownames(df), 5)                # take top 5 gene names
# })
# 
# # draw violin graph
# vln_list <- list()
# title_list <- list()
# 
# for (clust in clusters) {
#   genes <- top_markers[[clust]]
#   print(paste("Cluster", clust, "top markers:", paste(genes, collapse = ", ")))
#   vln_list[[as.character(clust)]] <- make_v_plot(obj = data_obj, feature = genes, nPerRow = 5, x_text = element_text(hjust = 1))
#   title_list[[as.character(clust)]] <- title_plot(paste0("Cluster ", clust), 18)
# }
# 
# # combine and output
# combined <- NULL
# for (clust in clusters) {
#   if (is.null(combined)) {
#     combined <- title_list[[as.character(clust)]] / vln_list[[as.character(clust)]]
#   } else {
#     combined <- combined / title_list[[as.character(clust)]] / vln_list[[as.character(clust)]]
#   }
# }
# 
#   png(paste0(out_path, "/", data_name, "cluster_violins.png"), width = 3000, height = 4000, res = 150)
# print(combined)
# dev.off()
