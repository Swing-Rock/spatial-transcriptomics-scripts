source("scripts/st_helper.R")
library(Seurat)
library(harmony)

merge_dataset <- function(data_obj, out_path, verbose_output = FALSE) {
  #merge samples
  cat("=============== MERGING SAMPLES ===============")
  merged_obj <- merge(
    x = data_obj[[1]],              # first object
    y = data_obj[-1],               # rest of the objects
    add.cell.ids = names(data_obj), # optional prefixes for cell names
    project = "MergedProject"
  )
  
  cat("---------- running standard workflow ----------")
  merged_object <- FindVariableFeatures(merged_obj, assay = "SCT", selection.method = "vst", nfeatures = 3000)
  #Standard workflow for running PCA clustering and UMAP
  merged_object <- RunPCA(merged_object, assay = "SCT", reduction.name = "pca.SCT", verbose = verbose_output)
  merged_object <- FindNeighbors(merged_object, assay = "SCT", reduction = "pca.SCT", dims = 1:15, verbose = verbose_output)
  merged_object <- FindClusters(merged_object, cluster.name = "seurat_cluster.SCT", resolution = 0.5, verbose = verbose_output)
  merged_object <- RunUMAP(merged_object, reduction = "pca.SCT", reduction.name = "umap.SCT", return.model = T, dims = 1:15, verbose = verbose_output)
  
  cat("---------- Running Harmony ----------")
  merged_object <- RunHarmony(
    object = merged_object,
    group.by.vars = "orig.ident",      
    assay.use = "SCT",
    reduction = "pca.SCT",
    reduction.save = "harmony.SCT",
    theta = 8, 
    verbose = verbose_output
  )
  
  cat("---------- Reclustering ----------")
  merged_object <- FindNeighbors(merged_object, reduction = "harmony.SCT", dims = 1:15, verbose = verbose_output) #notice that the harmony.SCT reduction is used here
  merged_object <- FindClusters(merged_object, resolution = 0.5, cluster.name = "seurat_cluster.harmony.SCT", verbose = verbose_output)
  merged_object <- RunUMAP(merged_object, reduction = "harmony.SCT", reduction.name = "umap.harmony.SCT", dims = 1:15, verbose = verbose_output) #notice that the harmony.SCT reduction is used here. The new reduction name for the UMAP labels has been named umap.harmony.SCT
  
  
  #visualization
  p1 <- DimPlot(merged_object, reduction = "umap.SCT", group.by = "orig.ident", label = TRUE)
  p2 <- DimPlot(merged_object, reduction = "umap.harmony.SCT", group.by = "orig.ident", label = TRUE)
  t1 <- title_plot("Before Harmony")
  t1 <- title_plot("After Harmony")
    
    
  combined <- t1/ p1 / t2/ p2
    
  png(paste0(out_path, "/", "merged_cluster_plots.png"), width = 1000, height = 1000, res = 150)
  print(combined)
  dev.off()
  
  return (merged_obj)
  
}