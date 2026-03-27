merged_obj <- merge(
  x = data1, 
  y = list(data2, data3),  # list of other objects to merge
  add.cell.ids = c("Sample1", "Sample2", "Sample3"),  # optional, prefixes for cell names
  project = "MergedProject"  # optional project name
)

merged_object <- FindVariableFeatures(merged_obj, assay = "SCT", selection.method = "vst", nfeatures = 3000)



#Standard workflow for running PCA clustering and UMAP
merged_object <- RunPCA(merged_object, assay = "SCT", reduction.name = "pca.SCT")
merged_object <- FindNeighbors(merged_object, assay = "SCT", reduction = "pca.SCT", dims = 1:15)
merged_object <- FindClusters(merged_object, cluster.name = "seurat_cluster.SCT", resolution = 0.5)
merged_object <- RunUMAP(merged_object, reduction = "pca.SCT", reduction.name = "umap.SCT", return.model = T, dims = 1:15)

DimPlot(merged_object, reduction = "umap.SCT", group.by = "orig.ident", label = TRUE)

merged_object <- RunHarmony(
  object = merged_object,
  group.by.vars = "orig.ident",      
  assay.use = "SCT",
  reduction = "pca.SCT",
  reduction.save = "harmony.SCT",
  theta = 8
)
merged_object <- FindNeighbors(merged_object, reduction = "harmony.SCT", dims = 1:15) #notice that the harmony.SCT reduction is used here
merged_object <- FindClusters(merged_object, resolution = 0.5, cluster.name = "seurat_cluster.harmony.SCT")
merged_object <- RunUMAP(merged_object, reduction = "harmony.SCT", reduction.name = "umap.harmony.SCT", dims = 1:15) #notice that the harmony.SCT reduction is used here. The new reduction name for the UMAP labels has been named umap.harmony.SCT

DimPlot(merged_object, reduction = "harmony.SCT", group.by = "orig.ident", label = TRUE)