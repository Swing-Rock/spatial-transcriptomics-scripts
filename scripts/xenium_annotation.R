source("visium_scripts/st_helper.R")
library(Seurat)
library(patchwork)
library(SingleR)
library(celldex)

# Cluster and annotate a single Xenium sample.
#
# Because Xenium is single-cell resolution, we annotate at the cell level
# (single_cell = TRUE) or at the cluster level (single_cell = FALSE),
# mirroring the behavior in st_clustering.R.
#
# INPUT:
#   (Seurat) data_obj      — post-QC / SCTransform Xenium Seurat object
#   (str)    data_name     — sample identifier for output files
#   (str)    out_path      — directory for PNG output
#   (bool)   single_cell   — TRUE: label each cell individually;
#                            FALSE: cluster first, then label each cluster
# OUTPUT: Seurat object with SingleR labels in metadata / active ident

run_xenium_annotation <- function(data_obj, data_name, out_path,
                                  single_cell = TRUE, verbose_output = FALSE) {

  cat("=============== Starting annotation for:", data_name, "===============\n")

  ref <- HumanPrimaryCellAtlasData()

  # ── Dimensionality reduction ───────────────────────────────────────────────
  cat("Running PCA...\n")
  data_obj <- RunPCA(data_obj, assay = "SCT", verbose = verbose_output)

  if (!single_cell) {
    cat("--- Cluster-level labeling: finding clusters first ---\n")
    data_obj <- FindNeighbors(data_obj, reduction = "pca", dims = 1:30, verbose = verbose_output)
    data_obj <- FindClusters(data_obj, resolution = 0.5, verbose = verbose_output)
    cat("Cluster sizes:\n")
    print(table(data_obj@active.ident))
  }

  data_obj <- RunUMAP(data_obj, reduction = "pca", dims = 1:30, verbose = verbose_output)

  # ── SingleR annotation ────────────────────────────────────────────────────
  cat("Running SingleR annotation...\n")
  test_data <- GetAssayData(data_obj, assay = "SCT", layer = "data")

  if (single_cell) {
    pred <- SingleR(
      test   = test_data,
      ref    = ref,
      labels = ref$label.main
    )
    data_obj$SingleR_label <- pred$labels
    cat("Cell type distribution:\n")
    print(sort(table(data_obj$SingleR_label), decreasing = TRUE))

    # UMAP colored by cell type
    p_umap <- DimPlot(data_obj, reduction = "umap", group.by = "SingleR_label", label = TRUE)
    png(
      file.path(out_path, paste0(data_name, "_annotated_umap.png")),
      width = 1500, height = 1000, res = 150
    )
    print(p_umap)
    dev.off()

    # Spatial (image) plot colored by cell type
    p_spatial <- ImageDimPlot(data_obj, group.by = "SingleR_label", size = 0.4)
    png(
      file.path(out_path, paste0(data_name, "_annotated_spatial.png")),
      width = 1400, height = 1100, res = 150
    )
    print(p_spatial)
    dev.off()

  } else {
    # Cluster-level: one label per cluster
    pred <- SingleR(
      test     = test_data,
      ref      = ref,
      labels   = ref$label.main,
      clusters = Idents(data_obj)
    )
    cluster_labels <- pred$labels
    names(cluster_labels) <- levels(data_obj)
    data_obj <- RenameIdents(data_obj, cluster_labels)

    # Store label in metadata for downstream access
    data_obj$SingleR_label <- as.character(Idents(data_obj))

    cat("Cluster label assignments:\n")
    print(cluster_labels)

    p_clusters    <- DimPlot(data_obj, reduction = "umap", label = TRUE)
    p_annotated   <- DimPlot(data_obj, reduction = "umap", label = TRUE)
    combined      <- p_clusters / p_annotated

    png(
      file.path(out_path, paste0(data_name, "_annotated_cluster_umap.png")),
      width = 900, height = 1400, res = 150
    )
    print(combined)
    dev.off()

    p_spatial <- ImageDimPlot(data_obj, size = 0.4)
    png(
      file.path(out_path, paste0(data_name, "_annotated_cluster_spatial.png")),
      width = 1400, height = 1100, res = 150
    )
    print(p_spatial)
    dev.off()
  }

  cat("=============== Annotation done for:", data_name, "===============\n")
  return(data_obj)
}
