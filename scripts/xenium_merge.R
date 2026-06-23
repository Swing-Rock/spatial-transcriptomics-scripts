source("visium_scripts/st_helper.R")
library(Seurat)
library(harmony)
library(patchwork)

# Merge a named list of annotated Xenium Seurat objects and integrate with Harmony.
#
# INPUT:
#   (list)  data_obj  — named list of Seurat objects (post-QC + annotation)
#   (str)   out_path  — directory for PNG output
# OUTPUT: merged Seurat object with both pre- and post-Harmony reductions

merge_xenium <- function(data_obj, out_path, verbose_output = FALSE) {

  cat("=============== Merging Xenium samples ===============\n")

  # ── Merge ─────────────────────────────────────────────────────────────────
  merged_obj <- merge(
    x            = data_obj[[1]],
    y            = data_obj[-1],
    add.cell.ids = names(data_obj),
    project      = "MergedXenium"
  )

  # After merging, SCT layers need to be rejoined and re-run.
  # PrepSCTFindMarkers / JoinLayers won't recalculate residuals across samples,
  # so we re-run SCTransform on the merged object to get a shared feature space.
  cat("--- Re-running SCTransform on merged object ---\n")
  merged_obj <- SCTransform(merged_obj, assay = "Xenium", verbose = verbose_output)

  # ── Standard workflow (before batch correction) ───────────────────────────
  cat("--- Running PCA / neighbors / clusters / UMAP (pre-Harmony) ---\n")
  merged_obj <- FindVariableFeatures(merged_obj, assay = "SCT",
                                     selection.method = "vst", nfeatures = 3000)
  merged_obj <- RunPCA(merged_obj, assay = "SCT",
                       reduction.name = "pca.SCT", verbose = verbose_output)
  merged_obj <- FindNeighbors(merged_obj, assay = "SCT",
                               reduction = "pca.SCT", dims = 1:15,
                               verbose = verbose_output)
  merged_obj <- FindClusters(merged_obj, resolution = 0.5,
                              cluster.name = "seurat_cluster.SCT",
                              verbose = verbose_output)
  merged_obj <- RunUMAP(merged_obj, reduction = "pca.SCT",
                        reduction.name = "umap.SCT",
                        return.model = TRUE, dims = 1:15,
                        verbose = verbose_output)

  # ── Harmony batch correction ──────────────────────────────────────────────
  cat("--- Running Harmony ---\n")
  merged_obj <- RunHarmony(
    object         = merged_obj,
    group.by.vars  = "orig.ident",
    assay.use      = "SCT",
    reduction      = "pca.SCT",
    reduction.save = "harmony.SCT",
    theta          = 8,
    verbose        = verbose_output
  )

  # ── Re-cluster on Harmony embedding ──────────────────────────────────────
  cat("--- Re-clustering on Harmony embedding ---\n")
  merged_obj <- FindNeighbors(merged_obj, reduction = "harmony.SCT",
                               dims = 1:15, verbose = verbose_output)
  merged_obj <- FindClusters(merged_obj, resolution = 0.5,
                              cluster.name = "seurat_cluster.harmony.SCT",
                              verbose = verbose_output)
  merged_obj <- RunUMAP(merged_obj, reduction = "harmony.SCT",
                        reduction.name = "umap.harmony.SCT",
                        dims = 1:15, verbose = verbose_output)

  # ── Visualization ─────────────────────────────────────────────────────────
  cat("--- Generating merge/integration plots ---\n")

  p1 <- DimPlot(merged_obj, reduction = "umap.SCT",
                group.by = "orig.ident", label = TRUE)
  p2 <- DimPlot(merged_obj, reduction = "umap.harmony.SCT",
                group.by = "orig.ident", label = TRUE)
  p3 <- DimPlot(merged_obj, reduction = "umap.harmony.SCT",
                group.by = "SingleR_label", label = TRUE)

  t1 <- title_plot("Before Harmony (by sample)", 17)
  t2 <- title_plot("After Harmony (by sample)", 17)
  t3 <- title_plot("After Harmony (by cell type)", 17)

  combined <- t1 / p1 / t2 / p2 / t3 / p3

  png(
    file.path(out_path, "xenium_merged_cluster_plots.png"),
    width = 1000, height = 2600, res = 150
  )
  print(combined)
  dev.off()

  cat("=============== Finished merging Xenium samples ===============\n")
  return(merged_obj)
}
