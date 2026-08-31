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
  # ── Strip per-sample SCT artifacts before merging ─────────────────────────
  # SCTransform is re-run on the merged object, so per-sample SCT assays,
  # reductions, and graphs are wasted memory. Removing them before merge
  # significantly reduces peak RAM usage.
  # Assays to keep; everything else is dropped before the merge to cut
  # peak RAM. BlankCodeword / ControlCodeword / ControlProbe are QC-only
  # and not used downstream. SCT is re-run on the merged object.
  # FOV images hold per-transcript coordinates and are the single largest
  # memory consumer (~40-50 % of each object); they are not needed for
  # clustering / Harmony integration.
  ASSAYS_TO_KEEP <- "Xenium"

  cat("--- Stripping per-sample SCT/reductions/images to free memory ---\n")
  sample_names <- names(data_obj)
  for (nm in sample_names) {
    obj <- data_obj[[nm]]

    # Coerce to character: Seurat v5 + Bioconductor packages can cause
    # Assays() / Reductions() to return S4 objects that break %in% / match.
    assay_names     <- as.character(SeuratObject::Assays(obj))
    reduction_names <- as.character(SeuratObject::Reductions(obj))

    # Drop every assay except Xenium counts
    for (a in setdiff(assay_names, ASSAYS_TO_KEEP)) {
      cat(a, " ")
      if (a != "SCT"){
        obj[[a]] <- NULL
      }
      cat("\n")
    }
    # Drop all reductions and graphs
    for (r in reduction_names)   obj[[r]] <- NULL
    for (g in names(obj@graphs)) obj[[g]] <- NULL
    # Drop spatial image data (FOV transcript coords / cell polygons)
    # — largest single component, not needed for Harmony integration

    data_obj[[nm]] <- obj
    gc()
  }
  # ── Merge ─────────────────────────────────────────────────────────────────
  merged_obj <- merge(
    x            = data_obj[[1]],
    y            = data_obj[-1],
    add.cell.ids = names(data_obj),
    project      = "MergedXenium"
  )
  gc()

  
  # 
  # # ── Sequential merge ──────────────────────────────────────────────────────
  # # Merge one sample at a time and immediately free it from the list so the
  # # garbage collector can reclaim memory before the next iteration.
  # # Peak RAM = size(merged so far) + size(next sample), rather than holding
  # # the entire list plus the full merged object simultaneously.
  # #
  # # Cell barcodes are prefixed with the sample name up-front (one object at
  # # a time) so that cross-sample collisions cannot occur during the sequential
  # # merge and add.cell.ids can be left NULL in each merge() call.
  # cat("--- Merging samples sequentially ---\n")
  # 
  # # Prefix and pop the first sample
  # obj1 <- data_obj[[sample_names[1]]]
  # obj1 <- RenameCells(obj1, new.names = paste0(sample_names[1], "_", Cells(obj1)))
  # data_obj[[sample_names[1]]] <- NULL
  # gc()
  # 
  # merged_obj <- obj1
  # rm(obj1); gc()
  # 
  # for (nm in sample_names[-1]) {
  #   cat("  merging:", nm, "\n")
  #   next_obj <- data_obj[[nm]]
  #   next_obj <- RenameCells(next_obj,
  #                           new.names = paste0(nm, "_", Cells(next_obj)))
  #   data_obj[[nm]] <- NULL
  #   gc()
  # 
  #   merged_obj <- merge(merged_obj, y = next_obj, project = "MergedXenium")
  #   rm(next_obj); gc()
  # }
  # rm(data_obj)
  # gc()

  # After merging, SCT layers need to be rejoined and re-run.
  # PrepSCTFindMarkers / JoinLayers won't recalculate residuals across samples,
  # so we re-run SCTransform on the merged object to get a shared feature space.
  cat("--- Re-running SCTransform on merged object ---\n")
  merged_obj <- SCTransform(merged_obj, assay = "Xenium", verbose = verbose_output)

  # ── Standard workflow (before batch correction) ───────────────────────────
  cat("--- Running PCA / neighbors / clusters / UMAP (pre-Harmony) ---\n")
  cat("Running FindVariableFeatures \n")
  merged_obj <- FindVariableFeatures(merged_obj, assay = "SCT",
                                     selection.method = "vst", nfeatures = 3000)
  cat("Running RunPCA \n")
  merged_obj <- RunPCA(merged_obj, assay = "SCT",
                       reduction.name = "pca.SCT", verbose = verbose_output)
  cat("Running FindNeighbors \n")
  merged_obj <- FindNeighbors(merged_obj, assay = "SCT",
                               reduction = "pca.SCT", dims = 1:15,
                               verbose = verbose_output)
  cat("Running FindClusters \n")
  merged_obj <- FindClusters(merged_obj, resolution = 0.5,
                              cluster.name = "seurat_cluster.SCT",
                              verbose = verbose_output)
  cat("Running RunUMAP \n")
  merged_obj <- RunUMAP(merged_obj, reduction = "pca.SCT",
                        reduction.name = "umap.SCT",
                        return.model = TRUE, dims = 1:15,
                        verbose = verbose_output)
  gc()
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
  gc()
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
  gc()

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
