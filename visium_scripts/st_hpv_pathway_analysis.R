# ============================================================
# st_hpv_pathway_analysis.R
# Differential expression + GO/KEGG pathway enrichment
# between HPV+ and HPV- malignant spots (Visium)
# ============================================================
# Dependencies: Seurat, clusterProfiler, org.Hs.eg.db, enrichplot, ggplot2
# Install if needed:
#   BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))
# ============================================================

library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)

# ── Config ────────────────────────────────────────────────────────────────────

# Cell type label column in merged_data@meta.data

# Which labels count as malignant — adjust to match your annotation
MALIGNANT_LABELS <- c("Epithelial_cells", "Keratinocytes")

# FDR and fold-change thresholds for DE
FDR_CUTOFF <- 0.05
LFC_CUTOFF <- 0.5      # log2 fold-change

# Number of top pathways to plot
N_PATHWAYS <- 20


# ── Step 1: Annotate HPV status ───────────────────────────────────────────────

add_hpv_status <- function(seurat_obj, hpv_plus_samples) {
  cat("----------adding hpv labels----------")
  seurat_obj$hpv_status <- ifelse(
    seurat_obj$orig.ident %in% hpv_plus_samples,
    "HPV_pos", "HPV_neg"
  )
  cat("HPV status distribution:\n")
  print(table(seurat_obj$hpv_status, seurat_obj$orig.ident))
  return(seurat_obj)
}


# ── Step 2: Subset to malignant spots ────────────────────────────────────────
subset_malignant <- function(seurat_obj, deconv, malignant_labels = NULL) {
  cat("----------subsetting malignant spots----------")
  seurat_obj <- PrepSCTFindMarkers(seurat_obj)
  
  if (deconv == FALSE){
    cell_type_col <- "SingleR_label"
    # Check what labels actually exist first
    all_labels <- unique(seurat_obj@meta.data[[cell_type_col]])
    cat("All labels found:\n")
    print(all_labels)
    
    # Use direct metadata indexing instead of expression=
    cells_keep <- rownames(seurat_obj@meta.data)[
      seurat_obj@meta.data[[cell_type_col]] %in% malignant_labels
    ]
    
    cat("Cells matching malignant labels:", length(cells_keep), "\n")
    
    if (length(cells_keep) == 0) {
      stop("No cells matched MALIGNANT_LABELS — check the label names above and update MALIGNANT_LABELS")
    }
    
  }else{
    prop_col  <- "RCTD_Malignant"
    threshold  <- 0.5
    
    props <- seurat_obj@meta.data[[prop_col]]
    cells_keep <- rownames(seurat_obj@meta.data)[!is.na(props) & props >= threshold]
    
    cat("Malignant spots (RCTD_Malignant >=", threshold, "):", length(cells_keep), "\n")
    
    if (length(cells_keep) == 0) {
      stop("No spots passed the malignant threshold — try lowering threshold")
    }
  }
  
  
  
  malignant <- subset(seurat_obj, cells = cells_keep)
  cat("Malignant spots by HPV status:\n")
  print(table(malignant$hpv_status))
  return(malignant)
}


# ── Step 3: Differential expression (HPV+ vs HPV-) ───────────────────────────

run_de <- function(malignant, out_path) {
  cat("----------running DE----------")
  
  Idents(malignant) <- "hpv_status"
  
  de_results <- FindMarkers(
    malignant,
    ident.1         = "HPV_pos",
    ident.2         = "HPV_neg",
    assay           = "Spatial",
    layer           = "data",    # ← "slot" → "layer" in Seurat v5
    test.use        = "wilcox",
    min.pct         = 0.1,
    logfc.threshold = LFC_CUTOFF
  )

  de_results$gene <- rownames(de_results)
  de_results <- de_results[de_results$p_val_adj < FDR_CUTOFF, ]
  de_results <- de_results[order(de_results$avg_log2FC, decreasing = TRUE), ]

  write.csv(de_results, file.path(out_path, "hpv_de_results.csv"), row.names = TRUE)
  cat("Significant DE genes:", nrow(de_results), "\n")
  cat("  Upregulated in HPV+:", sum(de_results$avg_log2FC > 0), "\n")
  cat("  Upregulated in HPV-:", sum(de_results$avg_log2FC < 0), "\n")

  # Volcano plot
  de_results$direction <- ifelse(de_results$avg_log2FC > 0, "HPV+", "HPV-")
  p <- ggplot(de_results, aes(x = avg_log2FC, y = -log10(p_val_adj), color = direction)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("HPV+" = "#E64B35", "HPV-" = "#4DBBD5")) +
    geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(FDR_CUTOFF), linetype = "dashed", color = "grey50") +
    labs(title = "HPV+ vs HPV- DE (Malignant Spots)",
         x = "Log2 Fold Change", y = "-Log10 Adjusted P-value",
         color = "Enriched in") +
    theme_classic()
  ggsave(file.path(out_path, "hpv_volcano.png"), p, width = 7, height = 5, dpi = 150)

  return(de_results)
}


# ── Step 4: GO enrichment ─────────────────────────────────────────────────────

run_go <- function(de_results, out_path) {
  cat("----------running go enrichment----------")
  
  # Convert gene symbols to Entrez IDs
  gene_df <- bitr(de_results$gene, fromType = "SYMBOL",
                  toType = "ENTREZID", OrgDb = org.Hs.eg.db)

  # Separate HPV+ and HPV- upregulated genes
  hpv_pos_genes <- de_results$gene[de_results$avg_log2FC > 0]
  hpv_neg_genes <- de_results$gene[de_results$avg_log2FC < 0]

  run_enrichgo <- function(genes, label) {
    entrez <- gene_df$ENTREZID[gene_df$SYMBOL %in% genes]
    if (length(entrez) < 5) {
      cat("  Not enough genes for GO enrichment in:", label, "\n")
      return(NULL)
    }
    ego <- enrichGO(
      gene          = entrez,
      OrgDb         = org.Hs.eg.db,
      ont           = "BP",          # Biological Process
      pAdjustMethod = "BH",
      pvalueCutoff  = FDR_CUTOFF,
      readable      = TRUE
    )
    if (is.null(ego) || nrow(ego@result) == 0) {
      cat("  No significant GO terms for:", label, "\n")
      return(NULL)
    }

    write.csv(ego@result, file.path(out_path, paste0("go_", label, ".csv")), row.names = FALSE)

    # Dotplot
    p <- dotplot(ego, showCategory = N_PATHWAYS, title = paste("GO BP —", label)) +
      theme(axis.text.y = element_text(size = 8))
    ggsave(file.path(out_path, paste0("go_dotplot_", label, ".png")),
           p, width = 8, height = 8, dpi = 150)

    return(ego)
  }

  ego_pos <- run_enrichgo(hpv_pos_genes, "HPV_pos")
  ego_neg <- run_enrichgo(hpv_neg_genes, "HPV_neg")

  return(list(HPV_pos = ego_pos, HPV_neg = ego_neg))
}


# ── Step 5: KEGG enrichment ───────────────────────────────────────────────────

run_kegg <- function(de_results, out_path) {
  cat("----------running kegg enrichment----------")
  
  gene_df <- bitr(de_results$gene, fromType = "SYMBOL",
                  toType = "ENTREZID", OrgDb = org.Hs.eg.db)

  hpv_pos_genes <- de_results$gene[de_results$avg_log2FC > 0]
  hpv_neg_genes <- de_results$gene[de_results$avg_log2FC < 0]

  run_enrichkegg <- function(genes, label) {
    entrez <- gene_df$ENTREZID[gene_df$SYMBOL %in% genes]
    if (length(entrez) < 5) return(NULL)

    kk <- enrichKEGG(
      gene          = entrez,
      organism      = "hsa",
      pAdjustMethod = "BH",
      pvalueCutoff  = FDR_CUTOFF
    )
    if (is.null(kk) || nrow(kk@result) == 0) {
      cat("  No significant KEGG pathways for:", label, "\n")
      return(NULL)
    }

    write.csv(kk@result, file.path(out_path, paste0("kegg_", label, ".csv")), row.names = FALSE)

    p <- dotplot(kk, showCategory = N_PATHWAYS, title = paste("KEGG —", label)) +
      theme(axis.text.y = element_text(size = 8))
    ggsave(file.path(out_path, paste0("kegg_dotplot_", label, ".png")),
           p, width = 8, height = 8, dpi = 150)

    return(kk)
  }

  kk_pos <- run_enrichkegg(hpv_pos_genes, "HPV_pos")
  kk_neg <- run_enrichkegg(hpv_neg_genes, "HPV_neg")

  return(list(HPV_pos = kk_pos, HPV_neg = kk_neg))
}


# ── Step 6: Spatial pathway scoring (GSVA / module scoring) ──────────────────

run_spatial_scoring <- function(seurat_obj, go_results, out_path) {
  cat("----------running spatial scoring----------")
  
  # Extract top GO gene sets for HPV+ and HPV- enriched programs
  get_geneset <- function(ego, n_terms = 5) {
    if (is.null(ego)) return(list())
    top_terms <- head(ego@result$geneID[ego@result$p.adjust < FDR_CUTOFF], n_terms)
    gene_lists <- lapply(top_terms, function(x) unlist(strsplit(x, "/")))
    names(gene_lists) <- head(ego@result$Description[ego@result$p.adjust < FDR_CUTOFF], n_terms)
    return(gene_lists)
  }

  pos_sets <- get_geneset(go_results$HPV_pos)
  neg_sets <- get_geneset(go_results$HPV_neg)

  # Score each spot using Seurat AddModuleScore
  score_and_plot <- function(gene_sets, label, obj) {
    for (term_name in names(gene_sets)) {
      genes <- gene_sets[[term_name]]
      genes <- genes[genes %in% rownames(obj)]
      if (length(genes) < 3) next

      safe_name <- gsub("[^A-Za-z0-9]", "_", substr(term_name, 1, 30))
      score_col  <- paste0(label, "_", safe_name)

      obj <- AddModuleScore(obj, features = list(genes),
                            name = score_col, assay = "SCT")
      # Seurat appends "1" to the name
      actual_col <- paste0(score_col, "1")

      n_images <- length(Images(obj))
      n_col    <- min(n_images, 3)
      n_row    <- ceiling(n_images / n_col)
      fig_w    <- n_col * 5
      fig_h    <- n_row * 4 + 1.5   # extra for title

      p <- SpatialFeaturePlot(obj, features = actual_col,
                              pt.size.factor = 1.3, ncol = n_col) +
        labs(title = paste(label, "—", term_name)) +
        theme(
          legend.position = "right",
          plot.title      = element_text(size = 14, face = "bold"),
          strip.text      = element_text(size = 10)
        )
      ggsave(file.path(out_path, paste0("spatial_score_", score_col, ".png")),
             p, width = fig_w, height = fig_h, dpi = 150)
    }
    return(obj)
  }

  seurat_obj <- score_and_plot(pos_sets, "HPVpos", seurat_obj)
  seurat_obj <- score_and_plot(neg_sets, "HPVneg", seurat_obj)

  return(seurat_obj)
}
