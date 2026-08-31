source("visium_scripts/st_helper.R")
library(Seurat)
library(patchwork)

# INPUT: 
#   (str) data folder containing the .h5 file, OR a path to a single .rda file
#   (str) name of the data set 
#   (bool) if graphical output is desired (true by default) 
#   (str) output path for graphical output (same as data folder by default)
# OUTPUT: SeuratObject of .h5 after qc
# make sure corresponding libraries are already downloaded
# GSM8655157 - Load from MTX + coordinates



detect_data_format <- function(data_path) {
  files <- list.files(data_path)
  
  if (any(grepl("\\.mtx(\\.gz)?$", files))) {
    return("mtx")
  } else if (any(grepl("filtered_feature_bc_matrix\\.h5$", files))) {
    return("spaceranger_h5")
  } else if (any(grepl("\\.h5ad$", files))) {
    return("h5ad")
  } else {
    return("unknown")
  }
}

library(rhdf5)
library(Matrix)

write_10x_h5 <- function(counts, ensembl_ids, gene_symbols, output_path, genome = "GRCh38") {
  
  # Ensure sparse column-compressed format (required for indices/indptr)
  counts <- as(counts, "CsparseMatrix")
  
  if (file.exists(output_path)) file.remove(output_path)
  h5createFile(output_path)
  
  h5createGroup(output_path, "matrix")
  h5createGroup(output_path, "matrix/features")
  
  # Core matrix components (CSC format)
  h5write(counts@x, output_path, "matrix/data")
  h5write(counts@i, output_path, "matrix/indices")       # 0-based row indices
  h5write(counts@p, output_path, "matrix/indptr")        # column pointers
  h5write(c(nrow(counts), ncol(counts)), output_path, "matrix/shape")
  
  # Barcodes and features
  h5write(colnames(counts), output_path, "matrix/barcodes")
  h5write(ensembl_ids, output_path, "matrix/features/id")
  h5write(gene_symbols, output_path, "matrix/features/name")
  h5write(rep("Gene Expression", nrow(counts)), output_path, "matrix/features/feature_type")
  h5write(rep(genome, nrow(counts)), output_path, "matrix/features/genome")
  
  cat("✓ Wrote 10x-style H5 to:", output_path, "\n")
  cat("  Genes:", nrow(counts), "| Spots:", ncol(counts), "\n")
}

load_mtx_spatial <- function(mtx_path, coords_path, scale_path, image_path, features_path) {
  library(Matrix)
  library(Seurat)
  library(jsonlite)
  library(png)
  
  mtx <- readMM(mtx_path)
  
  # Load real gene symbols from features file
  features_df <- read.delim(features_path, header = FALSE)
  gene_names <- features_df[, 2]  # Gene symbols (column 2)
  
  counts <- as(mtx, "CsparseMatrix")
  rownames(counts) <- gene_names
  
  # Sparse-native aggregation of duplicate gene symbols
  unique_genes <- unique(gene_names)
  # Build a sparse indicator matrix: rows = unique genes, cols = original genes
  group_idx <- match(gene_names, unique_genes)
  indicator <- sparseMatrix(i = group_idx, 
                            j = seq_along(gene_names), 
                            x = 1,
                            dims = c(length(unique_genes), length(gene_names)))
  
  # Multiply to sum duplicate rows (stays sparse the whole time)
  counts <- indicator %*% counts
  rownames(counts) <- unique_genes
  counts <- as(counts, "CsparseMatrix")
  colnames(counts) <- paste0("spot_", seq_len(ncol(counts)))
  
  # After your aggregation step, counts has unique gene symbols as rownames
  write_10x_h5(
    counts = counts,
    ensembl_ids = rownames(counts),   # or map back to first ENSG per symbol if you kept that mapping
    gene_symbols = rownames(counts),
    output_path = file.path(dirname(mtx_path), "filtered_feature_bc_matrix.h5")
  )
  
  # 2. Load spatial coordinates
  coords_df <- read.csv(coords_path, row.names = 1, header = FALSE)
  colnames(coords_df) <- c("in_tissue", "array_row", "array_col", "pxl_row", "pxl_col")
  
  # 3. Load scale factors
  scale_factors <- fromJSON(scale_path)
  
  # 4. Load image
  img <- readPNG(image_path)
  
  # 5. Create Seurat object
  obj <- CreateSeuratObject(counts = counts, assay = "Spatial")
  
  # Add spatial metadata
  obj@meta.data <- cbind(obj@meta.data, coords_df)
  
  # # Register the image
  # if (!is.null(img)) {
  #   image_obj <- Image(image = img, scale.factors = scale_factors)
  #   obj@images[["image"]] <- image_obj
  # }
  
  # Register image (manually, since we're not using Load10X_Spatial)
  # This is tricky—you may need to use a VisiumV1/V2 image object
  
  return(obj)
}

run_qc <- function(data_path, data_name, graph_output = TRUE, out_path = data_path, normalize = TRUE, verbose_output = FALSE) {
  
  cat("===============Starting QC for sample", data_name, "===============\n")
  
  #set crucial vars
  dataset_name <- data_name
  graphical_output <- graph_output
  output_path <- out_path
  
  #qc
  cat("Detecting data format...\n")
  format <- detect_data_format(data_path)
  cat("Format detected:", format, "\n")
  
  # ===== LOAD DATA BASED ON FORMAT =====
  if (format == "mtx") {
    cat("Loading MTX format...\n")
    
    # Find MTX files in current sample folder
    mtx_file <- list.files(data_path, pattern = ".*matrix\\.mtx(\\.gz)?$", full.names = TRUE)[1]
    coords_file <- list.files(data_path, pattern = ".*tissue_positions.*\\.csv(\\.gz)?$", full.names = TRUE)[1]
    scale_file <- list.files(data_path, pattern = ".*scalefactors.*\\.json(\\.gz)?$", full.names = TRUE)[1]
    image_file <- list.files(data_path, pattern = ".*tissue_hires_image\\.(png|jpg)(\\.gz)?$", full.names = TRUE)[1]
    
    # Find shared annotation files in parent GSE folder
    parent_dir <- dirname(data_path)
    features_file <- list.files(parent_dir, pattern = ".*features\\.tsv(\\.gz)?$", full.names = TRUE)[1]
    barcodes_file <- list.files(parent_dir, pattern = ".*barcodes\\.tsv(\\.gz)?$", full.names = TRUE)[1]
    
    if (is.na(mtx_file) || is.na(coords_file) || is.na(scale_file)) {
      stop("MTX sample missing required files in:", data_path)
    }
    if (is.na(features_file)) {
      stop("Features file not found in parent GSE folder:", parent_dir)
    }
    
    seuObj <- load_mtx_spatial(mtx_file, coords_file, scale_file, image_file, features_file)
    
  } else if (format == "spaceranger_h5") {
    cat("Loading Space Ranger H5 format...\n")
    h5_path <- file.path(data_path, "filtered_feature_bc_matrix.h5")
    seuObj <- Load10X_Spatial(data_path, filename = "filtered_feature_bc_matrix.h5")
    
  } else if (format == "h5ad") {
    cat("Loading H5AD format...\n")
    seuObj <- load_h5ad_visium(data_path)
    
  } else {
    stop("No valid input file found in: ", data_path,
         "\n  Expected: MTX format (matrix.mtx + tissue_positions.csv + scalefactors.json),\n",
         "            Space Ranger format (filtered_feature_bc_matrix.h5),\n",
         "            or H5AD format (.h5ad file)")
  }   
  
  
  seuObj[["percent.mt"]] <- PercentageFeatureSet(object = seuObj, pattern = "^MT-")
  seuObj[["percent.ribo"]] <- PercentageFeatureSet(seuObj, pattern = "^RP[SL]")
  # Spots with 0 total counts produce NaN (0/0) — replace with 0 so VlnPlot doesn't crash
  seuObj$percent.mt[is.na(seuObj$percent.mt)]     <- 0
  seuObj$percent.ribo[is.na(seuObj$percent.ribo)] <- 0
  sub_seuObj <- subset(seuObj, subset = nFeature_Spatial < 7500 & nFeature_Spatial > 200 & nCount_Spatial < 50000 & nCount_Spatial > 250 & percent.mt < 15 & percent.ribo < 40)
  #normalize
  if(normalize == TRUE){
    cat("normalizing\n")
    norm_sub_seuObj <- SCTransform(sub_seuObj, assay = "Spatial", verbose = verbose_output)
  }
  
  #plotting
  has_images <- length(Images(seuObj)) > 0
  if (graphical_output && !has_images){
    cat("skipping spatial feature plots for", data_name,
        "- no @images slot present (likely stripped from this .rda). ",
        "Violin-plot QC metrics still saved.\n")
  }
  if (graphical_output){
    cat("plotting\n")
    row_len = 4
    features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo")
    
    #generate each graph
    p1 <- make_v_plot(seuObj, features, row_len)
    p2 <- make_v_plot(sub_seuObj, features, row_len)
    
    t1 <- title_plot("Raw Data")
    t2 <- title_plot("After QC Processing")
    
    combined <- t1 / p1 / t2 / p2
    if (has_images) {
      s2 <- make_sf_plot(sub_seuObj, features, row_len)
      combined <- t1 / p1 / t2 / p2 / s2
    }

    if (normalize == TRUE){
      p3 <- make_v_plot(norm_sub_seuObj, features, row_len)
      t3 <- title_plot("After Normalizing")

      if (has_images) {
        s3 <- make_sf_plot(norm_sub_seuObj, features, row_len)
        combined <- combined / t3 / p3 / s3 # "/" stacks vertically; "|" stacks horizontally
      } else {
        combined <- combined / t3 / p3
      }
    }
    
    
    # Save to PNG
    png(paste0(output_path, "/", dataset_name, "_qc_plots.png"), width = 2000, height = 3000, res = 150)  # adjust size as needed
    print(combined)
    dev.off()
  
  } 
  #rm(seuObj, sub_seuObj, combined, p1, p2, p3, s2, s3, t1, t2, t3)
  cat("=============== QC done for sample", data_name, "===============\n")
  if (normalize == TRUE){
    return(norm_sub_seuObj)
  } else{
    return(sub_seuObj)
  }
  
}