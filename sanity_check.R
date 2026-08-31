 data_path_60 <- "CRC_project/GSE326101/GSM9622960"
 
   # Check for Zone.Identifier files first
   zone_files_60 <- list.files(data_path_60, pattern = "Zone\\.Identifier$", recursive = TRUE, full.names = TRUE)
   if (length(zone_files_60) > 0) { file.remove(zone_files_60); cat("Removed", length(zone_files_60), "Zone.Identifier files\n") }
   
     seuObj60 <- Load10X_Spatial(data_path_60, filename = "filtered_feature_bc_matrix.h5")
     seuObj60[["percent.mt"]]   <- PercentageFeatureSet(seuObj60, pattern = "^MT-")
     seuObj60[["percent.ribo"]] <- PercentageFeatureSet(seuObj60, pattern = "^RP[SL]")
     seuObj60$percent.mt[is.na(seuObj60$percent.mt)]     <- 0
     seuObj60$percent.ribo[is.na(seuObj60$percent.ribo)] <- 0
     
       meta60 <- seuObj60@meta.data
       features <- c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo")
       cat("=== Total spots:", nrow(meta60), "===\n\n")
       cat("=== Quantiles ===\n")
       print(round(sapply(meta60[, features], quantile, probs = c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1)), 2))