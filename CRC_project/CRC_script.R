source("visium_scripts/st_data_qc.R")
library(edgeR)
library(MAP)


# define path for everything 
# data folder
data_folder = "CRC_project/GSE226997" 
# graphical output folder
picture_folder = "CRC_project/GSE226997/output"
dir.create(picture_folder, showWarnings = FALSE, recursive = TRUE)
# data checkpoint folder
cache_folder_og  = file.path(picture_folder, "cache") 
dir.create(cache_folder_og, showWarnings = FALSE, recursive = TRUE)

#define operation parameters
#if there's multiple samples in the folder that need to be merged
batch_processing <- TRUE 

start_time <- Sys.time()

#----------all the functions that actually runs the code are here----------

single_sample_qc <- function(data_name){
    data <- run_qc(data_folder, data_name, TRUE, picture_folder, FALSE)
    return (data)
}

# Define patched sigscore using new GSVA API
sigscore_patched <- function(expr, sigDir) {
  tmp <- NULL
  exprrow <- read.delim(expr, header = TRUE, check.names = FALSE, row.names = 1)
  
  for (i in 1:ncol(exprrow)) {
    expsub <- as.matrix(exprrow[, i])
    row.names(expsub) <- row.names(exprrow)
    
    # New GSVA >=1.50 API
    score <- GSVA::gsva(GSVA::ssgseaParam(expsub, gmt_set, normalize = FALSE), verbose = FALSE)
    tmp <- cbind(tmp, score)
  }
  
  colnames(tmp) <- colnames(exprrow)
  es.m <- as.data.frame(cbind(Signature = rownames(tmp), round(tmp, 3)))
  write.table(es.m, paste0(sigDir, "/tmp/ss_ssgsea.txt"),
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  es.m <- read.delim(paste0(sigDir, "/tmp/ss_ssgsea.txt"), header = TRUE, check.names = FALSE)
  maxn <- 10; minn <- 1; score <- NULL
  
  for (i in 1:(ncol(es.m) - 1)) {
    tmp2 <- (maxn - minn) / (max(es.m[, i + 1]) - min(es.m[, i + 1])) *
      (es.m[, i + 1] - min(es.m[, i + 1])) + minn
    score <- rbind(score, tmp2)
  }
  
  score.m <- cbind(Signature = as.character(es.m[, 1]), round(t(score), 2))
  colnames(score.m) <- colnames(es.m)
  write.table(score.m, paste0(sigDir, "/tmp/ssgsea_score1-10.txt"),
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  return(score.m)
}

impute_gene <- function(gene_name, input){
  sample_medians <- apply(input, 2, median, na.rm = TRUE)
  
  cat("Imputing", gene_name, "using per-sample median expression values\n")
  cat("Range of imputed values:", round(range(sample_medians), 4), "\n")
  
  # Build the new row using per-sample medians
  new_row <- data.frame(matrix(sample_medians, nrow = 1, ncol = ncol(input)))
  colnames(new_row) <- colnames(input)
  rownames(new_row) <- gene_name
  
  map_input_fixed <- rbind(input, new_row)
  cat("Fixed matrix dimensions:", nrow(map_input_fixed), "genes x",
      ncol(map_input_fixed), "samples\n")
  
  # Overwrite the input file
  write.table(map_input_fixed, file.path(cache_folder, "MAP_input.txt"),
              sep = "\t", quote = FALSE, col.names = NA)
  cat("MAP_input.txt updated with", gene_name, "\n")
  
  return(map_input_fixed)
}



#---------- end of functions, start main----------
if (batch_processing){
  
  # ==== QC ====
    folders <- list.files(data_folder, full.names = TRUE)
    #initiate variables
    n <- 0
    file_list <- list()
    
    # Loop through folders
    for (folder in folders) {
      # Get folder name only
      folder_name <- basename(folder)
      # Skip the "output" folder
      if (folder_name == basename(picture_folder)) {
        next
      }
      file_list[[folder_name]] <- folder   # key = file name, value = absolute path
      n = n + 1
    }
    cat(n, "SAMPLES WILL BE PROCESSED\n")
    data_list <- list()
    cat("no qc data found, running qc")
    for (name in names(file_list)) {
      cache_folder = paste0(cache_folder_og, "/", name)
      dir.create(cache_folder, showWarnings = FALSE, recursive = TRUE)
      
      cat("\n\n::::::::::::::: Processing sample:", name, " :::::::::::::::\n")
      data <- run_qc(file_list[[name]], name, TRUE, picture_folder, FALSE)
      data_list[[name]] <- data
      data_list[[name]]$orig.ident <- name
      rm(data)
    
    #==== aggregate data for pseudobulking ====
    SeuratObject_sum <- AggregateExpression(
      data_list[[name]],
      group.by = "orig.ident",
      normalization.method = NULL,
      verbose = TRUE
    )
    SeuratObject_sum <- as.matrix(SeuratObject_sum$Spatial) #The output matrix data we want is within "Spatial" of the "SeuratObject_sum", so we invoke that using the $ operator.
    
    if ((dim(data_list[[name]])[1] == dim(SeuratObject_sum)[1]) && (dim(SeuratObject_sum)[2] == 1)){
      print(paste0("Pseudobulking successful, total number of genes: ", dim(SeuratObject_sum)[1]))
    }
    
    
    # ==== Log-normalization via EdgeR ====
    edgeR <- DGEList(counts = SeuratObject_sum)
    keep <- filterByExpr(edgeR, min.count = 50) #this function returns a list of all genes ABOVE the threshold set by "min.count"
    edgeR$counts[!keep, ] <- 0 #this will set all rows (gene names) that are NOT part of the "keep" list (thus any under 50) to 0
    edgeR_cpm <- cpm(edgeR, log = TRUE) #setting log = TRUE gives us log2(CPM+2). The +2 is a default offset in edgeR to prevent log(0), which will throw an error
    
    
    # ===== Prepare results as an input .txt file for MAP =====
    #reference genes list
    MAP_genes <- c("LY6G6D", "CYP2W1", "TNNC2", "CTTNBP2", "NKD1", "CAB39L", "MLH1", "EPM2AIP1",
                   "SHROOM4", "RNF43", "PRR15", "ATP9A", "H2AFJ", "FARP1", "TCF7", "MAPRE3",
                   "ZMYND8", "DDX27", "TGFBR2", "PIWIL4", "FECH", "DOCK5", "TYMS", "HPSE",
                   "ASPHD2", "AGR2", "GFI1", "RPL22L1", "RAB27B", "GNLY", "DUSP4")
    
    #verify 31 genes of interest are in the matrix
    gene_existance = MAP_genes %in% rownames(edgeR_cpm)
    all_genes_exist = TRUE
    missing_gene = list()
    
    for (i in 1:length(gene_existance)){
      if (gene_existance[i] == FALSE) {
        print(paste("Gene", MAP_genes[i], "not found in rownames(edgeR_cpm)"))
        missing_gene <- append(missing_gene, MAP_genes[i])
        all_genes_exist = FALSE
      }
    }
    
    write.table(edgeR_cpm, 
                file.path(cache_folder, "MAP_input.txt"),
                sep = "\t", 
                quote = FALSE, 
                row.names = TRUE, 
                col.names = NA)
    
    if(all_genes_exist){
      print("all MAP_genes are found in rownames(edgeR_cpm)")
    } else{
      map_input <- read.table(file.path(cache_folder, "MAP_input.txt"), 
                              header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
      
      for (g in missing_gene){
        print(paste("imputing", g))
        map_input <- impute_gene(g, map_input)
      }
      rm(map_input)
    }
    
    
    
    # ===== RUNNING MAP =====
    gmt_set <- get("gmt_set", envir = asNamespace("MAP"))
    assignInNamespace("sigscore", sigscore_patched, ns = "MAP")
    cat("sigscore patched successfully.\n")
    
    runMAP(file.path(cache_folder, "MAP_input.txt"), cache_folder)
    }

} else { 
  data_name = basename(data_folder)
  dataObj = single_sample_qc(data_name)
}






end_time <- Sys.time()
cat("PROGRAM START TIME:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("PROGRAM END TIME:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")

total_secs <- as.numeric(difftime(end_time, start_time, units = "secs"))
# Convert to hours, minutes, seconds
hours <- floor(total_secs / 3600)
minutes <- floor((total_secs %% 3600) / 60)
seconds <- round(total_secs %% 60, 2)
# Print nicely
cat("TOTAL RUN TIME: ",
    hours, "h", 
    minutes, "m", 
    seconds, "s\n")
cat("TOTAL RUN TIME: ",
    hours, "h", 
    minutes, "m", 
    seconds, "s\n")

cat("\nALL DONE :D\n")