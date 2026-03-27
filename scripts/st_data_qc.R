source("st_helper.R")
library(Seurat)
library(patchwork)

# INPUT: 
#   (str) data folder containing the .h5 file
#   (str) name of the data set 
#   (bool) if graphical output is desired (true by default) 
#   (str) output path for graphical output (same as data folder by default)
# OUTPUT: SeuratObject of .h5 after qc
# DATA FOLDER FORMAT:
#   (Folder) Sample name
#   -> (File) (SampleName)_filtered_feature_bc_matrix.h5
#   -> (Folder) (SampleName)
#   ->-> (File) tissue_lowres_image.png
#   ->-> (File) scalefactors_json.json
#   ->-> (File) tissue_positions_list.csv
# make sure corresponding libraries are already downloaded

run_qc <- function(data_path, data_name, graph_output = TRUE, out_path = data_path) {
  
  print("----------Starting QC----------")
  
  #set crucial vars
  dataset_name <- data_name
  data.dir <- data_path
  graphical_output <- graph_output
  output_path <- out_path
  
  #qc
  print("running qc")
  seuObj <- Load10X_Spatial(data.dir, filename="GSM8633891_filtered_feature_bc_matrix.h5")
  seuObj[["percent.mt"]] <- PercentageFeatureSet(object = seuObj, pattern = "^MT-")
  seuObj[["percent.ribo"]] <- PercentageFeatureSet(seuObj, pattern = "^RP[SL]")
  sub_seuObj <- subset(seuObj, subset = nFeature_Spatial < 7500 & nFeature_Spatial > 200 & nCount_Spatial < 50000 & nCount_Spatial > 250 & percent.mt < 15 & percent.ribo < 40)
  norm_sub_seuObj <- SCTransform(sub_seuObj, assay = "Spatial")
  
  #plotting
  if (graphical_output){
    print("plotting")
    row_len = 4
    features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo")
    
    #generate each graph
    p1 <- make_v_plot(seuObj, features, row_len)
    p2 <- make_v_plot(sub_seuObj, features, row_len)
    s2 <- make_sf_plot(sub_seuObj, features, row_len)
    p3 <- make_v_plot(norm_sub_seuObj, features, row_len)
    s3 <- make_sf_plot(norm_sub_seuObj, features, row_len)
    
    #create title for each graph
    t1 <- title_plot("Raw Data")
    t2 <- title_plot("After QC Processing")
    t3 <- title_plot("After Normalizing")
    
    # Combine the three plots vertically
    combined <- t1 / p1 / t2 / p2 / s2 / t3 / p3 / s3 # "/" stacks vertically; "|" stacks horizontally
    
    # Save to PNG
    png(paste0(output_path, "/", dataset_name, "_qc_plots.png"), width = 2000, height = 3000, res = 150)  # adjust size as needed
    print(combined)
    dev.off()
  
  } 
  rm(seuObj, sub_seuObj, combined, p1, p2, p3, s2, s3, t1, t2, t3)
  print("----------QC done!----------")
  return(norm_sub_seuObj)
}

