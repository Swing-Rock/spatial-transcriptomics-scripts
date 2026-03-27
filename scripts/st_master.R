source("st_data_qc.R")
source("st_clustering.R")

data_folder = "C:/Users/srswi/OneDrive/Desktop/ongkeko lab works/spatial transcriptomics/GSM8633891"
data_name = "GSM8633891"
picture_folder = "C:/Users/srswi/OneDrive/Desktop/ongkeko lab works/spatial transcriptomics/GSM8633891/pngs"

data <- run_qc(data_folder, data_name, TRUE, picture_folder)
clustering(data, data_name, picture_folder)


