# spatial-transcriptomics-scripts
Ongkeko Lab spatial transcriptomics data processing pipeline

## st_master 
- define path and name variables
- calls other scripts

## st_data_qc
- quality control thresholds:
  - nFeature_Spatial less than 200 or greater than 7,500
  - nCount_Spatial less than 250 or greater than 50,000
  - percent.mt > 15%
  - percent.ribo > 40%
- normalize data using SCTransform

## st_clustering_qc
- uses Human Primary Cell Atlas Data as reference, compares target cell against reference to identify cell type
- **single cell labeling**: each cell is labeled a cell type before dimensional reduction
- **cluster labeling**: each cell is assigned to a cluster before each cluster is given a cell type
    - FindNeighbors dims = 1:30
    - FindCLusters resolution = 0.5

## st_helper
- graphing functions for formatting 

## flow chart of script dependencies


<img width="918" height="806" alt="image" src="https://github.com/user-attachments/assets/7f8e183e-714f-4d6d-aa4a-b10f5dab0ae2" />
https://lucid.app/lucidchart/9ca23dbf-cc43-4b66-b452-ce79e9b11e3e/edit?viewport_loc=-547%2C-325%2C3078%2C1476%2C0_0&invitationId=inv_80ce0e20-4ee8-4184-8f4a-6fc41391dc04

## file organization
Make sure each sample's file structure is organized as follows:
- (Folder) Sample name 
  - (File) (SampleName)_filtered_feature_bc_matrix.h5
  - (Folder) (SampleName)
    - (File) tissue_lowres_image.png
    - (File) scalefactors_json.json
    - (File) tissue_positions_list.csv


lots of thanks and credits to Alfred Kao for providing the code for workflow and many helps
