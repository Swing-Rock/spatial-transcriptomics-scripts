# Different Files in a Dataset
This file describes the different files in GSE281978 (Visium) and GSE300147 (Xenium)

## 10x Genomics Visium Spatial Transcriptomics of HNSCC
  - nFeature_Spatial less than 200 or greater than 7,500
  - nCount_Spatial less than 250 or greater than 50,000
  - percent.mt > 15%
  - percent.ribo > 40%
- normalize data using SCTransform  

## Xenium Spatial Transcriptomics of HNSCC
- **cell_feature_matrix.h5**: main gene expression matrix that contains Genes × cells count matrix, Cell barcodes, and Gene names
  - ⭐ This is the most important file! 
- **cells.parquet**: Metadata for each cell that contains Cell IDs, Coordinates (x, y), Cell area, and QC metrics
  - Links expression data to spatial position
  - Used for Adding metadata to Seurat, Filtering cells, and Spatial plotting
- **cell_boundaries.parquet**: Polygon outlines of each cell that contains the exact shape of each cell and the coordinates of boundaries
  - *You need this to draw real cell shapes on tissue amd do spatial analysis!*
- **nucleus_boundaries.parquet**: similar to cell_boundaries.parquet but for labeling region of nucleus 
- **transcripts.parquet**: Raw transcript-level data where Each row = one RNA molecule, containing gene name, xy positions, and which cell it belongs to
  - This has the most detailed data that's used for Subcellular analysis, Spatial gradients, Custom gene counting, and Rebuilding expression matrix manually
- **morphology.ome.tif**: High-resolution histology image that shows Actual tissue structure with Cells, nuclei, and morphology
  - Used for: Overlaying gene expression on tissue, Visual validation, Identifying tumor regions
