# Different Files in a Dataset
This file describes the different files in GSE281978 (Visium) and GSE300147 (Xenium)

## 10x Genomics Visium Spatial Transcriptomics of HNSCC
- ⭐ **filtered_feature_bc_matrix**: main gene expression matrix that contains Genes × spots matrix, Spot barcodes, Counts per gene per spot
- 🖼️ **tissue_positions_list.csv**: Maps spots to coordinates that contains Spot barcode, x/y coordinates, and Whether the spot is on tissue or not
  - Used for SpatialFeaturePlot and Mapping expression to tissue
- 🖼️ **tissue_hires_image.png**: High-resolution tissue image Used for Overlay gene expression on tissue
- 🖼️ **tissue_lowres_image.png**: Low-resolution version of the tissue image Used for Quick visualization
- 🖼️ **scalefactors_json.json**: Scaling info between image + coordinates thqat Convert pixel coordinates to spot coordinates
  - it Ensures spots align correctly on the image
- **aligned_fiducials.jpg**: Alignment reference markers that Shows fiducial spots used during imaging
  - Mostly used for Quality control, Not used in most analyses
- **detected_tissue_image.jpg**: Binary tissue mask that Shows which areas are tissue vs background
  - Helps filter out non-tissue spots

## Xenium Spatial Transcriptomics of HNSCC
- ⭐ **cell_feature_matrix.h5**: main gene expression matrix that contains Genes × cells count matrix, Cell barcodes, and Gene names
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
