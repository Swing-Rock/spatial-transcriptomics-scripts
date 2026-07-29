import scanpy as sc
import numpy as np
import pandas as pd
import cell2location
import torch
import squidpy as sq


#tutorial: https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_tutorial.html

# Device configuration
USE_GPU = torch.cuda.is_available()
if USE_GPU:
    print(f"GPU detected: {torch.cuda.get_device_name(0)}", flush=True)
    accelerator = "gpu"
else:
    print("No GPU detected, falling back to CPU", flush=True)
    accelerator = "cpu"

adata = sq.read.visium(
    path="CRC_project/GSE226997/GSM7089855",
    counts_file="filtered_feature_bc_matrix.h5",
    library_id="GSE226997"
)

#Writes AnnData object into .h5ad format
adata.write_h5ad("CRC_project/GSE226997/GSM7089855/processed.h5ad")

print("Step1 read data", flush=True)
#1 Read data
inf_aver = pd.read_csv( #trained model posteriors from previous step
    "CRC_project/inf_aver.csv",
    index_col=0 #so that gene names read as rownames instead of own col
)
adata_vis = sc.read_h5ad("CRC_project/GSE226997/GSM7089855/processed.h5ad")
print(adata_vis.var.columns.tolist(), flush=True) 

#2 filter out mitochondrial genes ("technical artifacts")
adata_vis.var_names_make_unique() #some genes may have non unique names so make them unique
print("Step2 filter MT genes", flush=True)
# find mitochondria-encoded (MT) genes
adata_vis.var['MT_gene'] = [gene.startswith('MT-') for gene in adata_vis.var_names]

# remove MT genes for spatial mapping (keeping their counts in the object)
adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT_gene'].values].X.toarray()
adata_vis = adata_vis[:, ~adata_vis.var['MT_gene'].values]


print("Step3 find shared genes", flush=True)
#3 find shared genes between reference and visium sample
# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

print("Step4 prepare visium for model", flush=True)
#4 prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key=None)

print("Step5 create model", flush=True)
#5 create model
mod = cell2location.models.Cell2location(
    adata_vis,
    cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=8,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=200
)

print("Step6 train model", flush=True)
#6 train model
mod.train(max_epochs=30000,
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1,
          accelerator=accelerator
         )
#saving the model data at designated directory
mod.save("CRC_project/GSE226997/GSM7089855/sample1_model", overwrite=True)
print("model saved in directory", flush=True)

print("Step export h5ad", flush=True)
#7 export posterior and add onto adata_vis object
adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs}
)
#saving anndata file with results (overwrites originally loaded h5ad)
adata_vis.write_h5ad("CRC_project/GSE226997/GSM7089855/GSM7089855_deconvolved.h5ad") #writes the anndata into h5ad
print("h5ad file with model results saved in directory", flush=True)

print("Step8 extract cell type abundances", flush=True)
#8 extract mean cell abundance per spot as dataframe
abundance_df = adata_vis.obsm['means_cell_abundance_w_sf']#takes posterior values from model
abundance_df = pd.DataFrame( #stores values in pandas dataframe
    abundance_df,
    index=adata_vis.obs_names,
    columns=adata_vis.uns['mod']['factor_names']
)

abundance_df_5 = adata_vis.obsm['q05_cell_abundance_w_sf']#takes 5% quantile posterior values from model
abundance_df_5 = pd.DataFrame( #stores values in pandas dataframe
    abundance_df_5,
    index=adata_vis.obs_names,
    columns=adata_vis.uns['mod']['factor_names']
)

# save to csv #saves this dataframe to a csv
abundance_df.to_csv("CRC_project/GSE226997/GSM7089855/GSM7089855_abundance.csv")
abundance_df_5.to_csv("CRC_project/GSE226997/GSM7089855/GSM7089855_conservative_abundance.csv")
print(f"Abundance table shape: {abundance_df.shape}", flush=True)
print("Done", flush=True)
