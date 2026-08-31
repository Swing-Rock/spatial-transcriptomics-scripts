import os
import sys
import json
import shutil
import scanpy as sc
import numpy as np
import pandas as pd
import cell2location
import torch
import squidpy as sq
import matplotlib.image as mpimg

#tutorial: https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_tutorial.html


# ── Configuration ──────────────────────────────────────────────────────────────
INPUT_FOLDER    = "CRC_project/GSE283052"           # folder containing sample subdirectories
INF_AVER_PATH   = "CRC_project/inf_aver.csv"        # path to reference cell type signatures
OUTPUT_FOLDER   = os.path.join(INPUT_FOLDER, "output", "deconvolved")

# Model hyperparameters
N_CELLS_PER_LOCATION = 8      # expected average cell abundance per spot
DETECTION_ALPHA      = 200    # controls normalisation of within-experiment RNA detection variation
MAX_EPOCHS           = 30000  # training epochs
N_POSTERIOR_SAMPLES  = 1000   # samples drawn when exporting posterior
# ──────────────────────────────────────────────────────────────────────────────

input_folder  = INPUT_FOLDER
inf_aver_path = INF_AVER_PATH
output_folder = OUTPUT_FOLDER
os.makedirs(output_folder, exist_ok=True)


# ── Helper: recover raw integer counts for Cell2location ──────────────────────
def recover_raw_counts(adata):
    """
    Cell2location requires raw integer counts in adata.X.
    If adata.X is normalized, try adata.raw → common layer names, in that order.
    Spatial metadata (obsm, uns) is preserved across the conversion.
    """
    import scipy.sparse as sp

    # Sample a small block to check whether values are already integers
    sample = adata.X[:10, :10]
    if sp.issparse(sample):
        sample = sample.toarray()
    if np.all(sample == np.floor(sample)) and sample.min() >= 0:
        print("  [counts] adata.X already contains integer counts. No recovery needed.", flush=True)
        return adata

    print("  [counts] adata.X contains non-integer values — Cell2location needs raw counts.", flush=True)

    # 1. Try adata.raw (most common in processed h5ad files)
    if adata.raw is not None:
        print("  [counts] Recovering from adata.raw ...", flush=True)
        raw = adata.raw.to_adata()
        # adata.raw may cover more genes; keep only those in the processed object
        shared_genes = raw.var_names.intersection(adata.var_names)
        raw = raw[:, shared_genes].copy()
        # Carry over spatial coordinates and images
        for key in adata.obsm.keys():
            raw.obsm[key] = adata.obsm[key]
        raw.uns.update(adata.uns)
        print(f"  [counts] Recovered {raw.n_vars} genes x {raw.n_obs} spots from adata.raw.", flush=True)
        return raw

    # 2. Try common layer names
    for layer in ("counts", "raw_counts", "raw"):
        if layer in adata.layers:
            print(f"  [counts] Recovering from adata.layers['{layer}'] ...", flush=True)
            adata.X = adata.layers[layer]
            return adata

    print("  WARNING: No raw counts source found (adata.raw is None, no 'counts' layer).", flush=True)
    print("           Proceeding with current adata.X — training will likely fail.", flush=True)
    return adata
# ──────────────────────────────────────────────────────────────────────────────

# ── Helper: load Visium from h5ad + spatial/ folder ───────────────────────────
def load_h5ad_visium(sample_path, sample_name):
    """
    Load a Visium sample from an .h5ad file when filtered_feature_bc_matrix.h5
    is not available. Attaches spatial images and scale factors from the
    spatial/ subfolder so the AnnData object is equivalent to what
    sq.read.visium() produces.
    """
    h5ad_files = [f for f in os.listdir(sample_path) if f.endswith(".h5ad")]
    if not h5ad_files:
        raise FileNotFoundError(f"No .h5ad file found in: {sample_path}")

    h5ad_path = os.path.join(sample_path, h5ad_files[0])
    print(f"  [h5ad] Loading from: {os.path.basename(h5ad_path)}", flush=True)
    adata = sc.read_h5ad(h5ad_path)

    # Attach spatial images + scale factors from the spatial/ folder
    spatial_dir = os.path.join(sample_path, "spatial")
    if os.path.isdir(spatial_dir):
        scalefactors_path = os.path.join(spatial_dir, "scalefactors_json.json")
        with open(scalefactors_path) as fh:
            scalefactors = json.load(fh)

        images = {}
        for res in ("hires", "lowres"):
            img_path = os.path.join(spatial_dir, f"tissue_{res}_image.png")
            if os.path.exists(img_path):
                images[res] = mpimg.imread(img_path)

        adata.uns["spatial"] = {
            sample_name: {
                "images": images,
                "scalefactors": scalefactors,
            }
        }
        print(f"  [h5ad] Spatial images attached ({list(images.keys())})", flush=True)
    else:
        print("  [h5ad] Warning: no spatial/ folder found; images not attached.", flush=True)

    return adata
# ──────────────────────────────────────────────────────────────────────────────

# ── Device configuration ───────────────────────────────────────────────────────
USE_GPU = torch.cuda.is_available()
if USE_GPU:
    print(f"GPU detected: {torch.cuda.get_device_name(0)}", flush=True)
    accelerator = "gpu"
else:
    print("No GPU detected, falling back to CPU", flush=True)
    accelerator = "cpu"

# ── Load reference signatures once ────────────────────────────────────────────
print(f"Loading inf_aver from: {inf_aver_path}", flush=True)
inf_aver_all = pd.read_csv(inf_aver_path, index_col=0)

# ── Discover sample directories ────────────────────────────────────────────────
# Accept folders with either a standard 10x .h5 file or an .h5ad file.
def _is_valid_visium_dir(folder_path):
    if os.path.exists(os.path.join(folder_path, "filtered_feature_bc_matrix.h5")):
        return True
    return any(f.endswith(".h5ad") for f in os.listdir(folder_path))

sample_dirs = sorted([
    d for d in os.listdir(input_folder)
    if os.path.isdir(os.path.join(input_folder, d))
    and d != "output"
    and _is_valid_visium_dir(os.path.join(input_folder, d))
])

if not sample_dirs:
    print(f"No valid Visium samples found in {input_folder}. "
          "Each sample folder must contain filtered_feature_bc_matrix.h5 or a .h5ad file.", flush=True)
    sys.exit(1)

print(f"\nFound {len(sample_dirs)} sample(s): {sample_dirs}", flush=True)

# ── Process each sample ────────────────────────────────────────────────────────
for sample_name in sample_dirs:
    sample_path = os.path.join(input_folder, sample_name)
    sample_output = os.path.join(output_folder, sample_name)
    os.makedirs(sample_output, exist_ok=True)

    print(f"\n{'='*60}", flush=True)
    print(f"Processing sample: {sample_name}", flush=True)
    print(f"{'='*60}", flush=True)

    # Step 1: Read Visium data
    print("Step1 read data", flush=True)
    h5_path = os.path.join(sample_path, "filtered_feature_bc_matrix.h5")
    if os.path.exists(h5_path):
        # Standard 10x HDF5 format
        spatial_dir = os.path.join(sample_path, "spatial")
        hires_path  = os.path.join(spatial_dir, "tissue_hires_image.png")
        lowres_path = os.path.join(spatial_dir, "tissue_lowres_image.png")
        if not os.path.exists(hires_path) and os.path.exists(lowres_path):
            shutil.copy(lowres_path, hires_path)
            print("  [info] tissue_hires_image.png missing — copied lowres as proxy", flush=True)
        adata = sq.read.visium(
            path=sample_path,
            counts_file="filtered_feature_bc_matrix.h5",
            library_id=sample_name
        )
    else:
        # Fall back to h5ad
        adata = load_h5ad_visium(sample_path, sample_name)
    # Ensure raw integer counts before saving / modelling
    #adata = recover_raw_counts(adata)

    processed_h5ad = os.path.join(sample_output, "processed.h5ad")
    adata.write_h5ad(processed_h5ad)
    adata_vis = sc.read_h5ad(processed_h5ad)
    print(adata_vis.var.columns.tolist(), flush=True)

    # Step 2: Filter mitochondrial genes
    print("Step2 filter MT genes", flush=True)
    adata_vis.var_names_make_unique()
    adata_vis.var['MT_gene'] = [gene.startswith('MT-') for gene in adata_vis.var_names]
    adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT_gene'].values].X.toarray()
    adata_vis = adata_vis[:, ~adata_vis.var['MT_gene'].values]

    # Step 3: Find shared genes with reference
    print("Step3 find shared genes", flush=True)
    inf_aver = inf_aver_all.copy()
    intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
    adata_vis = adata_vis[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()

    # Step 4: Prepare AnnData for model
    print("Step4 prepare visium for model", flush=True)
    cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key=None)

    # Step 5: Create model
    print("Step5 create model", flush=True)
    mod = cell2location.models.Cell2location(
        adata_vis,
        cell_state_df=inf_aver,
        N_cells_per_location=N_CELLS_PER_LOCATION,
        detection_alpha=DETECTION_ALPHA
    )

    # Step 6: Train model
    print("Step6 train model", flush=True)
    mod.train(
        max_epochs=MAX_EPOCHS,
        batch_size=None,
        train_size=1,
        accelerator=accelerator
    )
    model_dir = os.path.join(sample_output, "model")
    mod.save(model_dir, overwrite=True)
    print(f"Model saved to {model_dir}", flush=True)

    # Step 7: Export posterior
    print("Step7 export posterior", flush=True)
    adata_vis = mod.export_posterior(
        adata_vis, sample_kwargs={'num_samples': N_POSTERIOR_SAMPLES, 'batch_size': mod.adata.n_obs}
    )
    deconvolved_h5ad = os.path.join(sample_output, f"{sample_name}_deconvolved.h5ad")
    adata_vis.write_h5ad(deconvolved_h5ad)
    print(f"Deconvolved h5ad saved to {deconvolved_h5ad}", flush=True)

    # Step 8: Extract and save cell type abundances
    print("Step8 extract cell type abundances", flush=True)
    abundance_df = pd.DataFrame(
        adata_vis.obsm['means_cell_abundance_w_sf'].values,
        index=adata_vis.obs_names,
        columns=adata_vis.uns['mod']['factor_names']
    )
    abundance_df_5 = pd.DataFrame(
        adata_vis.obsm['q05_cell_abundance_w_sf'].values,
        index=adata_vis.obs_names,
        columns=adata_vis.uns['mod']['factor_names']
    )

    abundance_csv = os.path.join(sample_output, f"{sample_name}_abundance.csv")
    abundance_df_5_csv = os.path.join(sample_output, f"{sample_name}_conservative_abundance.csv")
    abundance_df.to_csv(abundance_csv)
    abundance_df_5.to_csv(abundance_df_5_csv)

    print(f"Abundance table shape: {abundance_df.shape}", flush=True)
    print(f"Results saved to {sample_output}", flush=True)

print(f"\nAll samples complete. Results in: {output_folder}", flush=True)
