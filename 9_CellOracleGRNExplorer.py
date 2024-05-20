#!/usr/bin/env python
# coding: utf-8

# ## Import libraries

import os
import sys
import copy
import glob
import time
import shutil

import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import anndata as ad

from tqdm.auto import tqdm

import celloracle as co
from celloracle.applications import Pseudotime_calculator
co.__version__


# ## Load scRNA-seq data

# Load the raw data and metadata from TSV files
raw_data_file_path = "/media/cdn-bc/RAID/StudentProjects/Alaric/outputs/raw_data.tsv"
metadata_file_path = "/media/cdn-bc/RAID/StudentProjects/Alaric/outputs/metadata.tsv"

raw_data_df = pd.read_csv(raw_data_file_path, sep="\t", index_col=0)
metadata_df = pd.read_csv(metadata_file_path, sep="\t", index_col=0)

# Create an AnnData object with processed data
adata = ad.AnnData(X=raw_data_df, obs=metadata_df)

# Copy the raw data into a layer called "raw_count"
adata.layers["raw_count"] = adata.X.copy()

# Dimensionalty reduction
# Perform PCA on the dataset
sc.tl.pca(adata)

# Perform t-SNE
sc.tl.tsne(adata, n_pcs=50)  # Using the first 50 principal components for t-SNE

# Save result as a parquet
df = pd.read_table('/media/cdn-bc/RAID/StudentProjects/Alaric/outputs/GRN_MATRIX.tsv')
df.to_parquet("/media/cdn-bc/RAID/StudentProjects/Alaric/outputs/GRN_MATRIX.parquet")

# Load TF info 
base_GRN = pd.read_parquet('/media/cdn-bc/RAID/StudentProjects/Alaric/outputs/GRN_MATRIX.parquet')

# ## Make Oracle object

# Instantiate Oracle object
oracle = co.Oracle()

# Check data in anndata
print("Metadata columns :", list(adata.obs.columns))
print("Dimensional reduction: ", list(adata.obsm.keys()))

# Unscaled mRNA count for the input of Oracle object.
adata.X = adata.layers["raw_count"].copy()

# Assign a placeholder value to NaN cells
# Add 'Unknown' as a new category to the 'cell.subtype' column
# adata.obs['cell.subtype'] = adata.obs['cell.subtype'].cat.add_categories('Unknown') #Error running 
# Fill NaN values with 'Unknown'
adata.obs['cell.subtype'].fillna('Unknown', inplace=True)

# Save the Anndata object to a file
adata.write("/media/cdn-bc/RAID/StudentProjects/Alaric/outputs/adata_celloracle.h5ad")
# Load the Anndata object from the file
adata = ad.read_h5ad("/media/cdn-bc/RAID/StudentProjects/Alaric/outputs/adata_celloracle.h5ad")

# Instantiate Oracle object.
oracle.import_anndata_as_raw_count(adata=adata,cluster_column_name="cell.subtype",embedding_name="X_tsne")

# Load base-GRN data into oracle object
oracle.import_TF_data(TF_info_matrix=base_GRN)

# KNN imputation
# CellOracle uses the same strategy as velocyto for visualizing cell transitions. This process requires KNN imputation in advance.
# 
# For the KNN imputation, we first need to calculate and select PCs.
# 
# ### PCA
# Perform PCA
oracle.perform_PCA()

# Select important PCs
plt.plot(np.cumsum(oracle.pca.explained_variance_ratio_)[:100])
n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
plt.axvline(n_comps, c="k")
plt.show()
print(n_comps)
n_comps = min(n_comps, 50)

# ### KNN imputation
# Estimate the optimal number of nearest neighbors for KNN imputation.
n_cell = oracle.adata.shape[0]
print(f"cell number is :{n_cell}")

k = int(0.025*n_cell)
print(f"Auto-selected k is :{k}") # 125

oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8, b_maxl=k*4, n_jobs=4)

# Save oracle object.
oracle.to_hdf5("/media/cdn-bc/RAID/StudentProjects/Alaric/outputs/CellOracle.celloracle.oracle")
# Loading oracle object
#oracle = co.load_hdf5("/media/cdn-bc/RAID/StudentProjects/Alaric/outputs/CellOracle.celloracle.oracle")


# ## GRN calculation
# The next step constructs a cluster-specific GRN for all clusters.
# 
# - You can calculate GRNs with the `get_links` function, and it will return the results as a `Links` object.
# The `Links` object stores the inferred GRNs and the corresponding metadata. Most network structure analysis is performed with the `Links` object.
# - A GRN will be calculated for each cluster/sub-group. In the example below, we construct GRN for each unit of the "SCT_snn_res_1" clustering.

# Check clustering data
sc.pl.tsne(oracle.adata, color="cell.type", palette=sns.color_palette("hls", 128),legend_loc='on data')

# Get GRNs
links = oracle.get_links(cluster_name_for_GRN_unit="cell.subtype", alpha=10, verbose_level=10)

# ## Network preprocessing
# ### Filter network edges 
# 
# Using the base GRN, CellOracle constructs the GRN models as a list of directed edges between a TF and its target genes.
# We need to remove the weak edges or insignificant edges before doing network structure analysis.
# 
# We filter the network edges as follows.
# 
# 1. Remove uncertain network edges based on the p-value.
#  
# 2. Remove weak network edge. In this tutorial, we keep the top 2000 edges ranked by edge strength.
# 
# 
# The raw network data is stored in the `links_dict` attribute, while the filtered network data is stored in the `filtered_links` attribute. 

links.filter_links(p=0.001, weight="coef_abs", threshold_number=2000)

# ## Degree distribution 
# In the first step, we examine the network degree distribution.
# 
# Network degree, which is the number of edges for each node, is one of the important metrics used to investigate the network structure (https://en.wikipedia.org/wiki/Degree_distribution).
# 
# Please keep in mind that the degree distribution may change depending on the filtering threshold.

# ## Calculate netowork score
# 
# Next, we calculate several network scores.
# Although previous version of celloracle, `links.get_score()` required R packages for the netowrk score calculation, the function was replaced with new function, `links.get_network_scores()`, which does NOT require any R package. This new function can be available with celloracle >= 0.10.0.
# The old function, `links.get_score()` is still available, but it will be removed in the future version.

# Calculate network scores. 
links.get_network_score()

# The score is stored as a attribute `merged_score`.
links.merged_score.head()
links.merged_score

# Save Links object.
links.to_hdf5(file_path="/media/cdn-bc/RAID/StudentProjects/Alaric/outputs/cell_oraclelinks_scaleup.celloracle.links")
# Load Links object.
links = co.load_hdf5(file_path="/media/cdn-bc/RAID/StudentProjects/Alaric/outputs/cell_oraclelinks_scaleup.celloracle.links")

## In silico gene perturbation with GRNs
# Load file.
oracle = co.load_hdf5("/media/cdn-bc/RAID/StudentProjects/Alaric/outputs/CellOracle.celloracle.oracle")
oracle

# ## Load inferred GRNs
# Now, we will use these GRNs for the perturbation simulations. First, we will import the GRNs from the `Links` object.
# Attention!! Please use the function below when you use your data.
links = co.load_hdf5(file_path="/media/cdn-bc/RAID/StudentProjects/Alaric/outputs/cell_oraclelinks_scaleup.celloracle.links")

# # Make predictive models for simulation
# Here, we will need to fit the ridge regression models again. This process will take less time than the GRN inference in the previous notebook, because we are using the filtered GRN models.

links.filter_links()
oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=True)

# # In silico TF perturbation analysis
# Next, we will simulate the TF perturbation effects on cell identity to investigate its potential functions and regulatory mechanisms. Please see the CellOracle paper for more details on scientific rationale.
# ## Check gene expression pattern.

a = oracle.adata.layers['imputed_count']

# Creating a DataFrame from the layer
df = pd.DataFrame(data=a, index=oracle.adata.obs_names, columns=oracle.adata.var_names)

# Saving the DataFrame to a CSV file
df.to_csv('/media/cdn-bc/RAID/StudentProjects/Alaric/outputs/imputed_scaled_count.csv')


# Calculate future gene expression after perturbation.
# 
# - You can use any gene expression value in the in silico perturbations, but please avoid extremely high values that are far from the natural gene expression range. The upper limit allowed is twice the maximum gene expression.
# 
# To simulate gene KO, we will set this gene expression to 0.
genes_to_ko = ["daf-7", "ins-6", "tph-1", "daf-9"]
for goi in genes_to_ko:
    
    tmporacle = oracle.copy()
    # oracle.simulate_shift() performs the KO and alters internal state
    tmporacle.simulate_shift(perturb_condition={goi: 0.0}, n_propagation=3)
    
    
    # Extract the dataframe from the "simulated_count_layer"
    altered_expression_df =tmporacle.adata.to_df(layer="simulated_count")
    
    altered_expression_df['KO_gene'] = goi  # Tagging the dataframe with the KO gene
    
    # Save to a CSV file
    output_path = f'/media/cdn-bc/RAID/StudentProjects/Alaric/outputs/KO_{goi}_altered_expression.csv'
    altered_expression_df.to_csv(output_path)
    
    print(f"Saved altered expression data for {goi} knockout to {output_path}")


# Prepare the perturbation condition for multiple gene KOs
perturb_conditions = {gene: 0.0 for gene in genes_to_ko}

# Perform the combined knockout simulation
oracle.simulate_shift(perturb_condition=perturb_conditions, n_propagation=3)

# Extract the dataframe from the "simulated_count" layer for the combined KO
combined_altered_expression_df = oracle.adata.to_df(layer="simulated_count")

# You might want to tag the dataframe with a combined identifier if needed
combined_altered_expression_df['KO_genes'] = 'combined_KO'

# Save the combined KO altered expression data to a CSV file
combined_altered_expression_df.to_csv('/media/cdn-bc/RAID/StudentProjects/Alaric/outputs/combined_KO_altered_expression.csv')


