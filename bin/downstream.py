#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
#from turtle import color
import scanpy as sc
from os.path import exists
import pandas as pd
import muon as mu
from muon import prot as pt
import matplotlib.pyplot as plt
from plot_utils import new_plot
import numpy as np

def main(
	muon_ori: Path
):
	# Read the muon data here
	rna_expr = mu.read(str(muon_ori / "rna"))
	rna_expr.X = rna_expr.layers['spliced']
	adt_expr = mu.read(str(muon_ori / "adt"))
	
	print("Construct muon object with RNA and ADT expression data.")
	mdata_raw = mu.MuData({'rna': rna_expr, 'adt':adt_expr})
	print(mdata_raw)
	
	mdata_rna = mdata_raw.mod['rna']
	mdata_adt = mdata_raw.mod['adt']
	
	## ADT normalization: CLR method
	print("Normalizing ADT data...")
	pt.pp.clr(mdata_adt)

	## RNA normalization
	# For now just use the same method as the original pipeline
	# Do we need to filter the cells and genes? 
	# Can filter out cells that do not pass QC
	print("Normalizing RNA data...")
	sc.pp.normalize_total(mdata_rna, target_sum=1e4)
	sc.pp.log1p(mdata_rna)

	## Downstream analysis for ADT
	print("Performing downstream analysis for ADT...")
	sc.tl.pca(mdata_adt)
	sc.pp.neighbors(mdata_adt)
	for axis in [0, 1]:
		neighbor_counts = np.array((mdata_adt.obsp['distances'] > 0).sum(axis=axis)).flatten()
		mdata_adt.obs.loc[:, f'neighbor_counts_ax{axis}'] = neighbor_counts
		cells_with_neighbors = (neighbor_counts > 0).astype(int)
		mdata_adt.obs.loc[:, f'has_neighbors_ax{axis}'] = cells_with_neighbors
	mdata_adt.obs.loc[:, 'has_neighbors'] = mdata_adt.obs.loc[:, 'has_neighbors_ax0'] * mdata_adt.obs.loc[:, 'has_neighbors_ax1']
	mu.pp.filter_obs(mdata_adt, 'has_neighbors_ax0', lambda x: x > 0)
	sc.tl.leiden(mdata_adt)
	sc.tl.umap(mdata_adt)
	with new_plot():
		sc.pl.umap(mdata_adt, color="leiden", legend_loc="on data")
		plt.savefig("leiden_cluster_adt.pdf")
	
	## Downstream analysis for RNA
	print("Performing downstream analysis for RNA...")
    	
	# find highly variable genes
	sc.pp.highly_variable_genes(mdata_rna, min_mean=0.02, max_mean=4, min_disp=0.5)
	print("Found", np.sum(mdata_rna.var.highly_variable), "highly variable genes.")
    
	# scaling the data
	# save log-normalised counts
	mdata_rna.raw = mdata_rna
	print("Scaling the log-normalised counts to zero mean and unit variance...")
	sc.pp.scale(mdata_rna, max_value=10)
	
	# pca and neighbourhood graph
	print("Performing PCA...")
	sc.tl.pca(mdata_rna)
	print("Constructing neighborhood graph...")
	sc.pp.neighbors(mdata_rna, n_neighbors=10, n_pcs=20)
	
	for axis in [0, 1]:
		neighbor_counts = np.array((mdata_rna.obsp['distances'] > 0).sum(axis=axis)).flatten()
		mdata_rna.obs.loc[:, f'neighbor_counts_ax{axis}'] = neighbor_counts
		cells_with_neighbors = (neighbor_counts > 0).astype(int)
		mdata_rna.obs.loc[:, f'has_neighbors_ax{axis}'] = cells_with_neighbors
	mdata_rna.obs.loc[:, 'has_neighbors'] = mdata_rna.obs.loc[:, 'has_neighbors_ax0'] * mdata_rna.obs.loc[:, 'has_neighbors_ax1']
	mu.pp.filter_obs(mdata_rna, 'has_neighbors_ax0', lambda x: x > 0)
	# clustering
	print("Performing leiden clustering with the computed neighbourhood graph...")
	sc.tl.leiden(mdata_rna, resolution=.75)
	sc.tl.umap(mdata_rna, spread=1., min_dist=.5, random_state=11)
	with new_plot():
		sc.pl.umap(mdata_rna, color="leiden", legend_loc="on data")
		plt.savefig("leiden_cluster_rna.pdf")
	# print(mdata_adt)
	# print(mdata_rna)
	## Multi-omics integration
	print("Performing multi-omics integration...")
	mdata_raw.update()
    # if we filter the cells at the RNA QC step, subset them in the protein modality
	print("Intersecting the common cells...")
	mu.pp.intersect_obs(mdata_raw)
	print(mdata_raw)
	# multiplex clustering
	mu.pp.neighbors(mdata_raw)
	mu.tl.umap(mdata_raw)
	mu.tl.louvain(mdata_raw, resolution=[2, .1], random_state=1)
	mu.tl.leiden(mdata_raw, resolution=[2, .1], random_state=1)
	print(mdata_raw)
	with new_plot():
		mu.pl.embedding(mdata_raw, basis="X_umap", color="louvain")
		plt.savefig("louvain_cluster_combined.pdf")
	with new_plot():
		mu.pl.embedding(mdata_raw, basis="X_umap", color="leiden")
		plt.savefig("leiden_cluster_combined.pdf")

	
	# multi-omics factor analysis
	# generate an interpretable latent space for both RNA and ADT modalities
	# mdata_adt.var["highly_variable"] = True
	# mdata_raw.update()
	# mu.tl.mofa(mdata_raw, outfile="mofa.hdf5", n_factors=30)
	# print(mdata_raw.obsm["X_mofa"].shape)
	# sc.pp.neighbors(mdata_raw, use_rep="X_mofa")
	# sc.tl.umap(mdata_raw, random_state=1)
	# with new_plot():
        	# mu.pl.embedding(mdata_raw, basis="X_umap", color="louvain")
        	# plt.savefig("louvain_cluster_combined.pdf")
	#with new_plot():
        	#mu.pl.embedding(mdata_raw, basis="X_umap", color="leiden")
	        #plt.savefig("leiden_cluster_combined.pdf")

	mdata_raw.write("citeseq_downstream.h5mu")
	
if __name__ == "__main__":
	p = ArgumentParser()
	p.add_argument("--muon_ori", type=Path)

	args = p.parse_args()

	main(
		args.muon_ori,
	)
