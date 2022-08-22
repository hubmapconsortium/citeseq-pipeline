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
	muon_dir: Path
):
	rna_expr = mu.read("mudata_raw.h5mu/rna")
	rna_expr.X = rna_expr.layers['spliced']
	adt_expr = mu.read("mudata_raw.h5mu/adt")
	print("Construct muon object with RNA and ADT expression data.")
	mdata_raw = mu.MuData({'rna': rna_expr, 'adt':adt_expr})
	print(mdata_raw)
	
	mdata_rna = mdata_raw.mod['rna']
	mdata_adt = mdata_raw.mod['adt']
	
	## ADT normalization: CLR method
	print("Normalizing ADT data...")
	pt.pp.clr(mdata_adt)

	## RNA normalization
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
	mdata_adt.obs.loc[:, 'has_neighbors'] = mdata_adt.obs.loc[:, 'has_neighbors_ax1']
	mu.pp.filter_obs(mdata_adt, 'has_neighbors', lambda x: x > 0)

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
	sc.pp.neighbors(mdata_rna)

	for axis in [0, 1]:
		neighbor_counts = np.array((mdata_rna.obsp['distances'] > 0).sum(axis=axis)).flatten()
		mdata_rna.obs.loc[:, f'neighbor_counts_ax{axis}'] = neighbor_counts
		cells_with_neighbors = (neighbor_counts > 0).astype(int)
		mdata_rna.obs.loc[:, f'has_neighbors_ax{axis}'] = cells_with_neighbors
	mdata_rna.obs.loc[:, 'has_neighbors'] = mdata_rna.obs.loc[:, 'has_neighbors_ax1']
	mu.pp.filter_obs(mdata_rna, 'has_neighbors', lambda x: x > 0)

	# clustering
	print("Performing leiden clustering with the computed neighbourhood graph...")
	sc.tl.leiden(mdata_rna)
	sc.tl.umap(mdata_rna)
	with new_plot():
		sc.pl.umap(mdata_rna, color="leiden", legend_loc="on data")
		plt.savefig("leiden_cluster_rna.pdf")

	## Multi-omics integration
	mdata_raw.update()
    # # if we filter the cells at the RNA QC step, subset them in the protein modality
	# mu.pp.intersect_obs(mdata_raw)
	print(mdata_raw)
	mdata_raw.write("citeseq_normalized.h5mu")
	## Multi-omics factor analysis
	mdata_adt.var["highly_variable"] = True
	mdata_raw.update()
	mu.tl.mofa(mdata_raw, outfile="citeseq_mofa.hdf5", n_factors = 30)

	# multiplex clustering
	sc.pp.neighbors(mdata_raw['rna'])
	sc.pp.neighbors(mdata_raw['adt'])
	print(mdata_raw)
	sc.pp.neighbors(mdata_raw, use_rep="X_mofa", key_added="mofa")
	sc.tl.umap(mdata_raw, neighbors_key='mofa')
	sc.tl.leiden(mdata_raw, resolution=1.0, neighbors_key='mofa', key_added='leiden_wnn')
	# print(mdata_raw)
	with new_plot():
		sc.pl.umap(mdata_raw, color='leiden_wnn', legend_loc='on data')
		plt.savefig("leiden_cluster_combined.pdf")

	# mdata_raw.write("citeseq_downstream.h5mu")
	
if __name__ == "__main__":
	p = ArgumentParser()
	p.add_argument("--muon_dir", type=Path)

	args = p.parse_args()

	main(
		args.muon_dir,
	)
