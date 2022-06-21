#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
import scanpy as sc
import os
from os.path import exists
import pandas as pd
import muon as mu


def generate_barcode_dict(directory: Path):
	colnames = ['transformed','original']
	trans_info = pd.read_csv(directory, sep='\t', header=None, names=colnames)
	barcode_dict = dict(zip(trans_info.transformed, trans_info.original))
	print("Generated ", len(barcode_dict), " key-value pairs for RNA barcode transformation.")
	return barcode_dict

def transform_rna_barcode(barcode_dict: Path, rna_name: list):
	count = 0
	for i in range(len(rna_name)):
		barcode = rna_name[i]
		if barcode not in barcode_dict:
			count += 1
			continue
		else:
			rna_name[i] = barcode_dict[barcode]
	return rna_name, count

def main(
    out_dir: Path,
	trans_dir: Path,
):
    rna_path = str(out_dir) + "/raw_expr_rna.h5ad"
    adt_path = str(out_dir) + "/raw_expr_adt.h5ad"
    hto_path = str(out_dir) + "/raw_expr_hto.h5ad"
    
    rna_expr = sc.read_h5ad(rna_path)
    adt_expr = sc.read_h5ad(adt_path)
    hto_expr = sc.read_h5ad(hto_path)

    rna_name = list(rna_expr.obs_names)
    adt_name = list(adt_expr.obs_names)
    hto_name = list(hto_expr.obs_names)

	# if the '3M-february-2018.txt' file exist, perform transformation step
    trans_file_path = str(trans_dir) + "/3M-february-2018.txt"

    if exists(trans_file_path):
        print("Find 3M-february-2018.txt under given directory.")
        print("Performing transformation step of the cellular barcodes of RNA...")
        barcode_dict = generate_barcode_dict(trans_file_path)
        rna_expr.obs.index, count = transform_rna_barcode(barcode_dict, rna_name)
        rna_name = list(rna_expr.obs_names)
        print("A total of", count, "out of", len(rna_name), "RNA cell barcodes were not found in RNA barcode transformation dictionary." )
    common_cells = list(set(adt_name) & set(hto_name) & set(rna_name))
    print("There are", len(common_cells), "common cells in RNA, ADT and HTO experiments.")
    
    # Write h5ad file for RNA, ADT and HTO separately
    rna_expr = rna_expr[common_cells, :]
    adt_expr = adt_expr[common_cells, :]
    hto_expr = hto_expr[common_cells, :]
    print("Saving h5ad files filtered by", len(common_cells), "common cells in RNA, ADT and HTO experiments...")
    os.mkdir("anndata_out")
    rna_expr.write_h5ad("anndata_out/rna_filtered.h5ad")
    adt_expr.write_h5ad("anndata_out/adt_filtered.h5ad")
    hto_expr.write_h5ad("anndata_out/hto_filtered.h5ad")

if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("--out_dir", type=Path)
    p.add_argument("--trans_dir", type=Path)

    args = p.parse_args()

    main(
        args.out_dir,
		args.trans_dir,
	)

    

    




	