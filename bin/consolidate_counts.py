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
    trans_filename: str,
):
    rna_path = out_dir / "raw_expr_rna.h5ad"
    adt_path = out_dir / "raw_expr_adt.h5ad"
    hto_path = out_dir / "raw_expr_hto.h5ad"
    
    rna_expr = sc.read_h5ad(rna_path)
    rna_name = list(rna_expr.obs_names)
    adt_expr = sc.read_h5ad(adt_path)
    adt_name = list(adt_expr.obs_names)
    if exists(hto_path):
        hto_expr = sc.read_h5ad(hto_path)
        hto_name = list(hto_expr.obs_names)
    else:
        print("HTO expression matrix is not found. Will not combine HTO experiment.")
    
	# if the transformation file of given name exist, perform transformation step
    trans_file_path = trans_dir / trans_filename
    if exists(trans_file_path):
        print("Find ", trans_filename, " under given directory.")
        print("Performing transformation step of the cellular barcodes of RNA...")
        barcode_dict = generate_barcode_dict(trans_file_path)
        rna_expr.obs.index, count = transform_rna_barcode(barcode_dict, rna_name)
        rna_name = list(rna_expr.obs_names)
        print("A total of", count, "out of", len(rna_name), "RNA cell barcodes were not found in RNA barcode transformation dictionary." )
    
    if exists(hto_path):
        common_cells = list(set(adt_name) & set(hto_name) & set(rna_name))
        print("There are", len(common_cells), "common cells in RNA, ADT and HTO experiments.")
        mdata = mu.MuData({'rna': rna_expr[common_cells, :], 'adt': adt_expr[common_cells, :], 'hto': hto_expr[common_cells, :]})
        print("Saving MuData filtered by", len(common_cells), "common cells in RNA, ADT and HTO experiments...")
    else:
        common_cells = list(set(adt_name) & set(rna_name))
        print("There are", len(common_cells), "common cells in RNA and ADT experiments.")
        mdata = mu.MuData({'rna': rna_expr[common_cells, :], 'adt': adt_expr[common_cells, :]})
        print("Saving MuData filtered by", len(common_cells), "common cells in RNA and ADT experiments...")
    
    mdata.write("mudata_raw.h5mu")

if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("--out_dir", type=Path)
    p.add_argument("--trans_dir", type=Path)
    p.add_argument("--trans_filename", type=str)

    args = p.parse_args()

    main(
        args.out_dir,
		args.trans_dir,
        args.trans_filename
	)

    

    




	