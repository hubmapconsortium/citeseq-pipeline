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
    rna_file: Path,
    adt_file: Path,
    hto_file: Path,
	trans_dir: Path,
    trans_filename: str,
):
   
    rna_expr = sc.read_h5ad(rna_file)
    rna_name = list(rna_expr.obs_names)
    adt_expr = sc.read_h5ad(adt_file)
    adt_name = list(adt_expr.obs_names)
    if exists(hto_file):
        hto_expr = sc.read_h5ad(hto_file)
        hto_name = list(hto_expr.obs_names)
    else:
        print("HTO expression matrix is not found. Will not combine HTO experiment.")
    
	# if the transformation file of given name exist, perform transformation step
    if trans_dir != None:
        trans_file_path = trans_dir / trans_filename
        if exists(trans_file_path):
            print("Performing transformation step of the cellular barcodes of RNA...")
            barcode_dict = generate_barcode_dict(trans_file_path)
            rna_expr.obs.index, count = transform_rna_barcode(barcode_dict, rna_name)
            rna_name = list(rna_expr.obs_names)
            print("A total of", count, "out of", len(rna_name), "RNA cell barcodes were not found in RNA barcode transformation dictionary." )
        else:
            raise ValueError(trans_filename, " is not found under given directory.")
    
    if exists(hto_file):
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
    p.add_argument("--rna_file", type=Path)
    p.add_argument("--adt_file", type=Path)
    p.add_argument("--hto_file", type=Path)
    p.add_argument("--trans_dir", type=Path)
    p.add_argument("--trans_filename", type=str)

    args = p.parse_args()

    main(
        args.rna_file,
        args.adt_file,
        args.hto_file,
		args.trans_dir,
        args.trans_filename
	)

    

    




	