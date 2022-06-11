#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
import scanpy as sc
from os.path import exists
import pandas as pd
import muon as mu


def generate_barcode_dict(directory: Path):
	colnames = ['transformed','original']
	trans_info = pd.read_csv(directory, sep='\t', header=None, names=colnames)
	barcode_dict = dict(zip(trans_info.transformed, trans_info.original))
	print("Generated ", len(barcode_dict), " key-value pairs for RNA barcode transformation.")
	return barcode_dict

def transform_rna_barcode(barcode_dict: Path, rna_name: set):
	rna_name = list(rna_name)
	count = 0
	for i in range(len(rna_name)):
		barcode = rna_name[i]
		if barcode not in barcode_dict:
			count += 1
			continue
		else:
			rna_name[i] = barcode_dict[barcode]
	return set(rna_name), count

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

    rna_name = set(rna_expr.obs_names)
    adt_name = set(adt_expr.obs_names)
    hto_name = set(hto_expr.obs_names)

	# if the '3M-february-2018.txt' file exist, perform transformation step
    trans_file_path = str(trans_dir) + "/3M-february-2018.txt"

    if exists(trans_file_path):
        print("Find 3M-february-2018.txt under given directory.")
        print("Performing transformation step of the cellular barcodes of RNA...")
        barcode_dict = generate_barcode_dict(trans_file_path)
        rna_name, count = transform_rna_barcode(barcode_dict, rna_name)
        print("A total of", count, "out of", len(rna_name), "RNA cell barcodes were not found in RNA barcode transformation dictionary." )

    common_cells = list(adt_name & hto_name & rna_name)
    print("There are", len(common_cells), "common cells in RNA, ADT and HTO experiments.")
    
    mdata = mu.MuData({'rna': rna_expr[common_cells, :], 'adt': adt_expr[common_cells, :], 'hto': hto_expr[common_cells, :]})
    print("Saving MuData filtered by", len(common_cells), "common cells in RNA, ADT and HTO experiments...")
    mdata.write("mudata.h5mu")


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("--out_dir", type=Path)
    p.add_argument("--trans_dir", type=Path)

    args = p.parse_args()

    main(
        args.out_dir,
		args.trans_dir,
	)

    

    




	