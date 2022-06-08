#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
import scanpy as sc
# import anndata
import muon as mu


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("--out_dir", type=Path)

    args = p.parse_args()

    rna_path = str(args.out_dir) + "/raw_expr_rna.h5ad"
    adt_path = str(args.out_dir) + "/raw_expr_adt.h5ad"
    hto_path = str(args.out_dir) + "/raw_expr_hto.h5ad"
    
    rna_expr = sc.read_h5ad(rna_path)
    adt_expr = sc.read_h5ad(adt_path)
    hto_expr = sc.read_h5ad(hto_path)

    print(rna_expr)
    print(adt_expr)
    print(hto_expr)

    # print(rna_expr.obs_names)
    # print(adt_expr.obs_names)
    # print(hto_expr.obs_names)
    # print(rna_expr.var_names)
    # print(adt_expr.var_names)
    # print(hto_expr.var_names)
    # print(rna_expr.shape, adt_expr.shape, hto_expr.shape)
    mdata = mu.MuData({'rna': rna_expr, 'adt': adt_expr, 'hto': hto_expr})
    print(mdata)

    




	