#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from typing import Sequence

import pandas as pd
import scipy.io
from anndata import AnnData


def build_anndata(X, rows: Sequence[str], cols: Sequence[str], **kwargs) -> AnnData:
    """
    Helper to construct an AnnData object from a data matrix and
    specified rows/columns, with calls to pd.DataFrame for obs and
    var data structures
    """
    return AnnData(
        X=X,
        obs=pd.DataFrame(index=rows),
        var=pd.DataFrame(index=cols),
        **kwargs,
    )


def convert(input_dir: Path) -> AnnData:
    alevin_dir = input_dir / "alevin"

    with open(alevin_dir / "quants_mat_rows.txt") as f:
        cb_names = [line.strip() for line in f]

    with open(alevin_dir / "quants_mat_cols.txt") as f:
        gene_names = [line.strip() for line in f]

    print("Reading count matrix")
    raw_matrix = scipy.io.mmread(alevin_dir / "quants_mat.mtx.gz").tocsr()
    raw_labeled = build_anndata(
        X=raw_matrix,
        rows=cb_names,
        cols=gene_names,
    )

    return raw_labeled


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("--alevin_out_dir", type=Path)

    args = p.parse_args()
    raw = convert(args.alevin_out_dir)

    if "adt" in str(args.alevin_out_dir):
        filename = "raw_expr_adt.h5ad"
    elif "hto" in str(args.alevin_out_dir):
        filename = "raw_expr_hto.h5ad"

    raw.write_h5ad(filename)
    print("h5ad file " + filename + " written.")
