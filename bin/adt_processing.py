#!/usr/bin/env python3
from argparse import ArgumentParser
# from msilib.schema import File
from pathlib import Path
import scanpy as sc
# from os.path import exists
import pandas as pd
import muon as mu
from muon import prot as pt
import matplotlib.pyplot as plt

def main(
	mudata: Path
):
	adata = mu.read(str(mudata))
	prot = adata['adt']

	# CLR normalization
	pt.pp.clr(prot)

	# sc.tl.pca(prot)





if __name__ == "__main__":
	p = ArgumentParser()
	p.add_argument("--mudata", type=Path)

	args = p.parse_args()

	main(
		args.mudata,
	)