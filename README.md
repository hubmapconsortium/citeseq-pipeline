# HBM-CITEseq: HuBMAP CITE-seq pipeline

## Overview
HBM-CITEseq is an extension for HuBMAP scRNA-seq pipeline (https://github.com/hubmapconsortium/salmon-rnaseq) to process CITE-seq data. It is built on Alevin, Scanpy and Muon, and is implemented as a CWL workflow wrapping commamd-line tools encapsulated in Docker containers.

## Requirements
Running the pipeline requires a CWL workflow execution engine and container runtime; we recommend Docker and the `cwltool` reference implementation. `cwltool` is written in Python and can be installed into a sufficiently recent Python environment with `pip install cwltool`. Afterward, clone this repository, check out a tag, and invoke the pipeline as:

```
cwltool pipeline.cwl --fastq_dir_rna RNA_FASTQ_DIR --fastq_dir_adt ADT_FASTQ_DIR --adt_tsv ADT_BARCODE --assay ASSAY --threads THREADS
```

Supplementary info for input:  
+ HTO data is an optional input for HBM-CITEseq. Extend the command to input HTO raw data and barcode information:
```
--fastq_dir_hto HTO_FASTQ_DIR --hto_tsv HTO_BARCODE
```
+ In order to quantify ADT and HTO data, we need to index the feature barcodes. The indexing provided by salmon should base on a tab separated file (tsv format), where each line represents one reference, id and the reference sequence are separated by the tab. Example `adt.tsv` and `hto.tsv` files can be found under `example_data` folder.

+ If the feature barcoding protocol use TotalSeq B or C, the cellular barcodes of RNA and the feature barcodes of ADT/HTO might not exactly be the same. Transformation needs to be performed based on a mapping file.  
To download the mapping file: https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/translation/3M-february-2018.txt.gz  
More information from 10XGENOMICS: https://kb.10xgenomics.com/hc/en-us/articles/360031133451-Why-is-there-a-discrepancy-in-the-3M-february-2018-txt-barcode-whitelist-  
To perform the mapping, extend the original command by:
```
--trans_dir MAPPING_FILE_DIR --trans_filename MAPPING_FILE_NAME
```

## Possible issue
For the HTO quantification, all entries in the output HTO expression matrix might be zero. This is due to the malformed MTX file in HTO quantification output of Salmon. We've reported the issue here: https://github.com/COMBINE-lab/salmon/issues/791

## Environments


