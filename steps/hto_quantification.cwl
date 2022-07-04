#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.2
label: HTO quantification, FASTQ -> H5AD count matrix
inputs:
  fastq_dir_hto:
    label: "Directory containing FASTQ files"
    type: Directory
  hto_tsv:
    label: "HTO feature barcode"
    type: File
  threads:
    label: "Number of threads for Salmon"
    type: int
    default: 1
outputs:
  count_matrix_h5ad_hto:
    outputSource: hto_conversion/raw_expr_h5ad_hto
    type: File
    label: "Unfiltered count matrix from Alevin, converted to H5AD"

steps:
  hto_index:
    in:
      hto_tsv:
        source: hto_tsv
    out: [hto_index]
    run: hto_quantification/hto_index.cwl
  hto_quant:
    in:
      hto_index:
        source: hto_index/hto_index
      fastq_dir_hto:
        source: fastq_dir_hto
      threads:
        source: threads
    out: [hto_salmon_dir]
    run: hto_quantification/hto_quant.cwl
  hto_conversion:
    in:
      hto_salmon_dir:
        source: hto_quant/hto_salmon_dir
    out: [raw_expr_h5ad_hto]
    run: hto_quantification/hto_conversion.cwl
  