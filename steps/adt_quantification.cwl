#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.2
label: ADT quantification, FASTQ -> H5AD count matrix
inputs:
  fastq_dir_adt:
    label: "Directory containing FASTQ files"
    type: Directory
  adt_tsv:
    label: "ADT feature barcode"
    type: File
  threads:
    label: "Number of threads for Salmon"
    type: int
    default: 1
outputs:
  count_matrix_h5ad_adt:
    outputSource: adt_conversion/raw_expr_h5ad_adt
    type: File
    label: "Unfiltered count matrix from Alevin, converted to H5AD"

steps:
  adt_index:
    in:
      adt_tsv:
        source: adt_tsv
    out: [adt_index]
    run: adt_quantification/adt_index.cwl
  adt_quant:
    in:
      adt_index:
        source: adt_index/adt_index
      fastq_dir_adt:
        source: fastq_dir_adt
      threads:
        source: threads
    out: [adt_salmon_dir]
    run: adt_quantification/adt_quant.cwl
  adt_conversion:
    in:
      adt_salmon_dir:
        source: adt_quant/adt_salmon_dir
    out: [raw_expr_h5ad_adt]
    run: adt_quantification/adt_conversion.cwl
  