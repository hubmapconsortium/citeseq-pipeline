#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.0
label: CITE-seq pipeline using Salmon and Alevin (HuBMAP scRNA-seq pipeline)
requirements:
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}
inputs:
  fastq_dir_rna:
    label: "Directory containing RNA-seq FASTQ files"
    type: Directory
  fastq_dir_adt:
    label: "Directory containing ADT FASTQ files"
    type: Directory
  fastq_dir_hto:
    label: "Directory containing HTO FASTQ files"
    type: Directory
  assay:
    label: "scRNA-seq assay"
    type: string
  threads:
    label: "Number of threads for Salmon"
    type: int
    default: 1
  expected_cell_count:
    type: int?
outputs:
  count_matrix_h5mu:
    outputSource: consolidate_counts/expr_h5mu
    type: File
    label: "Consolidated expression per cell: gene expression, ADT, HTO"
steps:
  prep_adt_hto_sequences:
    in:
      fastq_dir_adt:
        source: fastq_dir_adt
      fastq_dir_hto:
        source: fastq_dir_hto
    out: [adt_tsv, hto_tsv]
    run: steps/prep_adt_hto_sequences.cwl
  adt_index:
    in:
      adt_tsv:
        source: prep_adt_hto_sequences/adt_tsv
      threads:
        source: threads
    out: [adt_index]
    run: steps/adt_index.cwl
  hto_index:
    in:
      hto_tsv:
        source: prep_adt_hto_sequences/hto_tsv
      threads:
        source: threads
    out: [hto_index]
    run: steps/hto_index.cwl
  adt_quant:
    in:
      adt_index:
        source: adt_index/adt_index
      threads:
        source: threads
    out: [adt_salmon_dir]
    run: steps/adt_quant.cwl
  hto_quant:
    in:
      hto_index:
        source: hto_index/hto_index
      threads:
        source: threads
    out: [hto_salmon_dir]
    run: steps/hto_quant.cwl
  consolidate_counts:
    in:
      adt_salmon_dir:
        source:
          adt_quant/adt_salmon_dir
      hto_salmon_dir:
        source:
          hto_quant/hto_salmon_dir
      rna_salmon_dir:
        source:
          rna_quant/rna_salmon_dir
    out:
      [expr_h5mu]
    run: steps/consolidate_counts.cwl
