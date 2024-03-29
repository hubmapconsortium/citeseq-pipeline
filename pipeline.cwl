#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.2
label: CITE-seq pipeline using Salmon and Alevin (HuBMAP scRNA-seq pipeline)
requirements:
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}
  InlineJavascriptRequirement: {}
  MultipleInputFeatureRequirement: {}
inputs:
  fastq_dir_rna:
    label: "Directory containing RNA-seq FASTQ files"
    type: Directory[]
  fastq_dir_adt:
    label: "Directory containing ADT FASTQ files"
    type: Directory
  fastq_dir_hto:
    label: "Directory containing HTO FASTQ files"
    type: Directory?
  adt_tsv:
    label: "ADT feature barcode"
    type: File
  hto_tsv:
    label: "HTO feature barcode"
    type: File?
  assay:
    label: "scRNA-seq assay"
    type: string
  threads:
    label: "Number of threads for Salmon"
    type: int
    default: 1
  expected_cell_count:
    type: int?
  keep_all_barcodes:
    type: boolean?
  trans_dir:
    label: "Directory of barcode transformation mapping file for feature barcoding protocol use TotalSeq B or C"
    type: Directory?
  trans_filename:
    label: "Filename of barcode transformation mapping file for feature barcoding protocol use TotalSeq B or C"
    type: string?
outputs:
  muon_original_h5mu:
    outputSource: consolidate_counts/muon_dir
    type: File
    label: "Consolidated expression per cell: gene expression, ADT, HTO (optional)"
  muon_processed_h5mu:
    outputSource: downstream_analysis/muon_processed
    type: File
    label: "Processed version of raw expression for each modality"
  mofa_model:
    outputSource: downstream_analysis/mofa_out
    type: File
    label: "Multi-omics factor analysis model"
  rna_embedding_result:
    outputSource: downstream_analysis/rna_embedding
    type: File
    label: "Leiden clustering result on rna modality"
  adt_embedding_result:
    outputSource: downstream_analysis/adt_embedding
    type: File
    label: "Leiden clustering result on adt modality"
  joint_embedding_result:
    outputSource: downstream_analysis/joint_embedding
    type: File
    label: "Leiden clustering result on joint modality"
steps:
  rna_quantification:
    in:
      fastq_dir:
        source: fastq_dir_rna
      assay:
        source: assay
      threads:
        source: threads
      expected_cell_count:
        source: expected_cell_count
      keep_all_barcodes:
        source: keep_all_barcodes
    out:
      - salmon_output
      - count_matrix_h5ad
      - raw_count_matrix
      - genome_build_json
    run: salmon-rnaseq/steps/salmon-quantification.cwl
  adt_quantification:
    in:
      fastq_dir_adt:
        source: fastq_dir_adt
      adt_tsv:
        source: adt_tsv
      threads:
        source: threads
    out: 
      - count_matrix_h5ad_adt
    run: steps/adt_quantification.cwl
  hto_quantification:
    in:
      fastq_dir_hto:
        source: fastq_dir_hto
      hto_tsv:
        source: hto_tsv
      threads:
        source: threads
    out: 
      - count_matrix_h5ad_hto
    when: $(inputs.fastq_dir_hto != null)
    run: steps/hto_quantification.cwl
  consolidate_counts:
    in:
      count_matrix_h5ad_rna:
        source:
          rna_quantification/count_matrix_h5ad
      count_matrix_h5ad_adt:
        source:
          adt_quantification/count_matrix_h5ad_adt
      count_matrix_h5ad_hto:
        source:
          hto_quantification/count_matrix_h5ad_hto
      trans_dir:
        source:
          trans_dir
      rna_salmon_dir:
        source:
          trans_filename
    out: [muon_dir]
    run: steps/consolidate_counts.cwl
  downstream_analysis:
    in:
      muon_dir:
        source:
          consolidate_counts/muon_dir
    out: 
      - muon_processed
      - mofa_out
      - rna_embedding
      - adt_embedding
      - joint_embedding
    run: steps/downstream.cwl
