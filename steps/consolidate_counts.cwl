cwlVersion: v1.1
class: CommandLineTool
label: Consolidate RNA, ADT and HTO counts
requirements:
  DockerRequirement:
      dockerPull: hubmap/citeseq_analysis:latest
baseCommand: /opt/consolidate_counts.py

inputs:
  count_matrix_h5ad_rna:
    type: File
    inputBinding:
      position: 0
      prefix: "--rna_file"
  count_matrix_h5ad_adt:
    type: File
    inputBinding:
      position: 1
      prefix: "--adt_file"
  count_matrix_h5ad_hto:
    type: File
    inputBinding:
      position: 2
      prefix: "--hto_file"
  transformation_dir:
    type: Directory?
    inputBinding:
      position: 3
      prefix: "--trans_dir"
  transformation_filename:
    type: string?
    inputBinding:
      position: 4
      prefix: "--trans_filename"
outputs:
  muon_dir:
    type: File
    outputBinding:
      glob: "mudata_raw.h5mu"