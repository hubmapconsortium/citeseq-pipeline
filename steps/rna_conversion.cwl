cwlVersion: v1.1
class: CommandLineTool
label: RNA output conversion
requirements:
  DockerRequirement:
      dockerPull: hubmap/citeseq_analysis:latest
baseCommand: /opt/alevin_to_anndata.py

inputs:
  rna_salmon_dir:
    type: Directory
    inputBinding:
      position: 0
      prefix: --alevin_out_dir
outputs:
  raw_expr_h5ad_dir:
    type: Directory
    outputBinding:
      glob: raw_expr_out