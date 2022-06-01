cwlVersion: v1.1
class: CommandLineTool
label: ADT output conversion
requirements:
  DockerRequirement:
      dockerPull: hubmap/citeseq_analysis:latest
baseCommand: /opt/alevin_to_anndata.py

inputs:
  adt_out_dir:
    type: Directory
    inputBinding:
      position: 0
      prefix: --alevin_out_dir
  # adt_name:
  #   type: string
  #   inputBinding:
  #     position: 1
  #     prefix: --name
outputs:
  raw_expr_h5ad:
    type: File
    outputBinding:
      glob: raw_expr.h5ad