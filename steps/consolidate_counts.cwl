cwlVersion: v1.1
class: CommandLineTool
label: Consolidate RNA, ADT and HTO counts
requirements:
  DockerRequirement:
      dockerPull: hubmap/citeseq_analysis:latest
baseCommand: /opt/consolidate_counts.py

inputs:
  raw_expr_h5ad_dir:
    type: Directory
    inputBinding:
      position: 0
      prefix: "--out_dir"
  transformation_dir:
    type: Directory?
    inputBinding:
      position: 1
      prefix: "--trans_dir"
  transformation_filename:
    type: string?
    inputBinding:
      position: 2
      prefix: "--trans_filename"
outputs:
  muon_original:
    type: File
    outputBinding:
      glob: "mudata_raw.h5mu"