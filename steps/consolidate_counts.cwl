cwlVersion: v1.1
class: CommandLineTool
label: Consolidate RNA, ADT and HTO counts
requirements:
  DockerRequirement:
      dockerPull: hubmap/citeseq_analysis:latest
baseCommand: /opt/consolidate_counts.py

inputs:
  rawexpr_dir:
    type: Directory
    inputBinding:
      position: 0
      prefix: "--out_dir"
  transformation_dir:
    type: Directory
    inputBinding:
      position: 1
      prefix: "--trans_dir"
outputs:
  muon_out_dir:
    type: File
    outputBinding:
      glob: mudata.h5mu