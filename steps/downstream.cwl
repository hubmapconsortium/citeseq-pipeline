cwlVersion: v1.1
class: CommandLineTool
label: Downstream analysis for RNA and ADT 
requirements:
  DockerRequirement:
      dockerPull: hubmap/citeseq_analysis:latest
baseCommand: /opt/downstream.py

inputs:
  muon_original:
    type: File
    inputBinding:
      position: 0
      prefix: "--muon_ori"
outputs:
  muon_processed:
    type: File
    outputBinding:
      glob: "citeseq_downstream.h5mu"
  