cwlVersion: v1.1
class: CommandLineTool
label: Extract ADT and HTO sequences for indexing
requirements:
  DockerRequirement:
      dockerPull: hubmap/citeseq_analysis:latest
  DockerGpuRequirement: {}
baseCommand: "/opt/prep_adt_hto_sequences.py"

inputs:
  adt_hto_metadata:
    type: File
    inputBinding:
      position: 0
outputs:
  adt_tsv:
    type: File
    outputBinding:
      glob: "adt.tsv"
  hto_tsv:
    type: File
    outputBinding:
      glob: "hto.tsv"
