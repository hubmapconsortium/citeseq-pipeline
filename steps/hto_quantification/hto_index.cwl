cwlVersion: v1.2
class: CommandLineTool
label: HTO indexing
requirements:
  DockerRequirement:
      dockerPull: hubmap/citeseq_analysis:latest
baseCommand: [salmon, index, --features, -k7, -i, hto_index]

inputs:
  hto_tsv:
    type: File
    inputBinding:
      position: 0
      prefix: -t
outputs:
  hto_index:
    type: Directory
    outputBinding:
      glob: hto_index