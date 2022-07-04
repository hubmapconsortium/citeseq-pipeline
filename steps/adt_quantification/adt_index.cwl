cwlVersion: v1.2
class: CommandLineTool
label: ADT indexing
requirements:
  DockerRequirement:
      dockerPull: hubmap/citeseq_analysis:latest
baseCommand: [salmon, index, --features, -k7, -i, adt_index]

inputs:
  adt_tsv:
    type: File
    inputBinding:
      position: 0
      prefix: -t
outputs:
  adt_index:
    type: Directory
    outputBinding:
      glob: adt_index