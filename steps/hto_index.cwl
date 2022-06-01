cwlVersion: v1.1
class: CommandLineTool
label: HTO indexing
requirements:
  DockerRequirement:
      dockerPull: hubmap/citeseq_analysis:latest
baseCommand: [salmon, index, --features, -k7]

inputs:
  hto_tsv:
    type: File
    inputBinding:
      position: 0
      prefix: -t
  out_dir:
    type: string
    inputBinding:
      position: 1
      prefix: -i
outputs:
  hto_index_dir:
    type: Directory
    outputBinding:
      glob: hto_index