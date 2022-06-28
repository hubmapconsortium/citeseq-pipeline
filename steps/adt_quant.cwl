cwlVersion: v1.1
class: CommandLineTool
label: ADT quantification
requirements:
  DockerRequirement:
      dockerPull: hubmap/citeseq_analysis:latest
baseCommand: /opt/adt_hto_salmon_wrapper.py

inputs:
  fastq_dir_adt:
    type: Directory
    inputBinding:
      position: 0
      prefix: "--fastq_dir"
  threads:
    type: int
    inputBinding:
      position: 1
      prefix: "--threads"
  adt_index:
    type: Directory
    inputBinding:
      position: 2
      prefix: "--index_dir"
outputs:
  adt_salmon_dir:
    type: Directory
    outputBinding:
      glob: alevin_adt