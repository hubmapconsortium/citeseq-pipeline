cwlVersion: v1.1
class: CommandLineTool
label: HTO quantification
requirements:
  DockerRequirement:
      dockerPull: hubmap/citeseq_analysis:latest
baseCommand: /opt/adt_hto_salmon_wrapper.py

inputs:
  fastq_dir_hto:
    type: Directory
    inputBinding:
      position: 0
      prefix: "--fastq_dir"
  threads:
    type: int
    inputBinding:
      position: 1
      prefix: "--threads"
  hto_index:
    type: Directory
    inputBinding:
      position: 2
      prefix: "--index_dir"
outputs:
  hto_salmon_dir:
    type: Directory
    outputBinding:
      glob: alevin_hto