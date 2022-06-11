cwlVersion: v1.1
class: CommandLineTool
label: HTO quantification
requirements:
  DockerRequirement:
      dockerPull: hubmap/citeseq_analysis:latest
baseCommand: /opt/adt_hto_salmon_wrapper.py

inputs:
  fastq_dir:
    type: Directory
    inputBinding:
      position: 0
      prefix: "--fastq_dir"
  threads:
    type: int
    inputBinding:
      position: 1
      prefix: "--threads"
  hto_index_dir:
    type: Directory
    inputBinding:
      position: 2
      prefix: "--index_dir"
  hto_name:
    type: string
    inputBinding:
      position: 3
      prefix: "--name"
#   UMI_dedup_flag:
#     type: boolean
# 	inputBinding:
# 	  position: 4
# 	  prefix: "--flag"
  # out_dir:
  #   type: Directory
  #   inputBinding:
  #     position: 3
  #     prefix: "--out_dir"
outputs:
  hto_out_dir:
    type: Directory
    outputBinding:
      glob: alevin_hto