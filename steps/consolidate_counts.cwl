cwlVersion: v1.1
class: CommandLineTool
label: Consolidate RNA, ADT and HTO counts
requirements:
  DockerRequirement:
      dockerPull: hubmap/citeseq_analysis:latest
baseCommand: /opt/salmon_wrapper.py

inputs:
  rawexpr_dir:
    type: Directory
    inputBinding:
      position: 0
	  prefix: "--out_dir"
#   adt_out_dir:
#     type: Directory
#     inputBinding:
#       position: 1
# 	  prefix: "--adt_out_dir"
#   hto_out_dir:
#     type: Directory
#     inputBinding:
#       position: 2
# 	  prefix: "--hto_out_dir"
outputs:
  combined_out_dir:
    type: Directory
    outputBinding:
      glob: combined_counts