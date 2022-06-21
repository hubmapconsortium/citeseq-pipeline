cwlVersion: v1.1
class: CommandLineTool
label: Processing RNA data for clustering and identifying marker genes
requirements:
  DockerRequirement:
      dockerPull: hubmap/citeseq_analysis:latest
baseCommand: /opt/adt_hto_salmon_wrapper.py