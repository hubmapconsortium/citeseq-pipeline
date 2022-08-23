cwlVersion: v1.1
class: CommandLineTool
label: Downstream analysis for RNA and ADT 
requirements:
  DockerRequirement:
      dockerPull: hubmap/citeseq_analysis:latest
baseCommand: /opt/downstream.py

inputs:
  muon_dir:
    type: File
    inputBinding:
      position: 0
      prefix: "--muon_dir"
outputs:
  muon_processed:
    type: File
    outputBinding:
      glob: "citeseq_normalized.h5mu"
  mofa_out:
    type: File
    outputBinding:
      glob: "citeseq_mofa.hdf5"
  rna_embedding:
    type: File
    outputBinding:
      glob: "leiden_cluster_rna.pdf"
  adt_embedding:
    type: File
    outputBinding:
      glob: "leiden_cluster_adt.pdf"
  joint_embedding:
    type: File
    outputBinding:
      glob: "leiden_cluster_combined.pdf"
  

  
  