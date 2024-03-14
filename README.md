# GitHub of the thesis project: "**Multi-omics study of isocitrate dehydrogenase mutations in acute myeloid leukemia**"

* The goal of this study was to integrate mutli omique data from patients cohort Methylomes, transcriptomes, proteomes and metabolomes into multilayer regulatory network and to analyse it.

Some of the folders correspond to exploratory and preliminary analysis.

All the datasets come from public data:

- Methylomes:
  - **GSE153349**
- Transcriptomes:
  - **GSE153349**
- Proteomes:
  - **PXD023201**
- Metabolomes:
  - **GSE153349** + COMPASS algorithm

## Description of the main folders

- Datasets:
  - Chromatine_part
    - Files that associate chromatine fragment from promoter capture Hi-C data to gene names
  - **Data_4_ARACNe**
    - Script to run ARACNe bootstrap on transcriptomics to generate GRN
    - GRN from different datasets
  - Drug_screening
    - Tables from different drug screening
  - **Metabolic_datasets**
    - Table that associate reactions to gene
  - **Proteomic**
    - Proteomics data from **PXD023201**
  - **Transcriptomics**
    - A lot of transcriptomics data from public data and Molm13/K562 and IDHm Toulouse cohort
  - \+ A lot of diverse data from diverse analysis


- Scripts:
  - Chromatin_part
    - Script using chromatin assortativity aspect
  - **Epigenomique**
    - Scripts analysing DNAmethylation from raw idat data or BMIQ \+ integration of promoter capture Hi-C data \+ functional analysis
  - In_vitro
    - Some invitro results from MOLM13/K562 cells
  - **Making_Phenotype_multilayer_network_scripts**
    - **Layer_x_making.Rmd**
      - Script converting datasets to networks
    - **Combine_networks.Rmd**
      - Combine the different layers and harmonise them by adding interlayer connections
