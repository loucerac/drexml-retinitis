# DRExM3L modelization of Retinitis Pigmentosa
> In this work we propose a unique approach that, starting from the set of RP disease affected genes, provides a comprehensive landscape of the molecular mechanisms of the disease along with its druggable space. The method establishes a mechanistic disease map as an actionable environment, and employs an explainable machine learning model, “Drug REpurposing using Mechanistic Models of signal transduction and eXplainable Machine Learning” (DRExM³L), to assess the influence of druggable molecules, like drug target proteins, over the disease environment. Our approach merges information from transcriptomics, pathway graphs, biological/clinical, and drug-target interactions databases, to generate an in-depth view of the disease. The novelty of this workflow lies in the integration of multiple data sources, reinforcing interpretability with biological knowledge while reducing the dimensionality of the datasets.

> Mechanistic demo [_here_](https://www.example.com).

## Table of Contents
* [Input Data](#input-data)
* [Screenshots](#screenshots)
* [Setup and Usage](#usage)
* [Contact](#contact)

## Input Data
- The Genotype-Tissue Expression (GTEx) RNA-Seq Data Gene read counts:
  File name: GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct
  File path: data/raw
  Source download: https://gtexportal.org/home/datasets

- The Anatomical Therapeutic Chemical (ATC) Classification table:
  File name: ATC.csv
  File path: data/raw
  Source download: https://bioportal.bioontology.org/ontologies/ATC - 2022AB CSV file

- The Human Phenotype Ontologies database (HPO):
  File name: hp.obo-v1.2-20190906
  File path: data/raw
  Source download: https://hpo.jax.org/app/data/ontology - September 6th, 2019 release

- The Human Phenotype Ontologies (HPO) annotations linking diseases and phenotypes
  File name: phenotype_annotation20190906.tab
  File path: data/raw
  Source download: https://hpo.jax.org/app/data/annotations - September 6th, 2019 release

- The Human Phenotype Ontologies (HPO) annotations linking genes and phenotypes
  File name: phenotype_to_genes20191010.txt
  File path: data/raw
  Source download: https://hpo.jax.org/app/data/annotations - September 6th, 2019 release

- HiPathia's list of physiological KEGG signaling pathways
  File name: physiological_paths.tsv
  File path: data/raw
  Source: http://hipathia.babelomics.org/
  
- HiPathia's list of physiological KEGG signaling pathways with GO/Uniprot functional annotations 
  File name: physPathsAnnot.tsv
  File path: data/raw
  Source: http://hipathia.babelomics.org/

## Screenshots
![DRExM3L model overview](./img/V4_graphical_abstract_RP_2023_rounded-Page-3.drawio.png)


## Setup and Usage
The project requires a working `conda` and a `GNU/Linux x64` system.

Install the dependencies with:

`bash install.sh`

If using a SLURM-based HPC, run the full analysis with:
`sbatch run.sbatch`

## Authors and contributors

This project is released under the MIT license.
Here are project's authors and contributors:

- Marina Esteban-Medina <marina.esteban@juntadeandalucia.es>
- Kinza Rian <kinza.rian@juntadeandalucia.es>
- Joaquin Dopazo <joaquin.dopazo@juntadeandalucia.es>
- Maria Peña-Chilet <mariapch84@gmail.com>
- Carlos Loucera <carlos.loucera@juntadeandalucia.es>
