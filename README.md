# DRExM3L modelization of Retinitis Pigmentosa
> In this work we propose a unique approach that, starting from the set of RP disease affected genes,  provides a comprehensive landscape of the molecular mechanisms of the disease along with its druggable space. The method establishes a mechanistic disease map as an actionable environment, and employs an explainable machine learning model, “Drug REpurposing using Mechanistic Models of signal transduction and eXplainable Machine Learning” (DRExM³L), to assess the influence of druggable molecules, like drug target proteins, over the disease environment. Our approach merges information from transcriptomics, pathway graphs, biological/clinical, and drug-target interactions databases, to generate an in-depth view of the disease. The novelty of this workflow lies in the integration of multiple data sources,  reinforcing interpretability with biological knowledge while reducing the dimensionality of the datasets.

> Live demo [_here_](https://www.example.com). <!-- If you have the project hosted somewhere, include the link here. -->

## Table of Contents
* [Input Data](#input-data)
* [Technologies Used](#technologies-used)
* [Features](#features)
* [Screenshots](#screenshots)
* [Setup](#setup)
* [Usage](#usage)
* [Project Status](#project-status)
* [Acknowledgements](#acknowledgements)
* [Contact](#contact)
<!-- * [License](#license) -->

## Input Data
- The Anatomical Therapeutic Chemical (ATC) Classification table:
  File name: ATC.csv
  File path: data/raw
  Source download: https://bioportal.bioontology.org/ontologies/ATC - 2022AB csv file

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

## Technologies Used
- Tech 1 - version 1.0
- Tech 2 - version 2.0
- Tech 3 - version 3.0


## Features
List the ready features here:
- Awesome feature 1
- Awesome feature 2
- Awesome feature 3


## Screenshots
![DRExM3L modelization overview](./V4_graphical_abstract_RP_2023_rounded-Page-3.drawio.png)


## Setup
What are the project requirements/dependencies? Where are they listed? A requirements.txt or a Pipfile.lock file perhaps? Where is it located?

Proceed to describe how to install / setup one's local environment / get started with the project.


## Usage
How does one go about using it?
Provide various use cases and code examples here.

`write-your-code-here`


## Project Status
Project is: _in progress_ / _complete_ / _no longer being worked on_. If you are no longer working on it, provide reasons why.


## Acknowledgements
Give credit here.
- This project was inspired by...
- This project was based on [this tutorial](https://www.example.com).
- Many thanks to...


## Contact
Created by [@flynerdpl](https://www.flynerd.pl/) - feel free to contact me!


<!-- Optional -->
<!-- ## License -->
<!-- This project is open source and available under the [... License](). -->

<!-- You don't have to include all sections - just the one's relevant to your project -->