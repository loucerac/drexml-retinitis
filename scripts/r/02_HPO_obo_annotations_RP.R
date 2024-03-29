#!/usr/bin/env RScript

#####################################################################################
######## ORGANIZE AND ANNOTATE HPO TERMS RELATED WITH GENES AND DISEASES ############
#####################################################################################

library(here)
library(data.table)
library(stringr)
library(dplyr)
library(ontologyIndex)
library(janitor)


## Load HPO OBO ontology file
## Obtained from https://hpo.jax.org/app/data/ontology 2019
file_HPO.obo <- "hp.obo-v1.2-20190906"
hpo <- get_ontology(file = here("data", "raw", file_HPO.obo), propagate_relationships = c("is_a", "part_of"))

##  Load phenotype_annotation.tab
## Obtained from https://hpo.jax.org/app/data/annotations version 2019
hpo_tab <- data.frame(fread(here("data", "raw", "phenotype_annotation20190906.tab"), fill = T, header = F, sep = "\t"), stringsAsFactors = T)

hpo_tab <- hpo_tab[, which(colSums(hpo_tab != "") != 0)] %>% remove_constant(.)
# remove repeated disease in capital letters
hpo_tab <- hpo_tab[, -11]

## v.2019
colnames(hpo_tab) <- c("db", "code", "disease", "hpo_qualifier", "hpo_id", "disease_resource", "evidence", "hpo_cc_onset", "hpo_freq", "aspect", "origin", "freq")


terms <- lapply(hpo_tab$hpo_id, function(x) get_term_property(hpo, "name", x)) %>%
  unlist(.) %>%
  stack(.) %>%
  mutate_all(as.character)
colnames(terms) <- c("term", "hpo_code")

hpo_tab$term_id <- terms$term


## Load phenotype_to_genes.txt file from HPOdb (downloaded in local)
## Obtained from https://hpo.jax.org/app/data/annotations version 2019
hpo2genes <- data.frame(read.delim(skip = 1, file = (here("data", "raw", "phenotype_to_genes20191010.txt")), fill = T, header = F, sep = "\t"), stringsAsFactors = T)
colnames(hpo2genes) <- c("hpo_id", "term_id", "entrez", "symbol")

hpo_genes_long <- data.frame(aggregate(cbind(hpo2genes$entrez) ~ as.character(hpo2genes$hpo_id), data = hpo2genes, FUN = paste), stringsAsFactors = F)
colnames(hpo_genes_long) <- c("HPO", "entrez_id")

genes_hpo_long <- data.frame(aggregate(cbind(as.character(hpo2genes$hpo_id)) ~ hpo2genes$entrez, data = hpo2genes, FUN = paste), stringsAsFactors = F)
colnames(genes_hpo_long) <- c("entrez_id", "HPO")


## TAG HPO phenotype_annotation.tab table with levels

# Tag all table

p <- hpo_tab$hpo_id

ancestors_p <- sapply(p, function(x) {
  get_term_property(ontology = hpo, property = "ancestors", term = x)
})

n_ancestors_p <- lapply(ancestors_p, function(x) {
  length(x)
})

p_levelstag <- stack(n_ancestors_p)
colnames(p_levelstag) <- c("level", "hpo_code")

p_levelstag$hpo_code <- as.character(p_levelstag$hpo_code)

hpo_tab$level <- p_levelstag$level # Add the level column to the hpo_tab table

hpo_tab <- hpo_tab[-base::grep("*obsolete*", hpo_tab$term_id), ]

## We only want the orphadata and P (P_abnormality) from the $aspect column of the table

## Filter only HPO associated with ORPHA diseases
hpo_or <- hpo_tab[hpo_tab$origin == "orphadata", ]
# erase duplicated info and empty columns
#TODO: always use column names (except in linear algebra)
# hpo_or_p <- hpo_or[hpo_or$aspect == "P", -c(1, 2, 6, 11)]
hpo_or_p <- hpo_or[hpo_or$aspect == "P", -c(1, 2, 11)]
#ERROR: 6 points to `disease_resource` which is later used in `modules/r/03_genes_circuits_ML.R`
# programatic way of removing duplicated columns (keeps the first)
#hpo_or <- hpo_or[!duplicated(as.list(hpo_or))]
# table(hpo_or_p$aspect)


## Add gene - HPO relations to the hpo_or_p table of ORPHA phenotypic_ab HPO codes.
## Tabla con ORPHAdis_HPO_pheno_ab_entrez_levels
hpo_or_p$entrezs <- hpo_genes_long$entrez_id[match(hpo_or_p$hpo_id, hpo_genes_long$HPO)]

saveRDS(hpo_or_p, file = here("data", "interim", "hpo_or_p_entrezs05122019.rds"))
