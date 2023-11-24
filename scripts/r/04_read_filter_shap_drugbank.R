#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
  db_name <- "genes_drugbank-v050108_mygene-20230120.tsv"
} else {
  db_name <- args[1] 
}

############################################
## Project: RP DRexML
## Script purpose:  Read and filter SHAP model results and drugbank db to plot and interpret.
## Date: 25.02.2023
## Author: Marina Esteban-Medina
#######################################

library("here")
library("rentrez")
library("hipathia")
library("utils")
library("stringr")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("dplyr")
library("tidyr")
library("openxlsx")
library("data.table")
library("scales")
library("NMF")
library("tibble")
library("ggplot2")
library("fmsb")
library("venn")
library("ggpolypath")

#### 0. Path setup ########

tables_folder <- here("results", "tables")
if(!dir.exists(tables_folder)){
  dir.create(tables_folder)
}

figures_folder <- here("results", "figures")
if(!dir.exists(figures_folder)){
  dir.create(figures_folder)
}

rds_folder <- here("results", "rds")
if(!dir.exists(rds_folder)){
  dir.create(rds_folder)
}

#### 1. Load SHAP model results and filter ########

data_folder = here("results", "ml")

## First of all Load pathways from hipathia R package
pathways <- hipathia::load_pathways("hsa")

## Load the relevance scores matrix with the threshold selection matrix for filtering ###
shap <- fread(file = file.path(data_folder,"shap_summary_symbol.tsv"), header = T) %>% as.data.frame()
rownames(shap)<- shap$circuit_name
shap <- shap[ ,-1]

shap_entrez <-  fread(file = file.path(data_folder,"shap_summary.tsv"), header = T) %>% as.data.frame()
rownames(shap_entrez)<- shap_entrez$V1
shap_entrez <- shap_entrez[ ,-1] 

threshold <- fread(file = file.path(data_folder,"shap_selection_symbol.tsv"), header = T) %>% as.data.frame() 
rownames(threshold) <- threshold$circuit_name
threshold <- threshold[ ,-1]

threshold_entrez <- fread(file = file.path(data_folder,"shap_selection.tsv"), header = T) %>% as.data.frame() 
rownames(threshold_entrez) <- threshold_entrez$V1
threshold_entrez <- threshold_entrez[ ,-1]

## Create tables of genes and circuits translations
genes_tr <- data.frame(entrez = colnames(shap_entrez), symbol = colnames(shap))
write.table(genes_tr, file.path(tables_folder, "KDT_genes_translate.tsv"), quote = F, sep = "\t", col.names = T, row.names = F)

circuits_tr <- data.frame(code = rownames(shap_entrez), name =  rownames(shap)) %>% 
  add_column(circuit = strsplit(.$code, "\\.") %>% sapply(., function(x){ ifelse (length(x) == 3, paste(paste(x[[1]],x[[2]], sep = "-"), x[[3]], sep =  "-"),
                                                                                  ifelse(length(x) == 4, paste(paste(paste(x[[1]],x[[2]], sep = "-"), x[[3]], sep =  "-"), x[[4]], sep = " "),
                                                                                         ifelse(length(x) == 5, paste(paste(paste(paste(x[[1]],x[[2]], sep = "-"), x[[3]], sep =  "-"), x[[4]], sep = " "), x[[5]], sep = " "),
                                                                                                ifelse(length(x) == 6, paste(paste(paste(paste(paste(x[[1]],x[[2]], sep = "-"), x[[3]], sep =  "-"), x[[4]], sep = " "), x[[5]], sep = " "), x[[6]], sep = " ")))))}))
write.table(circuits_tr,file.path(tables_folder,"circuits_translate.tsv"), quote = F, sep = "\t", col.names = T, row.names = F)

## Load the stability results for all circuits from the RP Map ### 
stability <- fread(file = file.path(data_folder, "stability_results.tsv")) %>% as.data.frame(.)
rownames(stability) <-stability$name
stability <- stability[-1,-1]
stability$cir_name <- circuits_tr$name[match(rownames(stability), circuits_tr$code)]
stability$cir_name <- rownames(shap) ## Add a column with the circuit names to stability so its easier to filter based on stability
stability[is.na(stability)] = 0

## Filter shap values from circuit with stability >= 0.4
shap_stable <- shap[rownames(shap) %in% stability$cir_name[stability$stability >= 0.4], ]
threshold_stable <- threshold[rownames(threshold) %in% stability$cir_name[stability$stability >= 0.4],]

shap_entrez_stable <- shap_entrez[rownames(shap_entrez) %in% rownames(stability)[stability$stability >= 0.4],]
threshold_entrez_stable <-  threshold_entrez[rownames(threshold_entrez) %in% rownames(stability)[stability$stability >= 0.4],]


#### 2.Filter only relevant targets  from the SHAP model results ####
# Subset  only the shap values which are relevant for at least 1 circuit in threshold
shap_stable[which(threshold_stable == 0, arr.ind = T)] <- 0 ## first we set to 0 the values of the non selected shap for those circuits
shap_relevant_stable <- shap_stable[,(apply(threshold_stable, 2, function(y) any(y == 1)))] ## Filter out the relevant KDTs
dim(shap_relevant_stable)

write.table(shap_relevant_stable,file.path(tables_folder,"shap_relevant_stable.tsv"), quote = F, sep = "\t", col.names = T, row.names = T)

## Subset  only the shap values which are relevant for at least 1 circuit in threshold matrix ###
shap_entrez_stable[which(threshold_entrez_stable == 0, arr.ind = T)] <- 0 
shap_entrez_relevant_stable <- shap_entrez_stable[,(apply(threshold_entrez_stable, 2, function(y) any(y == 1)))]
dim(shap_entrez_relevant_stable)

rownames(shap_entrez_relevant_stable)<- strsplit(rownames(shap_entrez_relevant_stable), "\\.") %>% sapply(., function(x){ ifelse (length(x) == 3, paste(paste(x[[1]],x[[2]], sep = "-"), x[[3]], sep =  "-"),
                                                                                                                                  ifelse(length(x) == 4, paste(paste(paste(x[[1]],x[[2]], sep = "-"), x[[3]], sep =  "-"), x[[4]], sep = " "),
                                                                                                                                         ifelse(length(x) == 5, paste(paste(paste(paste(x[[1]],x[[2]], sep = "-"), x[[3]], sep =  "-"), x[[4]], sep = " "), x[[5]], sep = " "),
                                                                                                                                                ifelse(length(x) == 6, paste(paste(paste(paste(paste(x[[1]],x[[2]], sep = "-"), x[[3]], sep =  "-"), x[[4]], sep = " "), x[[5]], sep = " "), x[[6]], sep = " ")))))})

targets_shap_rel_stable <- data.frame(  Gene_symbol = colnames(shap_relevant_stable),
                                        entrez = genes_tr$entrez[match(colnames(shap_relevant_stable), genes_tr$symbol)],
                                        Gene_name =  AnnotationDbi::mapIds(org.Hs.eg.db, keys, keys = colnames(shap_relevant_stable),keytype  = "SYMBOL", column = "GENENAME"),
                                        Entrez_dbi  =  AnnotationDbi::mapIds(org.Hs.eg.db, keys, keys = colnames(shap_relevant_stable), keytype = "SYMBOL", column = "ENTREZID"),
                                        stringsAsFactors = F)

all(targets_shap_rel_stable$Entrez_dbi == targets_shap_rel_stable$entrez) ## TRUE
any(is.na(targets_shap_rel_stable$Gene_name))

## Table of all relevant stable KDTs with their names ##
write.table(shap_entrez_relevant_stable, file.path(tables_folder,"shap_entrez_relevant_stable.tsv"), quote = F, sep = "\t", col.names = T, row.names = T)
write.xlsx(targets_shap_rel_stable[,c(1:3)], file = file.path(tables_folder, "relevant_targets_ml_stable.xlsx"))

#### 3. Filter DRUGBANK data for functional annotation of relevant targets #

## Load Drugbank database 
data_folder2 <- here("data/interim")
fname2 <- "drugbank-v050108.tsv"
fpath2 <-file.path(data_folder2,fname2)

drugbank_alltar <- read.delim(file = fpath2, sep = "\t" )
drugbank_curated <- drugbank_alltar
drugbank_curated$is_protein_group_target <- ifelse(drugbank_curated$is_protein_group_target == "False", FALSE, "ERROR")

## Load ammendments file with curated info to specify column "actions" for drug's pharmacological actions on targets
amendments <- read.delim(here("data", "raw", "amendments_drugActions_drugbank-v050108.tsv"))

## We will use data.table to update the values of column "actions" from the ammendments file info
setDT(amendments)
setDT(drugbank_curated)
drugbank_curated[amendments, on = c("drugbank_id", "uniprot_id", "known_action", "organism"), actions := i.actions] 

write.table(drugbank_curated, here("data","interim","drugbank-v050108_curated.tsv"))

##  Total KDT_drug_combinations
distinct(drugbank_alltar[, c(1,15)]) %>% dim(.) ## 26979
length(unique(drugbank_alltar$drugbank_id)) ## 7919 drugs
length(unique(drugbank_alltar$uniprot_id)) ## 5004 drugs

drugbank_app_action <- drugbank_alltar[-(base::grep("withdrawn",drugbank_alltar$groups)),] %>% .[(base::grep("^approved|approved,", .$groups)),] %>% .[(base::grep("target", .$category)),]%>%
  .[(base::grep("yes", .$known_action)),] %>% .[.$organism == "Humans",]### version 5.1.8 from parser 2022

dim(drugbank_app_action) ## 2690 combinations KDT-drug taken into account
length(unique(drugbank_app_action$uniprot_id)) ## 718 Uniprot IDs filtered 

## Load genes stable translator used to Unipro IDs
entrez_uniprot <- read.delim(file = file.path(data_folder2, db_name)) %>% .[.$uniprot_id %in% drugbank_app_action$uniprot_id , ]
dim(entrez_uniprot)

drugbank_app_action <- merge(drugbank_app_action, entrez_uniprot,) %>% .[-which(is.na(.$entrez_id)),] 
length(unique(drugbank_app_action$entrez_id)) 
all(targets_shap_rel_stable$entrez %in% drugbank_app_action$entrez_id) ## check 

saveRDS(drugbank_app_action, file.path(rds_folder, "drugbank_app_action_entrez.rds"))

drugbank_app_action$entrez_id <- as.character(drugbank_app_action$entrez_id)
drugbank_app_action_genes <- merge(drugbank_app_action, genes_tr, by.x = "entrez_id", by.y = "entrez" ) ## Add the symbol column from the stable translate table
drugbank_app_action_genes$actions[drugbank_app_action_genes$actions == ""] <- "Other" ## unify the not known drug actions to the same term "unknown"
drugbank_app_action_genes$actions[grep(drugbank_app_action$actions, pattern = "unknown|other/unknown")] <- "Other" ## unify the not known drug actions to the same term "unknown"
any(is.na(drugbank_app_action_genes$actions))
data.frame(table(drugbank_app_action_genes$actions))%>% write.table(., file = here("data","interim", "drug_actions.tsv"), sep = "\t", quote = F, col.names = T, row.names = F) ## Create drug effect table

drug_effects_translation <- fread(here("data","raw", "drug_actions_withSimplAction.csv")) ## Read edited table with simplified drug effect
drugbank_app_action_genes$simplified_action <- drug_effects_translation$Drug_eff_simpl[match(drugbank_app_action_genes$actions, drug_effects_translation$Var1)]

length(unique(drugbank_app_action_genes$entrez_id)) ## We have all included KDTs 711
write.xlsx(drugbank_app_action_genes, file = file.path(tables_folder, "supp_tabl3_drugbank518_filtered.xlsx"))

alldrug_byaction <- drugbank_app_action_genes[,c("name", "actions", "entrez_id")]
colnames(alldrug_byaction) <- c("drug", "drug_action", "KDT")  

## Add simplified drug effects
alldrug_byaction$Drug_effect <- drug_effects_translation$Drug_eff_simpl[match(alldrug_byaction$drug_action, drug_effects_translation$Var1)]
alldrug_byaction$Drug_effect <- factor(alldrug_byaction$Drug_effect, levels = unique(alldrug_byaction$Drug_effect))
length(unique(alldrug_byaction$drug)) ## 1410 unique drugs
alldrug_byaction$drugKDT <- paste0(alldrug_byaction$drug, alldrug_byaction$KDT)


df_venn <- alldrug_byaction[, c("drug","Drug_effect")]  %>% add_column(value = as.numeric(1)) %>% dplyr::group_by(drug, Drug_effect) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n >= 1L) %>% pivot_wider(names_from = Drug_effect, values_from = n )  %>% column_to_rownames("drug")

df_venn[is.na(df_venn)] <- 0
df_venn[df_venn>1] <- 1

## Do VennDiagram of the simplified drug effects with venn R pack - ALL DRUGS ASSESSED
png(filename = file.path(figures_folder,"venn_drugeffectsALL.png"), width = 8000, height = 8000, bg = F)
venn(df_venn, ilab=TRUE, zcolor = "style")
dev.off()

#### 4.  Filter and organize the drugs and their gene targets #
drugbank_effects_tar <- drugbank_app_action[which(drugbank_app_action$entrez_id %in% targets_shap_rel_stable$entrez), 
                                            which(colnames(drugbank_app_action) %in% c("uniprot_id","drugbank_id","name", "type" , "groups", "categories", "description", "actions", "entrez_id" )) ] 
drugbank_effects_tar <- merge(drugbank_effects_tar, genes_tr, by.x = "entrez_id", by.y = "entrez" ) ## Add the symbol column from the stable translate table

drugbank_effects_tar$actions[drugbank_effects_tar$actions == ""] <- "unknown" ## unify the not known drug actions to the same term "unknown"
drugbank_effects_tar$actions[drugbank_effects_tar$actions == "other"]  <- "unknown" ## unify the not known drug actions to the same term "unknown"

data.frame(table(drugbank_effects_tar$actions)) 
any(is.na(drugbank_effects_tar$actions))

saveRDS(drugbank_effects_tar, file.path(rds_folder, "drugbank_effects_tar.rds"))

## Create a data frame to do Venn diagram of the simplified drug effects - RELEVANT DRUGS FILTERED BY DREXM3L
drugs_rel_DF <- drugbank_effects_tar[, c("name", "actions", "symbol")]

drugs_rel_DF$Drug_effect <-drug_effects_translation$Drug_eff_simpl[match(drugs_rel_DF$actions, drug_effects_translation$Var1)] ## Add the simplified drug effects

df_venn_relevant <- drugs_rel_DF[, c("name","Drug_effect")]  %>% add_column(value = as.numeric(1)) %>% dplyr::group_by(name, Drug_effect) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n >= 1L) %>% pivot_wider(names_from = Drug_effect, values_from = n )  %>% column_to_rownames("name")

df_venn_relevant[is.na(df_venn_relevant)] <- 0
df_venn_relevant[df_venn_relevant>1] <- 1

## Do VennDiagram of the simplified drug effects with venn R pack - ALL DRUGS ASSESSED
png(filename = file.path(figures_folder,"venn_drugeffects_relevantDRUGS.png"), width = 8000, height = 8000, bg = F)
venn(df_venn_relevant, ilab=TRUE, zcolor = "style")
dev.off()


## Since a drug can have several effects depending on the KDT it targets, and KDTS can be targeted by many drugs, we will select the most common effect that drugs have on each KDT to then plot it .
drug_bygenes <- data.frame(aggregate(cbind(as.character(drugbank_effects_tar$actions)) ~ drugbank_effects_tar$symbol, data = drugbank_effects_tar , FUN = paste))
colnames(drug_bygenes) <- c("gene", "drug_action") 
drug_bygenes$drug_action <- sapply(drug_bygenes$drug_action, function(x) names(which.max(table(x)))) 

drug_ef <- data.frame( symbol = targets_shap_rel_stable$Gene_symbol, 
                       Drug_effect = drug_bygenes$drug_action [match( targets_shap_rel_stable$Gene_symbol, drug_bygenes$gene)]) 
any(is.na(drug_ef))
table(drug_ef$Drug_effect)

## Simplify the drug action vector to the selected categories: "Activator", "Ligand", "Other", "Inhibitor"
drug_ef$Drug_effect <-drug_effects_translation$Drug_eff_simpl[match(drug_ef$Drug_effect, drug_effects_translation$Var1)] ## Add the simplified drug effects

dim(drug_ef)
table(drug_ef$Drug_effect)

## See which function is the most predominant for a certain gene KDT
colnames(shap_relevant_stable)[!colnames(shap_relevant_stable) %in% drug_ef$symbol] ## Check that all KDTs have a drug effect
write.table(drug_ef, file.path(tables_folder, "drugEffects_KDT_relevant_stable_simpl.tsv"), quote = F, sep = "\t", col.names = T, row.names = F)
