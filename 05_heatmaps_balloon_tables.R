############################################
## Project: RP DRexML
## Script purpose:  Plot Heatmaps of relevan KDTs x stable circuits with annotations. Export KDTs tables. Plot hallmark ballon plot.
## Date: 25.02.2023
## Author: Marina Esteban-Medina
#######################################

if (!require(pacman, quietly = TRUE)){
  install.packages(pacman)
}else{
  library(pacman)
}

pacman::p_load("here", "rentrez","hipathia", "utils", "stringr", "AnnotationDbi", "org.Hs.eg.db",
               "dplyr","tidyr", "openxlsx", "data.table", "scales", "NMF", "tibble", "ggplot2", "fmsb")



if(!dir.exists(here("results/tables"))){
  dir.create(here("results/tables"))
}

if(!dir.exists(here("results/figures"))){
  dir.create(here("results/figures"))
}

if(!dir.exists(here("rds"))){
  dir.create(here("rds"))
}
### 1. Read SHAP and DRUGGANK filtered data ####

pathways <- hipathia::load_pathways("hsa")## First of all Load pathways from hipathia R package

## Load the relevance scores matrix with the threshold selection matrix for filtering ###
data_folder = here("results","ml")

shap <- fread(file = file.path(data_folder,"shap_summary_symbol.tsv"), header = T) %>% as.data.frame()
rownames(shap)<- shap$circuit_name
shap <- shap[ ,-1]

threshold <- fread(file = file.path(data_folder,"shap_selection_symbol.tsv"), header = T) %>% as.data.frame() 
rownames(threshold) <- threshold$circuit_name
threshold <- threshold[ ,-1]

## Load the filtered matrix with stable circuits and relevant KDTs
shap_relevant_stable <- read.delim(here("results/tables/shap_relevant_stable.tsv"))
shap_entrez_relevant_stable <- read.delim( here("results/tables/shap_entrez_relevant_stable.tsv"))

## Load the drug_eff DF from Drugbank DB
drugbank_effects_tar<- readRDS(here("rds", "drugbank_effects_tar.rds"))
drug_ef<- read.delim(here("results/tables/drugEffects_KDT_relevant_stable_simpl.tsv"))


### 2.  Prepare heatMap function with drug effect annotations using NMF R pack #
annot_testAheatmap <- function(matrix, annotations, title, anot_colors, colors_main = "-RdYlBu") {    
  NMF::aheatmap(matrix, color = colors_main, border_color = "white", annCol= annotations["Drug_effect"], Rowv = T, Colv = T, ## If we set Rowv and Colv to false it will cluster qwith hclust funct.
                cexCol = 4, cexRow = 10, main =  title, 
                annColors = anot_colors )
}


## Select colors for drug effects
annot.color.col <- list('Drug_effect'=c('green','lightblue', "purple", "black", "grey"))

### 3.  Heatmap with colors indicating the strength of the relevance ####

mat <- t(apply(shap, 1, function(x) x/max(abs(x)))) ## Rescale the score values so they are in a scale -1,1 on the whole matrix
mat[which(threshold == 0, arr.ind = T)] <- 0  ## set to 0 the non relevant ones
mat <- mat[rownames(shap_relevant_stable), colnames(shap_relevant_stable)]
dim(mat) ## 207 cir, 109 KDTs


## Plot the score HeatMap with the drug effect annotations on top
png(here("results","figures","heatmap_relevancesSHAP_RP2023_full_divMax_drugEffect.png"), height = 8000, width = 8000, res = 500)
annot_testAheatmap(matrix = mat, annotations = drug_ef[drug_ef$symbol %in% colnames(shap_relevant_stable),], anot_colors = annot.color.col,
                   title = "RELEVANCE SCORE OF RELEVANT KDTs x RP circuits \n Score of how relevant is each KDT for the activity of each circuit with sign")

dev.off()


### 4. Heatmap with yes/no indicating if a KDT is relevant ####

relevant_yesno <- shap_relevant_stable

## Subset  only the shap values which are relevant for at least 1 circuit in threshold
relevant_yesno[which(relevant_yesno != 0, arr.ind = T)] <- 1

## Plot the score HeatMap with the drug effect annotations on top
png(here("results", "figures","heatmap_yesnorelevancesSHAP_RP2023_drugEffect.png"), height = 8000, width = 8000, res = 500)
annot_testAheatmap(matrix = relevant_yesno, annotations = drug_ef[drug_ef$symbol %in% colnames(relevant_yesno),], anot_colors = annot.color.col,  colors_main = "-topo",
                   title = "RELEVANCE SCORE OF RELEVANT KDTs x RP circuits \n Score of how relevant is each KDT for the activity of each circuit with sign")

dev.off()


##############################################################################################
#########    TABLES PAPER 
###########################################################

### 1.  Table X1, cols : gene symbol (entrez), Drug name, drugbank ID, Action, Pharmacological action #
table_x1 <- data.frame('Gene' = paste0(drugbank_effects_tar$symbol, " (",drugbank_effects_tar$entrez_id, ")" ),
                       'Drug name' = drugbank_effects_tar$name,
                       'DrugBank ID' = drugbank_effects_tar$drugbank_id,
                       'Pharmacological.action' = drugbank_effects_tar$actions) 


length(unique(table_x1$Gene))
dim(table_x1)

## Unify the "unknown pharmacological action to "unkown"
table_x1$Pharmacological.action[grep("other|unknown", table_x1$Pharmacological.action)] <-  "unknown"
table_x1$Pharmacological.action[which(table_x1$Pharmacological.action == "") ] <-  "unknown"

write.xlsx(table_x1, file = here("results", "tables","Table_X1_genesDrugs.xlsx"))

## Genes experimentally validated
vali <- c("ALOX5", "ELOVL4", "GABRA1", "GRIN1", "SLC12A5", "GLRA2") ## 6
index_vali <- unlist(sapply(vali, function(x) grep(x, table_x1$Gene)))

table_x1_onlyValidated <- table_x1[index_vali, ]
length(unique(table_x1_onlyValidated$Gene)) ## check that we have all (6)

write.xlsx(table_x1_onlyValidated, file = here("results", "tables", "Table_X1_genesDrugs_onlyValidated.xlsx"))

### 3.  Table X2, cols : KEGG pathway, circuit eff, Function; per rel KDT: is Relevant, Value.
annotations_eff <- read.delim(here("data/raw/physPathsAnnot.tsv")) ##  GO annotations of effectors from GO db 2022

## First we read again the shap_entrez to have the values withour rescaling for the Heatmaps
shap_relevant_stable_notscaled <- shap_relevant_stable %>% add_column(circuit_code = rownames(shap_entrez_relevant_stable))

any(apply(shap_relevant_stable, 2, sum) == 0)

pivot_shapRel <- shap_relevant_stable_notscaled %>% rownames_to_column("circuit") %>% reshape2::melt(.) %>% 
  add_column(KEGG.Pathway = sapply(str_split(.$circuit, ": "), function(x) x[1]) ,.before = "circuit") %>%
  add_column(Circuit.effector = sapply(str_split(.$circuit, ": "), function(x) x[2]) ,.after = "circuit") %>% filter(., .$value != 0)

length(unique(pivot_shapRel$variable)) ## 109 OK
colnames(pivot_shapRel)<- gsub("variable", "KDT", colnames(pivot_shapRel)) %>%  gsub("value", "SHAP_score", .)

Hipathia_funct <- hipathia::get_pathways_annotations(pathway_names = unique(pivot_shapRel$circuit_code), metaginfo = pathways, dbannot = "uniprot", collapse = T)
colnames(Hipathia_funct) <- c("paths", "Hipathia_Uniprot")
pivot_shapRel <- merge(pivot_shapRel, Hipathia_funct, by.x = "circuit_code", by.y = "paths")


functions_anot <- data.frame()

for(i in pivot_shapRel$circuit_code ){
  
  functions_anot[i, "circuit"]<- i
  functions_anot[i, "GO_functions"] <- paste0(annotations_eff$Term[annotations_eff$pathway %in% i], collapse = ",")
  
}

pivot_shapRel <- merge(pivot_shapRel, functions_anot, by.x = "circuit_code", by.y = "circuit")

pivot_shapRel$Circuit.effector <- gsub("\\*", "", pivot_shapRel$Circuit.effector)


## Creat the per KDT data.frame to then cbind() them all

df_kdts <- list()

for(i in unique(pivot_shapRel$KDT) ){
  
  df <- pivot_shapRel[which(pivot_shapRel$KDT == i),]
  
  df_kdts[[i]] <- df
  
  
}

# saveRDS(df_kdts, here("rds", "Table_X2_circuits_functions_kdts.rds"))
openxlsx::write.xlsx(df_kdts , file = here("results", "tables","Table_X2_circuits_functions_kdts.xlsx"))


df_kdts_vali <- list()

for(i in unique(vali) ){
  
  df2 <- pivot_shapRel[which(pivot_shapRel$KDT == i),]
  
  df_kdts_vali[[i]] <- df2
  
  
}

openxlsx::write.xlsx(df_kdts_vali , file = here("results", "tables","Table_X2_circuits_functions_kdts_validated.xlsx"))


################################################################################################
##### PLOTS: SELECTION OF CIRCUIT FUNCTIONS FOR FUNCTIONAL MODULE ANALYSIS #######
#############################################################################################

## READ MANUALLY CURATED HALLMARKS 
hallmarks_circuits_anot <-  read.xlsx(xlsxFile =  "./data/final/RP_map_functions_MPC-annot.xlsx") %>% .[,-c(1,3,13,14,15)] %>% column_to_rownames(., "circuit_name")%>%
  apply(., 2, function(x) ifelse(is.na(x), F, T)) %>% as.data.frame(.)


## Here I am recounting to see how many hallmarks we have summarized by PATHWAY
n_cir <- data.frame(table(sapply(str_split(rownames(hallmarks_circuits_anot), ": "), function(x) x[1]))) # get the nº of circuit per pathway to obtain normalized Freq of Hallmark per pathway
n_hallmark <- data.frame(hallmark = colnames(hallmarks_circuits_anot), total = apply(hallmarks_circuits_anot, 2, sum))

table_pathways <- add_column(hallmarks_circuits_anot, "pathway" = sapply(str_split(rownames(hallmarks_circuits_anot), ": "), function(x) x[1])) %>% pivot_longer(cols = -pathway) %>% 
  mutate(value = case_when(value == "TRUE" ~ 1,
                           value == "FALSE" ~ 0))  %>% group_by(pathway, name) %>%  summarize(count = sum(value)) %>% .[-which(.$count == 0),] %>%
  add_column(n_cir = n_cir$Freq[match(.$pathway, n_cir$Var1)]) %>% ## add nº of circuit a pathway has in the RP map to normalized the plot
  add_column(total_hallmark = n_hallmark$total[match(.$name, n_hallmark$hallmark)])

table_pathways$percentage_Inpathway <- (table_pathways$count/table_pathways$n_cir)*100 ## Percentage of the hallmark normalized per nº of circuit the pathway has
table_pathways$percentage_perHallm <- (table_pathways$count/table_pathways$total_hallmark)*100 ## Percentage of the hallmark normalized by ttoal annotations fo that hallmark in the RP_map


length(unique(table_pathways$pathway)) # 40 pathways

## Make table pretty for : BALLOON  plot of Pathways and  9 hallmaarks of RP
table_pathways_pretty <- table_pathways  ## Make the names pretty for the plot
table_pathways_pretty$pathway[grep("pathway",table_pathways_pretty$pathway, invert = T)] <- paste0(table_pathways_pretty$pathway[grep("pathway",table_pathways_pretty$pathway, invert = T)], " signaling pathway")
table_pathways_pretty$name <- gsub("\\.", " ",table_pathways_pretty$name)
table_pathways_pretty$percentage_Inpathway <- round(table_pathways_pretty$percentage_Inpathway, digits = 0) ## Presence of the hallmark in the pathway: Ex 10Apoptosis_pathwayA / 5 circ_pathwayA
table_pathways_pretty$percentage_perHallm <- round(table_pathways_pretty$percentage_perHallm , digits = 0) ## Percentage of the hallmark from total annotations of that hallmark : Ex 10Apoptosis_pathwayA / 138_totalApopAnot  
colnames(table_pathways_pretty) <- c("Pathway", "Hallmark", "count", "n_cir","total_hallmark" ,"Presence", "Percentage") 
colnames(table_pathways) <- c("Pathway", "Hallmark", "count", "n_cir","total_hallmark" ,"Presence", "Percentage") 

## Balloon plot
balloon <- table_pathways_pretty %>% ggplot(aes(x = Hallmark, y = Pathway )) +
  geom_point(aes(size = n_cir, color = Percentage)) +
  scale_size_continuous(limits = c(1, max(table_pathways_pretty$n_cir)), range = c(1,8) ) +
  labs(x = NULL, y = NULL) +
  theme_light()+
  scale_color_gradient(low = "grey", high = "darkblue")+
  theme(legend.text = element_text(size = 16, family = "Helvetica"), legend.title = element_text(size = 22, family = "Helvetica"),
        axis.text.x = element_text(size = 15, family = "Helvetica", angle = 30, hjust = 1), 
        axis.text.y =  element_text(size = 15, family = "Helvetica"))+
  guides(fill = guide_legend(override.aes = list(size=8)))

png(here("results", "figures","balloonplot_pathway_hallmarks_onecolor_grad.png"), height = 8000, width = 8000, res = 500)
balloon
dev.off()

### Create the table of circuits and hallmarks in a long format
hallmarks_stable_circuits_anot <- hallmarks_circuits_anot[which(rownames(hallmarks_circuits_anot) %in% rownames(shap_relevant_stable)) ,] ## Binary col = Halmarks T/F, rows=stable_cirs 

table_hallmarks_circuit <- hallmarks_stable_circuits_anot %>% 
  rownames_to_column(., "circuit") %>%
  pivot_longer(cols = -circuit) %>% 
  mutate(value = case_when(value == "TRUE" ~ 1,
                           value == "FALSE" ~ 0))  %>% group_by(circuit, name) %>% .[-which(.$value == 0),-which(colnames(.)=="value")]
colnames(table_hallmarks_circuit) <- c("circuit", "Module")
table_hallmarks_circuit$Module <- gsub("\\.", " ",table_hallmarks_circuit$Module)

length(unique(table_hallmarks_circuit$circuit)) ## 203 /207 annotated


#### Create pivot table with SHAP + drugbank + functional annotations + modules ##########

## Export table with drugs and effects
pivot_shapRel_drugbank_functions <- pivot_shapRel %>% merge(.,table_pathways, by.x = "KEGG.Pathway", by.y = "Pathway") %>% merge(., drugbank_effects_tar[,c(10,1,2,3,4,9)], by.x = "KDT", by.y = "symbol")%>% add_column("Module"= table_hallmarks_circuit$Module[match(.$circuit, table_hallmarks_circuit$circuit)])
colnames(pivot_shapRel_drugbank_functions) <- gsub("name", "DRUG",colnames(pivot_shapRel_drugbank_functions) )

write.table(pivot_shapRel_drugbank_functions, file = here("results", "tables", "ALLpivot_cir_funct_KDT_shapScore_drug_table.tsv"), sep = "\t", quote = F, col.names = T, row.names = F)
saveRDS(pivot_shapRel_drugbank_functions, here("rds", "ALLpivot_cir_funct_KDT_shapScore_drug_table.rds"))

length(unique(pivot_shapRel_drugbank_functions$KDT))

## Export only validated 
pivot_shapRel_drugbank_functions_validated <- filter(pivot_shapRel_drugbank_functions, KDT %in% vali)
length(unique(pivot_shapRel_drugbank_functions_validated$KDT)) ## check (6)

write.table(pivot_shapRel_drugbank_functions_validated, file = here("results", "tables" ,"Validated_pivot_cir_funct_KDT_shapScore_drug_table.tsv"), sep = "\t", quote = F, col.names = T, row.names = F)
saveRDS(pivot_shapRel_drugbank_functions_validated, here("rds", "Validated_pivot_cir_funct_KDT_shapScore_drug_table.rds"))


#### GET KDTS x Hallmarks (in % of presence) for RADAR, HEATMAP and SPIDER plot ####

## KDTS ##
relevant_KDTbyfunc <- rownames_to_column(hallmarks_circuits_anot, "circuit") %>% merge(., distinct(pivot_shapRel_drugbank_functions[,c(1,4)])) %>% .[, -c(1)]  %>% pivot_longer(cols = -KDT) %>%
  mutate(value = case_when(value == "TRUE" ~ 1,
                           value == "FALSE" ~ 0))  %>% group_by(KDT, name) %>%
  summarize(count = sum(value)) %>% .[-which(.$count == 0),] %>%
  pivot_wider(names_from = name , values_from = count ) %>% column_to_rownames("KDT")## Here we transform the DF into a matrix of counts of hallmark per KDT
relevant_KDTbyfunc[is.na(relevant_KDTbyfunc)] <- 0
relevant_KDTbyfunc <- rbind( TOTAL = n_hallmark$total[match(colnames(relevant_KDTbyfunc),n_hallmark$hallmark)] ,relevant_KDTbyfunc) %>% 
  .[, c("Fatty.acid.and.lipid.metabolism", "Apoptosis", "Neuronal", "DNA.integrity" ,"Inflammatory.response", "Stress.response","Necrosis", "Development","Sensory.and.stimuli.transduction")] ## reorder according to hallmark circle

saveRDS(relevant_KDTbyfunc, file = here("rds", "relevant_KDTbyfunc.rds"))


## DRUGS ##
drug_byfunc <- rownames_to_column(hallmarks_circuits_anot, "circuit") %>% merge(., distinct(pivot_shapRel_drugbank_functions[,c(4,18)])) %>% .[, -c(1)]  %>% pivot_longer(cols = -DRUG) %>%
  mutate(value = case_when(value == "TRUE" ~ 1,
                           value == "FALSE" ~ 0))  %>% group_by(DRUG, name) %>%
  summarize(count = sum(value)) %>% .[-which(.$count == 0),] %>%
  pivot_wider(names_from = name , values_from = count ) %>% column_to_rownames("DRUG") ## Here we transform the DF into a matrix of counts of hallmark per KDT
drug_byfunc[is.na(drug_byfunc)]<- 0
drug_byfunc <- rbind( TOTAL = n_hallmark$total[match(colnames(drug_byfunc),n_hallmark$hallmark)] ,drug_byfunc) %>%
  .[, c("Fatty.acid.and.lipid.metabolism", "Apoptosis", "Neuronal", "DNA.integrity" ,"Inflammatory.response", "Stress.response","Necrosis", "Development","Sensory.and.stimuli.transduction")] ## reorder according to hallmark circle

saveRDS(drug_byfunc, file = here("rds", "drug_byfunc.rds"))
