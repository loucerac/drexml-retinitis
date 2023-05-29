############################################
## Project: RP DRexML
## Script purpose:  Plot Heatmaps of relevan KDTs x stable circuits with annotations. Export KDTs tables. Plot hallmark balloon plot.
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

### 1. Read SHAP and DRUGGANK filtered data ####

data_folder = here("results", "ml")

pathways <- hipathia::load_pathways("hsa")## First of all Load pathways from hipathia R package

## Load the relevance scores matrix with the threshold selection matrix for filtering ###
shap <- fread(file = file.path(data_folder,"shap_summary_symbol.tsv"), header = T) %>% as.data.frame()
rownames(shap)<- shap$circuit_name
shap <- shap[ ,-1]

threshold <- fread(file = file.path(data_folder,"shap_selection_symbol.tsv"), header = T) %>% as.data.frame() 
rownames(threshold) <- threshold$circuit_name
threshold <- threshold[ ,-1]

## Load the filtered matrix with stable circuits and relevant KDTs
shap_relevant_stable <- read.delim(file.path(tables_folder, "shap_relevant_stable.tsv"))
shap_entrez_relevant_stable <- read.delim(file.path(tables_folder, "shap_entrez_relevant_stable.tsv"))

## Load the drug_eff DF from Drugbank DB
drugbank_effects_tar<- readRDS(file.path(rds_folder, "drugbank_effects_tar.rds"))
drug_ef<- read.delim(file.path(tables_folder, "drugEffects_KDT_relevant_stable_simpl.tsv"))


### 2.  Prepare heatMap function with drug effect annotations using NMF R pack #
annot_testAheatmap <- function(matrix, annotations, title, anot_colors, colors_main = "-RdYlBu") {    
  NMF::aheatmap(matrix, color = colors_main, border_color = "white", annCol= annotations["Drug_effect"], Rowv = T, Colv = T, ## If we set Rowv and Colv to false it will cluster qwith hclust funct.
                cexCol = 4, cexRow = 10, main =  title, 
                annColors = anot_colors )
}

shap_scaled <- t(apply(shap, 1, function(x) x/max(abs(x)))) 
shap_scaled<- shap_scaled[rownames(shap_relevant_stable), colnames(shap_relevant_stable)]
  
png(file.path(figures_folder, "heatmap_SHAP_RP2023_full_abstract.png"), height = 8000, width = 8000, res = 500)
NMF::aheatmap(shap_scaled, color = "-RdYlBu", border_color = "white", Rowv = T, Colv = T, ## If we set Rowv and Colv to false it will cluster qwith hclust funct.
              cexCol = 4, cexRow = 10)
dev.off()


## Select colors for drug effects
annot.color.col <- list('Drug_effect'=c('green','lightblue', "purple", "black", "grey"))

### 3.  Heatmap with colors indicating the strength of the relevance ####

mat <- t(apply(shap, 1, function(x) x/max(abs(x)))) ## Rescale the score values so they are in a scale -1,1 on the whole matrix
mat[which(threshold == 0, arr.ind = T)] <- 0  ## set to 0 the non relevant ones
mat <- mat[rownames(shap_relevant_stable), colnames(shap_relevant_stable)]
dim(mat) ## 207 cir, 109 KDTs


## Plot the score HeatMap with the drug effect annotations on top
png(file.path(figures_folder, "heatmap_relevancesSHAP_RP2023_full_divMax_drugEffect.png"), height = 8000, width = 8000, res = 500)
annot_testAheatmap(matrix = mat, annotations = drug_ef[drug_ef$symbol %in% colnames(shap_relevant_stable),], anot_colors = annot.color.col,
                   title = "RELEVANCE SCORE OF RELEVANT KDTs x RP circuits \n Score of how relevant is each KDT for the activity of each circuit with sign")

dev.off()


### 4. Heatmap with yes/no indicating if a KDT is relevant ####

relevant_yesno <- shap_relevant_stable

## Subset  only the shap values which are relevant for at least 1 circuit in threshold
relevant_yesno[which(relevant_yesno != 0, arr.ind = T)] <- 1

## Plot the score HeatMap with the drug effect annotations on top
png(file.path(figures_folder, "heatmap_yesnorelevancesSHAP_RP2023_drugEffect.png"), height = 8000, width = 8000, res = 500)
annot_testAheatmap(matrix = relevant_yesno, annotations = drug_ef[drug_ef$symbol %in% colnames(relevant_yesno),], anot_colors = annot.color.col,  colors_main = "-topo",
                   title = "RELEVANCE SCORE OF RELEVANT KDTs x RP circuits \n Score of how relevant is each KDT for the activity of each circuit with sign")

dev.off()


########## STUDY OF KDT SCORES ##########

## Boxplots of top scored KDTs
scores_KDTS <- data.frame(score = apply(shap_relevant_stable, 2, function(x) mean(abs(x)))) %>% rownames_to_column("KDT") %>% .[order(.$score, decreasing = T),]

top30_boxplot_df <- data.frame(t(mat[, scores_KDTS$KDT[1:30]])) %>%  rownames_to_column("KDT") %>% pivot_longer(-KDT) %>% filter(value!=0) %>% 
  add_column(SIGN = ifelse(.$value > 0, "POSITVE", "NEGATIVE"))
top30_boxplot_df$value <- abs(top30_boxplot_df$value)
top30_boxplot_df$KDT <- reorder(top30_boxplot_df$KDT, -top30_boxplot_df$value, FUN = median) ## order boxplots to obtain decreasing medians
top30_boxplot_df$SIGN <- factor(top30_boxplot_df$SIGN, levels = c("POSITVE", "NEGATIVE"))

png(file.path(figures_folder, "boxplots_top30_kdt_score.png"), height = 8000, width = 12000, res = 500)
  ggplot(data =  top30_boxplot_df, aes(x = KDT, y = value,  fill= SIGN))+
    geom_boxplot(stat = "boxplot")+  
    theme_minimal()+
    theme(legend.title =  element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 40, family = "special", hjust = 1, colour = "black"),
          axis.text.x = element_text(size = 35, family = "helvetica",colour = "black", angle = 90, hjust = 1, vjust = 1),
          axis.title.x = element_blank(),
          legend.text = element_text(size = 35),
          legend.key.size = unit(2,"cm"))
dev.off()

## Assess the distribution of the KDT signs and scores
kdt_scores <- as.data.frame(t(mat)) %>% rownames_to_column(.,"KDTS") %>% pivot_longer(-KDTS) %>% 
  add_column("SIGN"= ifelse(.$value > 0, "POSITIVE", "NEGATIVE")) %>% add_column("path"= sapply(str_split(.$name, ": "), function(x) x[[1]])) %>% filter(.,value != 0)
kdt_scores$path[grep("pathway", kdt_scores$path, invert = T)] <- paste(kdt_scores$path[grep("pathway", kdt_scores$path, invert = T)], "pathway", sep = " ") ## Add the word pathway to the path names
# kdt_scores$value <-  abs(kdt_scores$value)
kdt_scores <- kdt_scores[order(abs(kdt_scores$value), decreasing = F), ]
kdt_scores$KDTS <- reorder(kdt_scores$KDTS, -abs(kdt_scores$value), FUN = sum, decreasing = T) ### Order the KDTs
  
png(file.path(figures_folder, "kdts_scores_stacked_vertical_onlyrelevant.png"), height = 12000, width = 8000, res = 500)
ggplot(data=kdt_scores, aes(x= KDTS, y= value, fill= SIGN))+
  geom_bar(stat="identity")+
  # geom_text(aes(label = percentage), hjust = -0.1, size = 2, family = "candara", colour = "#040f42" )+
  theme_classic()+
  theme(legend.title =  element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 15, family = "helvetica", hjust = 1, colour = "black"),
        axis.text.x = element_text(size = 20, family = "helvetica",colour = "black"),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 16))+
  coord_flip()+
   scale_y_continuous(limits = c(-100,60), breaks = seq(-100, 60, 10))+ ## When stacked
  # scale_y_continuous(limits = c(-100,100))+ ## when non stacked
  scale_fill_manual(values = c("#4575B4","#F46D43"))
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


write.xlsx(table_x1, file = file.path(tables_folder,"Table_X1_genesDrugs.xlsx"))

## Genes experimentally validated
vali <- c("ALOX5", "ELOVL4", "GABRA1", "GRIN1", "SLC12A5", "GLRA2") ## 6
index_vali <- unlist(sapply(vali, function(x) grep(x, table_x1$Gene)))

table_x1_onlyValidated <- table_x1[index_vali, ]
length(unique(table_x1_onlyValidated$Gene)) ## check that we have all (6)

write.xlsx(table_x1_onlyValidated, file = file.path(tables_folder, "Table_X1_genesDrugs_onlyValidated.xlsx"))

### 3.  Table X2, cols : KEGG pathway, circuit eff, Function; per rel KDT: is Relevant, Value.
annotations_eff <- read.delim(here("data","raw", "physPathsAnnot.tsv")) ##  GO annotations of effectors from GO db 2022

## First we read again the shap_entrez to have the values withour rescaling for the Heatmaps
shap_relevant_stable_notscaled <- shap_relevant_stable %>% add_column(circuit_code = rownames(shap_entrez_relevant_stable))
saveRDS(shap_relevant_stable_notscaled, file =  file.path(rds_folder,"shap_relevant_stable_notscaled_hallmarks.rds"))

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

saveRDS(df_kdts, file.path(rds_folder,"Table_X2_circuits_functions_kdts.rds"))
openxlsx::write.xlsx(df_kdts , file = file.path(tables_folder,"Table_X2_circuits_functions_kdts.xlsx"))


df_kdts_vali <- list()

for(i in unique(vali) ){
  
  df2 <- pivot_shapRel[which(pivot_shapRel$KDT == i),]
  
  df_kdts_vali[[i]] <- df2
  
  
}

openxlsx::write.xlsx(df_kdts_vali , file.path(tables_folder,"Table_X2_circuits_functions_kdts_validated.xlsx"))
