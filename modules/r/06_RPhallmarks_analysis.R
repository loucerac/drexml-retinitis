############################################
## Project: RP DRexML
## Script purpose:  Plot hallmark balloon plot. Create tables of hits per hallmark (drugs, relevan KDTs, hallmarks)
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

################################################################################################
##### PLOTS: SELECTION OF CIRCUIT FUNCTIONS FOR FUNCTIONAL MODULE ANALYSIS #######
#############################################################################################

data_folder = here("results", "ml")

### Load files###
shap_relevant_stable_notscaled <- readRDS(file =  file.path(rds_folder,"shap_relevant_stable_notscaled_hallmarks.rds")) ## Load shap matrix
shap_relevant_stable <- read.delim(file = file.path(tables_folder,"shap_relevant_stable.tsv")) ## Load the filtered matrix with stable circuits and relevant KDTs
drugbank_effects_tar<- readRDS(file.path(rds_folder, "drugbank_effects_tar.rds")) ## Load the drug_eff DF from Drugbank DB

vali <- c("ALOX5", "ELOVL4", "GABRA1", "GRIN1", "SLC12A5", "GLRA2") ## 6 validated KDTs

hallmarks_circuits_anot <-  read.xlsx(xlsxFile =  "./data/final/RP_map_functions_MPC-annot.xlsx") %>% 
  .[,-c(1,3,13,14,15)] %>%
  column_to_rownames(., "circuit_name") %>%
  apply(., 2, function(x) ifelse(is.na(x), F, T)) %>% 
  as.data.frame(.)

## Here I am recounting to see how many hallmarks we have summarized by PATHWAY ###
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
colnames(table_pathways_pretty) <- c("Pathway", "Hallmark", "count", "n circuits","total_hallmark" ,"Presence", "Percentage") 
colnames(table_pathways) <- c("Pathway", "Hallmark", "count", "n_cir","total_hallmark" ,"Presence", "Percentage") 

## Balloon plot
balloon <- table_pathways_pretty %>% ggplot(aes(x = Hallmark, y = Pathway )) +
  geom_point(aes(size = `n circuits`, color = Percentage)) +
  scale_size_continuous(limits = c(1, max(table_pathways_pretty$`n circuits`)), range = c(1,12), breaks = c(5,10,15,20,25,30,35)) +
  labs(x = NULL, y = NULL) +
  theme_light()+
  scale_color_gradient(low = "grey", high = "darkblue")+
  theme(legend.text = element_text(size = 16, family = "Helvetica"), legend.title = element_text(size = 22, family = "Helvetica"),
        axis.text.x = element_text(size = 20, family = "Helvetica", angle = 30, hjust = 1),
        axis.text.y =  element_text(size = 18, family = "Helvetica"))+
  guides(fill = guide_legend(override.aes = list(size=10)))

png(file.path(figures_folder,"balloonplot_pathway_hallmarks_onecolor_grad.png"), height = 8000, width = 8000, res = 500)
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
pivot_shapRel <- shap_relevant_stable_notscaled %>% rownames_to_column("circuit") %>% reshape2::melt(.) %>% 
  add_column(KEGG.Pathway = sapply(str_split(.$circuit, ": "), function(x) x[1]) ,.before = "circuit") %>%
  add_column(Circuit.effector = sapply(str_split(.$circuit, ": "), function(x) x[2]) ,.after = "circuit") %>% filter(., .$value != 0)
colnames(pivot_shapRel)<- gsub("variable", "KDT", colnames(pivot_shapRel)) %>%  gsub("value", "SHAP_score", .) ## change the names to be more descriptive

## Export table with drugs and effects
pivot_shapRel_drugbank_functions <- pivot_shapRel %>% merge(.,table_pathways, by.x = "KEGG.Pathway", by.y = "Pathway") %>% merge(., drugbank_effects_tar[,c(10,1,2,3,4,9)], by.x = "KDT", by.y = "symbol")%>% add_column("Module"= table_hallmarks_circuit$Module[match(.$circuit, table_hallmarks_circuit$circuit)])
colnames(pivot_shapRel_drugbank_functions) <- gsub("name", "DRUG", colnames(pivot_shapRel_drugbank_functions))

write.table(pivot_shapRel_drugbank_functions, file =file.path(tables_folder, "ALLpivot_cir_funct_KDT_shapScore_drug_table.tsv"), sep = "\t", quote = F, col.names = T, row.names = F)
saveRDS(pivot_shapRel_drugbank_functions, file.path(rds_folder, "ALLpivot_cir_funct_KDT_shapScore_drug_table.rds"))

length(unique(pivot_shapRel_drugbank_functions$KDT))# 109

## Export only validated 
pivot_shapRel_drugbank_functions_validated <- filter(pivot_shapRel_drugbank_functions, KDT %in% vali)
length(unique(pivot_shapRel_drugbank_functions_validated$KDT)) ## check (6)

write.table(pivot_shapRel_drugbank_functions_validated, file = file.path(tables_folder,"Validated_pivot_cir_funct_KDT_shapScore_drug_table.tsv"), sep = "\t", quote = F, col.names = T, row.names = F)
saveRDS(pivot_shapRel_drugbank_functions_validated, file.path(rds_folder, "Validated_pivot_cir_funct_KDT_shapScore_drug_table.rds"))


#### GET KDTS x Hallmarks (in % of presence) for RADAR, HEATMAP and SPIDER plot ####

## KDTS ##
relevant_KDTbyfunc <- rownames_to_column(hallmarks_circuits_anot, "circuit") %>% merge(., distinct(pivot_shapRel_drugbank_functions[, c("KDT", "circuit")])) %>% .[, -c(1)]  %>% pivot_longer(cols = -KDT) %>%
  mutate(value = case_when(value == "TRUE" ~ 1,
                           value == "FALSE" ~ 0))  %>% group_by(KDT, name) %>%
  summarize(count = sum(value)) %>% .[-which(.$count == 0),] %>%
  pivot_wider(names_from = name , values_from = count ) %>% column_to_rownames("KDT")## Here we transform the DF into a matrix of counts of hallmark per KDT
relevant_KDTbyfunc[is.na(relevant_KDTbyfunc)] <- 0
relevant_KDTbyfunc <- rbind( TOTAL = n_hallmark$total[match(colnames(relevant_KDTbyfunc),n_hallmark$hallmark)] ,relevant_KDTbyfunc) %>% 
  .[, c("Fatty.acid.and.lipid.metabolism", "Apoptosis", "Neuronal", "DNA.integrity" ,"Inflammatory.response", "Stress.response","Necrosis", "Development","Sensory.and.stimuli.transduction")] ## reorder according to hallmark circle

saveRDS(relevant_KDTbyfunc, file = file.path(rds_folder, "relevant_KDTbyfunc.rds"))


## DRUGS ##
drug_byfunc <- rownames_to_column(hallmarks_circuits_anot, "circuit") %>% merge(., distinct(pivot_shapRel_drugbank_functions[,c("circuit", "DRUG")])) %>% .[, -c(1)]  %>% pivot_longer(cols = -DRUG) %>%
  mutate(value = case_when(value == "TRUE" ~ 1,
                           value == "FALSE" ~ 0))  %>% group_by(DRUG, name) %>%
  summarize(count = sum(value)) %>% .[-which(.$count == 0),] %>%
  pivot_wider(names_from = name , values_from = count ) %>% column_to_rownames("DRUG") ## Here we transform the DF into a matrix of counts of hallmark per KDT
drug_byfunc[is.na(drug_byfunc)]<- 0
drug_byfunc <- rbind( TOTAL = n_hallmark$total[match(colnames(drug_byfunc),n_hallmark$hallmark)] ,drug_byfunc) %>%
  .[, c("Fatty.acid.and.lipid.metabolism", "Apoptosis", "Neuronal", "DNA.integrity" ,"Inflammatory.response", "Stress.response","Necrosis", "Development","Sensory.and.stimuli.transduction")] ## reorder according to hallmark circle

saveRDS(drug_byfunc, file = file.path(rds_folder, "drug_byfunc.rds"))


## Hallmarks vs Hallmarks through KDTs ##
ha_ha <- rownames_to_column(relevant_KDTbyfunc, "KDT") %>% merge(., distinct(pivot_shapRel_drugbank_functions[,c("KDT", "Module")])) %>% .[which(!is.na(.$Module)), -c(1)]  %>% pivot_longer(cols = -Module) %>%
  group_by(Module, name)%>%  summarize(count = sum(value))

ha_ha$name <- gsub("\\.", " ", ha_ha$name)

ha_ha <- ha_ha[ which(apply(ha_ha[,c(1,2)], 1, function(x) ifelse(length(unique(x))>1, TRUE, FALSE))),]%>%## Delete the same hallmark interactions
  pivot_wider(names_from = name , values_from = count ) %>% column_to_rownames("Module") %>%## Here we transform the DF into a matrix of counts of hallmark per KDT
  .[, c("Fatty acid and lipid metabolism", "Apoptosis", "Neuronal", "DNA integrity" ,"Inflammatory response", "Stress response","Necrosis", "Development","Sensory and stimuli transduction")] ## reorder according to hallmark circle

ha_ha[is.na(ha_ha)] <- 0
saveRDS(ha_ha, file = file.path(rds_folder, "hallmarks_byhallmarks.rds"))


### DO RADAR PLOT OF PERCENTAGE OF RP MAP COVERED BY EACH HALMARK ###
grid.col <- c("#FFFF33", "#E41A1C", "#984EA3","#4DAF4A", "#FF7F00", "#999999", "#A65628", "#377EB8","#F781BF") ## Set the colors in proper order (like in the spiderplot) from palette Set1
df_radar <- n_hallmark %>% add_column(RP_map = 226) %>% add_column(percent = .$total/226*100)
df_radar$hallmark<- gsub("\\.", " ", df_radar$hallmark) 
rownames(df_radar)<- gsub("\\.", " ", rownames(df_radar)) 
df_radar <- df_radar[colnames(ha_ha),]

df_radar$hallmark<- factor(df_radar$hallmark, levels = rev(colnames(ha_ha)))
grid.col <- setNames(grid.col, colnames(ha_ha))

gra <- ggplot(df_radar, aes(x = hallmark, y= percent, color = hallmark, fill = hallmark))
gra <- gra +
  coord_polar()+
  geom_col(fill = grid.col, color=grid.col) +
  theme_linedraw()+
  theme(legend.position = "none",axis.text.x = element_blank(),
        axis.ticks.y = element_blank(), panel.border = element_blank(),
        axis.text.y.left = element_blank())  +
  # geom_text(data=df_radar, mapping=aes(x=c(1:10), y=max(percent)+max(percent)/5, label= percent),
  #           size=20, angle=0, vjust=0.5, hjust=0.5, color = "gray10")+
  labs(x  =  NULL, y  =  NULL)

png(file.path(figures_folder,"radar_hallmarks_percentageRPmap.png"), height = 8000, width = 8000, res = 500)
gra
dev.off()

