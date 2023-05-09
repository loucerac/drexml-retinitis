############################################
## Project: RP DRexML
## Script purpose: Extrac ATC categories level 1 and level 2 to tag relevant resulting drugs vs all Drugbank DB
## Date: 25.02.2023
## Author: Marina Esteban-Medina
#######################################


if (!require(pacman, quietly = TRUE)){
  install.packages(pacman)
}else{
  library(pacman)
}

pacman::p_load("here", "stringr", "dplyr","tidyr", "openxlsx", "NMF", "ggnewscale", "ggplot2", "tibble", "gmodels", "RColorBrewer")

#### 1. READ TABLE OF CATEGEORIES FROM ATC  #### Downloaded from https://bioportal.bioontology.org/ontologies/ATC 2022AB CSV file

ATC_raw <- read.delim(here("data", "raw", "ATC.csv"), sep = ",") 
ATC_raw$Class.ID <- str_split(ATC_raw$Class.ID, "ATC/") %>% sapply(., function(x) x[2]) 
ATC_raw <- ATC_raw[which(!is.na(ATC_raw$ATC.LEVEL)), ] %>% .[order(.$Class.ID),]

## Contruct a DF with the first three ATC levels
ATC_df <- data.frame( atc_codes = ATC_raw$Class.ID[ATC_raw$ATC.LEVEL %in% c(1,2,3,4)], 
                      level = ATC_raw$ATC.LEVEL[ATC_raw$ATC.LEVEL %in% c(1,2,3,4)],
                      term = ATC_raw$Preferred.Label[ATC_raw$ATC.LEVEL %in% c(1,2,3,4)],
                      term.synonyms = ATC_raw$Synonyms[ATC_raw$ATC.LEVEL %in% c(1,2,3,4)])

## Load Drugbank database 
data_folder2 <- "./data/interim"
fname2 <- "drugbank-v050108.tsv"
fpath2 <-file.path(data_folder2,fname2)

drugbank_alltar <- read.delim(file = fpath2, sep = "\t" ) ## raw Database

drugbank_app_action <- drugbank_alltar[-(base::grep("withdrawn",drugbank_alltar$groups)),] %>% .[(base::grep("^approved|approved,", .$groups)),] %>% .[(base::grep("target", .$category)),]%>%
  .[(base::grep("yes", .$known_action)),] %>% .[.$organism == "Humans",]### version 5.1.8 from parser 2022 used to predict


drug_byfunc<- readRDS(file = here("rds", "drug_byfunc.rds")) ## relevant drugs

drugbank_relevants <- drugbank_app_action %>% .[which(.$name %in% rownames(drug_byfunc)), ] ## Filtering the relevant drugs from all the assessed

sum(unique(drugbank_relevants$atc_codes) %in% ATC_raw$Class.ID) ## 135 with direct ATC annotations
length(unique(drugbank_relevants$drugbank_id)) ## 284 unique relevant drugs

######## DO IT FOR ALL DRUGS ASSESSED #####
## Tag all drugs assessed with the ATC codes nomenclature
drugs_ATC <- data.frame(DRUG = drugbank_app_action$name, 
                        synonyms =  sapply(str_split(drugbank_app_action$name, " ") , function(x) x[1]),
                        ID = drugbank_app_action$drugbank_id, 
                        ATC_code =  drugbank_app_action$atc_codes,
                        drugbank_category1 = sapply(str_split(drugbank_app_action$categories, "\\|") , function(x) {x[1]}),
                        drugbank_category2 = sapply(str_split(drugbank_app_action$categories, "\\|") , function(x) {x[2]}),
                        drugbank_categoryLast = sapply(str_split(drugbank_app_action$categories, "\\|") , function(x) {x[length(x)]})) %>% distinct(.)

noATC_all <- which(drugs_ATC$ATC_code == "") ## 36 Relevant drugs without ATC codes
length(noATC_all)#194
drugs_ATC$ATC_code[is.na(drugs_ATC$ATC_code)] <- "--"


## Try to tag those that don't have an ATC code but are drugs derived from other compounds that do
drugs_ATC$ATC_code[noATC_all] <- drugs_ATC$ATC_code[match(drugs_ATC$synonyms[noATC_all],drugs_ATC$synonyms[-noATC_all])]
sum(is.na(drugs_ATC$ATC_code)) ## 163 now do not have an ATC code
drugs_ATC$ATC_code[is.na(drugs_ATC$ATC_code)] <- "--"

drugs_ATC$ATC_LEVEL1 <-  strsplit(drugs_ATC$ATC_code, "\\|") %>% sapply(., function(x) unique(unlist(substr(x, 1,1)))) ## Level 1 only has the letter
drugs_ATC$ATC_LEVEL1_TERM <- sapply(drugs_ATC$ATC_LEVEL1, function(x) ATC_df$term[match(x, ATC_df$atc_codes)])

drugs_ATC$ATC_LEVEL2 <- strsplit(drugs_ATC$ATC_code, "\\|") %>% sapply(., function(x) unique(unlist(substr(x, 1,3))))  ## Level 2 has the letter plus two apha-numerical
drugs_ATC$ATC_LEVEL2_TERM <- sapply(drugs_ATC$ATC_LEVEL2, function(x) ATC_df$term[match(x, ATC_df$atc_codes)])

drugs_ATC$ATC_LEVEL3 <- strsplit(drugs_ATC$ATC_code, "\\|") %>% sapply(., function(x) unique(unlist(substr(x, 1,4))))  ## Level 3 has the letter plus two numerical plus a letter
drugs_ATC$ATC_LEVEL3_TERM <- sapply(drugs_ATC$ATC_LEVEL3, function(x) ATC_df$term[match(x, ATC_df$atc_codes)])

drugs_ATC$ATC_LEVEL4 <- strsplit(drugs_ATC$ATC_code, "\\|") %>% sapply(., function(x) unique(unlist(substr(x, 1,5))))  ## Level 4 has the letter plus two numerical plus a letter plus two digits
drugs_ATC$ATC_LEVEL4_TERM <- sapply(drugs_ATC$ATC_LEVEL4, function(x) ATC_df$term[match(x, ATC_df$atc_codes)])

#### DO IT FOR ONLY RELEVANT DRUGS ####
## Tag each relevant drug with their category
drugs_ATC_rel <- data.frame(DRUG = drugbank_relevants$name, 
                            synonyms =  sapply(str_split(drugbank_relevants$name, " ") , function(x) x[1]),
                            ID = drugbank_relevants$drugbank_id, 
                            ATC_code =  drugbank_relevants$atc_codes,
                            drugbank_category1 = sapply(str_split(drugbank_relevants$categories, "\\|") , function(x) {x[1]}),
                            drugbank_category2 = sapply(str_split(drugbank_relevants$categories, "\\|") , function(x) {x[2]}),
                            drugbank_categoryLast = sapply(str_split(drugbank_relevants$categories, "\\|") , function(x) {x[length(x)]})) %>% distinct(.)

noATC <- which(drugs_ATC_rel$ATC_code == "") ## 36 Relevant drugs without ATC codes
length(noATC) ## 36 Relevant drugs without ATC codes
drugs_ATC_rel$ATC_code[is.na(drugs_ATC_rel$ATC_code)] <- "--"

drugs_ATC_rel$ATC_LEVEL1 <-  strsplit(drugs_ATC_rel$ATC_code, "\\|") %>% sapply(., function(x) unique(unlist(substr(x, 1,1)))) ## Level 1 only has the letter
drugs_ATC_rel$ATC_LEVEL1_TERM <- sapply(drugs_ATC_rel$ATC_LEVEL1, function(x) ATC_df$term[match(x, ATC_df$atc_codes)])

drugs_ATC_rel$ATC_LEVEL2 <- strsplit(drugs_ATC_rel$ATC_code, "\\|") %>% sapply(., function(x) unique(unlist(substr(x, 1,3))))  ## Level 2 has the letter plus two apha-numerical
drugs_ATC_rel$ATC_LEVEL2_TERM <- sapply(drugs_ATC_rel$ATC_LEVEL2, function(x) ATC_df$term[match(x, ATC_df$atc_codes)])

drugs_ATC_rel$ATC_LEVEL3 <- strsplit(drugs_ATC_rel$ATC_code, "\\|") %>% sapply(., function(x) unique(unlist(substr(x, 1,4))))  ## Level 3 has the letter plus two numerical plus a letter
drugs_ATC_rel$ATC_LEVEL3_TERM <- sapply(drugs_ATC_rel$ATC_LEVEL3, function(x) ATC_df$term[match(x, ATC_df$atc_codes)])

drugs_ATC_rel$ATC_LEVEL4 <- strsplit(drugs_ATC_rel$ATC_code, "\\|") %>% sapply(., function(x) unique(unlist(substr(x, 1,5))))  ## Level 4 has the letter plus two numerical plus a letter plus two digits
drugs_ATC_rel$ATC_LEVEL4_TERM <- sapply(drugs_ATC_rel$ATC_LEVEL4, function(x) ATC_df$term[match(x, ATC_df$atc_codes)])

### Table per level of ATC category ###

## Level 1 ##
all_leve1_bar <- table(unlist(drugs_ATC$ATC_LEVEL1_TERM))%>% data.frame(.) %>% .[order(.$Freq, decreasing = T),]
colnames(all_leve1_bar) <- c( "ATC", "Freq")
all_leve1_bar$ATC <- factor(all_leve1_bar$ATC, levels = all_leve1_bar$ATC[order(all_leve1_bar$Freq, decreasing = F)])
all_leve1_bar$ATC_code <- ATC_df$atc_codes[match(all_leve1_bar$ATC,ATC_df$term)]

leve1_bar <- table(unlist(drugs_ATC_rel$ATC_LEVEL1_TERM))%>% data.frame(.) %>% .[order(.$Freq, decreasing = T),]
colnames(leve1_bar) <- c( "ATC", "Freq")
leve1_bar$ATC <- factor(leve1_bar$ATC, levels = leve1_bar$ATC[order(leve1_bar$Freq, decreasing = F)])
leve1_bar$ATC_code <- ATC_df$atc_codes[match(leve1_bar$ATC,ATC_df$term)]

## Level 2 ##
all_leve2_bar <- table(unlist(drugs_ATC$ATC_LEVEL2_TERM))%>% data.frame(.) %>% .[order(.$Freq, decreasing = T),]
colnames(all_leve2_bar) <- c( "ATC", "Freq")
all_leve2_bar$ATC <- factor(all_leve2_bar$ATC, levels = all_leve2_bar$ATC[order(all_leve2_bar$Freq, decreasing = F)])
all_leve2_bar$ATC_code <- ATC_df$atc_codes[match(all_leve2_bar$ATC,ATC_df$term)]

leve2_bar <- table(unlist(drugs_ATC_rel$ATC_LEVEL2_TERM))%>% data.frame(.) %>% .[order(.$Freq, decreasing = T),]
colnames(leve2_bar) <- c( "ATC", "Freq")
leve2_bar$ATC <- factor(leve2_bar$ATC, levels = leve2_bar$ATC[order(leve2_bar$Freq, decreasing = F)])
leve2_bar$ATC_code <- ATC_df$atc_codes[match(leve2_bar$ATC,ATC_df$term)]

## Level 3 ##
all_leve3_bar <- table(unlist(drugs_ATC$ATC_LEVEL3_TERM))%>% data.frame(.) %>% .[order(.$Freq, decreasing = T),]
colnames(all_leve3_bar) <- c( "ATC", "Freq")
all_leve3_bar$ATC <- factor(all_leve3_bar$ATC, levels = all_leve3_bar$ATC[order(all_leve3_bar$Freq, decreasing = F)])
all_leve3_bar$ATC_code <- ATC_df$atc_codes[match(all_leve3_bar$ATC,ATC_df$term)]

leve3_bar <- table(unlist(drugs_ATC_rel$ATC_LEVEL3_TERM))%>% data.frame(.) %>% .[order(.$Freq, decreasing = T),]
colnames(leve3_bar) <- c( "ATC", "Freq")
leve3_bar$ATC <- factor(leve3_bar$ATC, levels = leve3_bar$ATC[order(leve3_bar$Freq, decreasing = F)])
leve3_bar$ATC_code <- ATC_df$atc_codes[match(leve3_bar$ATC,ATC_df$term)]


## Level 4 ##
all_leve4_bar <- table(unlist(drugs_ATC$ATC_LEVEL4_TERM))%>% data.frame(.) %>% .[order(.$Freq, decreasing = T),]
colnames(all_leve4_bar) <- c( "ATC", "Freq")
all_leve4_bar$ATC <- factor(all_leve4_bar$ATC, levels = all_leve4_bar$ATC[order(all_leve4_bar$Freq, decreasing = F)])
all_leve4_bar$ATC_code <- ATC_df$atc_codes[match(all_leve4_bar$ATC,ATC_df$term)]


leve4_bar <- table(unlist(drugs_ATC_rel$ATC_LEVEL4_TERM))%>% data.frame(.) %>% .[order(.$Freq, decreasing = T),]
colnames(leve4_bar) <- c( "ATC", "Freq")
leve4_bar$ATC <- factor(leve4_bar$ATC, levels = leve4_bar$ATC[order(leve4_bar$Freq, decreasing = F)])
leve4_bar$ATC_code <- ATC_df$atc_codes[match(leve4_bar$ATC,ATC_df$term)]


#################################################################
#### OVER REPRESENTATION ANALYASIS OF MY DRUGS CATEGORIES 
####### with Fisher's Exact test
################################################

data_fisherTest_level4 <- tidyr::unnest(drugs_ATC[, c("DRUG","ATC_LEVEL4")], cols = c(ATC_LEVEL4))  %>% add_column("Relevant" = .$DRUG %in% drugs_ATC_rel$DRUG) 
data_fisherTest_level3 <- tidyr::unnest(drugs_ATC[, c("DRUG","ATC_LEVEL3")], cols = c(ATC_LEVEL3))  %>% add_column("Relevant" = .$DRUG %in% drugs_ATC_rel$DRUG)
data_fisherTest_level2 <- tidyr::unnest(drugs_ATC[, c("DRUG","ATC_LEVEL2")], cols = c(ATC_LEVEL2))  %>% add_column("Relevant" = .$DRUG %in% drugs_ATC_rel$DRUG)
data_fisherTest_level1 <- tidyr::unnest(drugs_ATC[, c("DRUG","ATC_LEVEL1")], cols = c(ATC_LEVEL1))  %>% add_column("Relevant" = .$DRUG %in% drugs_ATC_rel$DRUG)

## Load ORA results from ATC levels 1,2,3,4 ###
ora_ATC <- read.delim(file =here("results", "tables","selected_drugs_atc_ora.tsv")) %>% filter(ora_bylevel_pval_adj<0.05) %>% 
  .[order(.$ ora_bylevel_pval_adj, decreasing = F ), ] %>% .[order(.$atc_level, decreasing = F), ]

### STACKED BARPLOTS OVEREPRESENTED FOR CATEGORIES 2, 3, 4 ####
## ATC 4 ###
ATC_4_together <- all_leve4_bar[all_leve4_bar$ATC_code %in% ora_ATC$atc_code, c("ATC", "Freq")] %>% add_column("Freq_relevants" = leve4_bar$Freq[match(.$ATC, leve4_bar$ATC)]) %>% .[which(!is.na(.$Freq_relevants)),]
max_lim4 <- max(ATC_4_together$Freq)
ATC4_rel_percentage <- data.frame(ATC = ATC_4_together$ATC , percent = paste0(round(x = ATC_4_together$Freq_relevants/ATC_4_together$Freq *100, digits = 0), " %"))## Create another Df to tag the percentage

ATC_4_together$Freq <- ATC_4_together$Freq - ATC_4_together$Freq_relevants  ## For all drugbank db drugs rest the ones that are relevant to have a correct length of the stacked plot
                                                                               ## showing the total of drugs when both plots are stacked
colnames(ATC_4_together) <- c("ATC_level4", "ALL DRUGBANK", "RELEVANT DRUGS")

ATC_4_together <- pivot_longer(ATC_4_together, cols = -ATC_level4 ) %>% 
  add_column(percentage = ATC4_rel_percentage$percent[ match(.$ATC_level4, ATC4_rel_percentage$ATC)]) 

ATC_4_together$percentage[which(ATC_4_together$name == "ALL DRUGBANK")] <- NA ## Delete labels of ALL and 0% to clean the plot

ATC_4_together <- ATC_4_together[-which(ATC_4_together$ATC_level4 %in% c("Antipropulsives")),] ## REMOVE

# png(here("results/figures/Relevant_ATC4_stacked.png"), height = 16000, width = 8000, res = 500)## FOR ALL RELEVANT
png(here("results/figures/ORA_Relevant_ATC4_stacked.png"), height = 3200, width = 14000, res = 500)## FOR ORA RELEVANT
# Stacked barplot with multiple groups
ggplot(data=ATC_4_together, aes(x= ATC_level4, y= value, fill= name))+
  geom_bar(stat="identity")+
  geom_text(aes(label = percentage), hjust = -0.1, size = 8, family = "candara", colour = "#040f42" )+
  theme_minimal()+
  theme(legend.title =  element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 25, family = "candara", hjust = 1),
        axis.text.x = element_text(size = 30),
        axis.title.x = element_blank())+
  coord_flip()+
  scale_y_continuous(expand = c(0, 0), limits = c(0,max_lim4), breaks = seq(0,max_lim4, 5))+
  # scale_x_discrete()+
  # scale_fill_brewer(palette="Paired")+
  scale_fill_manual(values = c("#FF99FF","#984EA3"))
dev.off()

## ATC 3 ###
ATC_3_together <- all_leve3_bar[all_leve3_bar$ATC_code %in% ora_ATC$atc_code, c("ATC", "Freq")] %>% add_column("Freq_relevants" = leve3_bar$Freq[match(.$ATC, leve3_bar$ATC)]) %>% .[which(!is.na(.$Freq_relevants)),]
max_lim3 <- max(ATC_3_together$Freq)
ATC3_rel_percentage <- data.frame(ATC = ATC_3_together$ATC , percent = paste0(round(x = ATC_3_together$Freq_relevants/ATC_3_together$Freq *100, digits = 0), " %"))## Create another Df to tag the percentage

ATC_3_together$Freq <- ATC_3_together$Freq - ATC_3_together$Freq_relevants  ## For all drugbank db drugs rest the ones that are relevant to have a correct length of the stacked plot
## showing the total of drugs when both plots are stacked
colnames(ATC_3_together) <- c("ATC_level3", "ALL DRUGBANK", "RELEVANT DRUGS")

ATC_3_together <- pivot_longer(ATC_3_together, cols = -ATC_level3 ) %>% 
  add_column(percentage = ATC3_rel_percentage$percent[ match(.$ATC_level3, ATC3_rel_percentage$ATC)]) 

ATC_3_together$percentage[which(ATC_3_together$name == "ALL DRUGBANK")] <- NA ## Delete labels of ALL and 0% to clean the plot

ATC_3_together <- ATC_3_together[-which(ATC_3_together$ATC_level3 %in% c("ANTIPROPULSIVES", "DRUGS FOR CONSTIPATION")),] ## REMOVE

png(here("results/figures/ORA_Relevant_ATC3_stacked.png"), height = 3000, width = 14000, res = 500)
# Stacked barplot with multiple groups
ggplot(data=ATC_3_together, aes(x= ATC_level3, y= value, fill= name))+
  geom_bar(stat="identity")+
  geom_text(aes(label = percentage), hjust = -0.1, size = 8, family = "candara", colour = "#040f42" )+
  theme_minimal()+
  theme(legend.title =  element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size =24, family = "candara", hjust = 1),
        axis.text.x = element_text(size = 30),
        axis.title.x = element_blank())+
  coord_flip()+
  scale_y_continuous(expand = c(0, 0), limits = c(0,max_lim3), breaks = seq(0,max_lim3, 20))+
  scale_fill_manual(values = c("#FF99FF","#984EA3"))
dev.off()

## ATC 2 ###
ATC_2_together <- all_leve2_bar[all_leve2_bar$ATC_code %in% ora_ATC$atc_code, c("ATC", "Freq")] %>% add_column("Freq_relevants" = leve2_bar$Freq[match(.$ATC, leve2_bar$ATC)]) %>% .[which(!is.na(.$Freq_relevants)),]
max_lim2 <- max(ATC_2_together$Freq)
ATC2_rel_percentage <- data.frame(ATC = ATC_2_together$ATC , percent = paste0(round(x = ATC_2_together$Freq_relevants/ATC_2_together$Freq *100, digits = 0), " %"))## Create another Df to tag the percentage

ATC_2_together$Freq <- ATC_2_together$Freq - ATC_2_together$Freq_relevants ## For all drugbank db drugs rest the ones that are relevant to have a correct length of the stacked plot showing the total of drugs when both plots are stacked
colnames(ATC_2_together) <- c("ATC_level2", "ALL DRUGBANK", "RELEVANT DRUGS")
ATC_2_together <- pivot_longer(ATC_2_together, cols = -ATC_level2 ) %>% add_column(percentage = ATC2_rel_percentage$percent[ match(.$ATC_level2, ATC2_rel_percentage$ATC)] )
ATC_2_together$percentage[which(ATC_2_together$name == "ALL DRUGBANK") ] <- NA ## Delete labels of ALL and 0% to clean the plot
ATC_2_together <- ATC_2_together[-which(ATC_2_together$ATC_level2 %in% c("ANTIPROPULSIVES", "DRUGS FOR CONSTIPATION")),] ## REMOVE


png(here("results/figures/ORA_Relevant_ATC2_stacked.png"), height = 2000, width = 14000, res = 500)
# Stacked barplot with multiple groups
ggplot(data=ATC_2_together, aes(x= ATC_level2, y= value, fill= name))+
  geom_bar(stat="identity")+
  geom_text(aes(label = percentage), hjust = -0.1, size = 8, family = "candara", colour = "#040f42" )+
  theme_minimal()+
  theme(legend.title =  element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 24, family = "candara", hjust = 1),
        axis.text.x = element_text(size = 30),
        axis.title.x = element_blank())+
  coord_flip()+
  # scale_y_continuous(expand = c(0, 0), limits = c(0, 150))+
  # scale_x_discrete()+
  # scale_fill_brewer(palette="Paired")+
  scale_fill_manual(values = c("#FF99FF","#984EA3"))
dev.off()


## ATC 1 ###
ATC_1_together <- all_leve1_bar[all_leve1_bar$ATC_code %in% ora_ATC$atc_code, c("ATC", "Freq")] %>% add_column("Freq_relevants" = leve1_bar$Freq[match(.$ATC, leve1_bar$ATC)]) %>% .[which(!is.na(.$Freq_relevants)),]
max_lim1 <- max(ATC_1_together$Freq)
ATC1_rel_percentage <- data.frame(ATC = ATC_1_together$ATC , percent = paste0(round(x = ATC_1_together$Freq_relevants/ATC_1_together$Freq *100, digits = 0), " %"))## Create another Df to tag the percentage

ATC_1_together$Freq <- ATC_1_together$Freq - ATC_1_together$Freq_relevants ## For all drugbank db drugs rest the ones that are relevant to have a correct length of the stacked plot showing the total of drugs when both plots are stacked
colnames(ATC_1_together) <- c("ATC_level1", "ALL DRUGBANK", "RELEVANT DRUGS")
ATC_1_together <- pivot_longer(ATC_1_together, cols = -ATC_level1 ) %>% add_column(percentage = ATC1_rel_percentage$percent[ match(.$ATC_level1, ATC1_rel_percentage$ATC)] )
ATC_1_together$percentage[which(ATC_1_together$name == "ALL DRUGBANK") ] <- NA ## Delete labels of ALL and 0% to clean the plot

png(here("results/figures/ORA_Relevant_ATC1_stacked.png"), height = 1000, width = 10000, res = 500)
# Stacked barplot with multiple groups
ggplot(data= ATC_1_together, aes(x= ATC_level1, y= value, fill= name))+
  geom_bar(stat="identity")+
  geom_text(aes(label = percentage), hjust = -0.1, size =8, family = "candara", colour = "#040f42" )+
  theme_minimal()+
  theme(legend.title =  element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 24, family = "candara", hjust = 1),
        axis.text.x = element_text(size = 30),
        axis.title.x = element_blank())+
  coord_flip()+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 250))+
  # scale_x_discrete()+
  # scale_fill_brewer(palette="Paired")+
  scale_fill_manual(values = c("#FF99FF","#984EA3"))
dev.off()

####################################################################################
#### ADD OVER REPRESENTATED DRUGS TO THE TABLES THAT CONTAINED RELEVANT KDTS #####

## Add the ATC_level 1 and 2 terms to the big table
pivot_shapRel_drugbank_functions <- readRDS(here("rds", "ALLpivot_cir_funct_KDT_shapScore_drug_table.rds"))
pivot_shapRel_drugbank_functions_ATC_1_2_3_4 <- merge(pivot_shapRel_drugbank_functions, data_fisherTest_level1[which(data_fisherTest_level1$Relevant == T) , c("DRUG","ATC_LEVEL1")]) %>%  
  merge(., data_fisherTest_level2[which(data_fisherTest_level2$Relevant == T) , c("DRUG","ATC_LEVEL2")]) %>%
  merge(., data_fisherTest_level3[which(data_fisherTest_level3$Relevant == T) , c("DRUG","ATC_LEVEL3")]) %>%
  merge(., data_fisherTest_level4[which(data_fisherTest_level4$Relevant == T) , c("DRUG","ATC_LEVEL4")])


## READ MANUALLY CURATED HALLMARKS 
hallmarks_circuits_anot <-  read.xlsx(xlsxFile =  "./data/final/RP_map_functions_MPC-annot.xlsx") %>% .[,-c(1,3,13,14,15)] %>% column_to_rownames(., "circuit_name")%>%
  apply(., 2, function(x) ifelse(is.na(x), F, T)) %>% as.data.frame(.)

## Here I am recounting to see how many hallmarks we have summarized by PATHWAY
n_cir <- data.frame(table(sapply(str_split(rownames(hallmarks_circuits_anot), ": "), function(x) x[1]))) # get the nº of circuit per pathway to obtain normalized Freq of Hallmark per pathway
n_hallmark <- data.frame(hallmark = colnames(hallmarks_circuits_anot), total = apply(hallmarks_circuits_anot, 2, sum))

## Add the info to create a table of counts of ATC ORA relevant categories from level 4,3,2,1 ##

## ATC_ ORA RELEVANT BY HALLMARKS ###
ATCora4_byfunc <- rownames_to_column(hallmarks_circuits_anot, "circuit") %>% merge(., distinct(pivot_shapRel_drugbank_functions_ATC_1_2_3_4[,c("circuit","ATC_LEVEL4")])) %>% .[, -c(1)]  %>% pivot_longer(cols = -ATC_LEVEL4) %>%
  mutate(value = case_when(value == "TRUE" ~ 1,
                           value == "FALSE" ~ 0))  %>% group_by(ATC_LEVEL4, name) %>%
  summarize(count = sum(value)) %>% .[-which(.$count == 0),] %>% ## Make it a DF of counts per hallmark
  pivot_wider(names_from = name , values_from = count ) %>%  ## Here we transform the DF into a matrix of counts of hallmark per KDT
  filter(ATC_LEVEL4 %in% ora_ATC$atc_code[ora_ATC$atc_level == 4]) %>% ## Filter the ones with ORA relevant FRD.pval
  add_column("ATC_LEVEL4_TERM"= ora_ATC$atc_name[match(.$ATC_LEVEL4, ora_ATC$atc_code)], .before = "ATC_LEVEL4") %>% .[,which(colnames(.)!= "ATC_LEVEL4")] 
ATCora4_byfunc$ATC_LEVEL4_TERM[which(ATCora4_byfunc$ATC_LEVEL4_TERM =="")] <- "Protein kinase inhibitors" ## Manually curated from WHOCC at https://www.whocc.no/news/proposed_new_classification_for_l01xe_protein_kinase_inhibit

ATCora4_byfunc <- column_to_rownames(.data = ATCora4_byfunc,"ATC_LEVEL4_TERM")
ATCora4_byfunc[is.na(ATCora4_byfunc)] <- 0
ATCora4_byfunc <- rbind( TOTAL = n_hallmark$total[match(colnames(ATCora4_byfunc),n_hallmark$hallmark)] ,ATCora4_byfunc) %>%
  .[, c("Fatty.acid.and.lipid.metabolism", "Apoptosis", "Neuronal", "DNA.integrity" ,"Inflammatory.response", "Stress.response","Necrosis", "Development","Sensory.and.stimuli.transduction")] ## reorder according to hallmark circle

saveRDS(ATCora4_byfunc, file = here("rds", "ATCora4_byfunc.rds"))

ATCora4_byfunc_percent <- apply(ATCora4_byfunc, 2, function(x) round(x/max(x)*100)) %>% data.frame(.)
colnames(ATCora4_byfunc_percent)<- gsub("\\.", " ", colnames(ATCora4_byfunc_percent))

png(filename = here("results/figures/ATC_level4_ora_byhallmark_perecent.png"),height = 8000 ,width = 8000, res = 500)
# NMF::aheatmap(t(ATCora4_byfunc_percent[-1,]), color = "-heat",  border_color = "white", Rowv = F, Colv = F, fontsize = 15)
pheatmap::pheatmap(ATCora4_byfunc_percent[-1,], angle_col = 45, color = rev(colorRampPalette(rev(brewer.pal(n = 7, name ="Purples")))(20)), fontsize = 15)
dev.off()
