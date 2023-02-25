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

pacman::p_load("here", "stringr", "dplyr","tidyr", "openxlsx", "NMF", "ggnewscale", "ggplot2", "tibble", "gmodels")

#### 1. READ TABLE OF CATEGEORIES FROM ATC  #### Downloaded from https://bioportal.bioontology.org/ontologies/ATC 2022AB CSV file

ATC_raw <- read.delim(here("data", "raw", "ATC.csv"), sep = ",") 
ATC_raw$Class.ID <- str_split(ATC_raw$Class.ID, "ATC/") %>% sapply(., function(x) x[2]) 
ATC_raw <- ATC_raw[which(!is.na(ATC_raw$ATC.LEVEL)), ] %>% .[order(.$Class.ID),]

## Contruct a DF with the first three ATC levels
ATC_df <- data.frame( atc_codes = ATC_raw$Class.ID[ATC_raw$ATC.LEVEL %in% c(1,2,3)], 
                      level = ATC_raw$ATC.LEVEL[ATC_raw$ATC.LEVEL %in% c(1,2,3)],
                      term = ATC_raw$Preferred.Label[ATC_raw$ATC.LEVEL %in% c(1,2,3)],
                      term.synonyms = ATC_raw$Synonyms[ATC_raw$ATC.LEVEL %in% c(1,2,3)])



## Load Drugbank database 
data_folder2 <- "./data/interim"
fname2 <- "drugbank-v050108.tsv"
fpath2 <-file.path(data_folder2,fname2)

drugbank_alltar <- read.delim(file = fpath2, sep = "\t" ) ## raw Database

drugbank_app_action <- drugbank_alltar[-(base::grep("withdrawn",drugbank_alltar$groups)),] %>% .[(base::grep("^approved|approved,", .$groups)),] %>% .[(base::grep("target", .$category)),]%>%
  .[(base::grep("yes", .$known_action)),] %>% .[.$organism == "Humans",]### version 5.1.8 from parser 2022 used to predict


drug_byfunc<- readRDS(file = here("rds", "drug_byfunc.rds")) ## relevant drugs

drugbank_relevants <- drugbank_alltar[-(base::grep("withdrawn",drugbank_alltar$groups)),] %>% .[(base::grep("^approved|approved,", .$groups)),] %>% .[(base::grep("target", .$category)),]%>%
  .[(base::grep("yes", .$known_action)),] %>% .[.$organism == "Humans",] %>% .[which(.$name %in% rownames(drug_byfunc)), ] ## Filtering the relevant drugs from all the assessed

sum(unique(drugbank_relevants$atc_codes) %in% ATC_raw$Class.ID) ## 135 with direct ATC annotations
length(unique(drugbank_relevants$drugbank_id)) ## 284 unique relevant drugs


######## DO IT FOR ALL DRUGS ASSESSED
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

## Try to tag those that don't have an ATC code but are drugs derived from other compounds that do
drugs_ATC$ATC_code[noATC_all] <- drugs_ATC$ATC_code[match(drugs_ATC$synonyms[noATC_all],drugs_ATC$synonyms[-noATC_all])]
sum(is.na(drugs_ATC$ATC_code)) ## 163 now do not have an ATC code
drugs_ATC$ATC_code[is.na(drugs_ATC$ATC_code)] <- "--"

drugs_ATC$ATC_LEVEL1 <-  strsplit(drugs_ATC$ATC_code, "\\|") %>% sapply(., function(x) unique(unlist(substr(x, 1,1)))) ## Level 1 only has the letter
drugs_ATC$ATC_LEVEL1_TERM <- sapply(drugs_ATC$ATC_LEVEL1, function(x) ATC_df$term[match(x, ATC_df$atc_codes)])

drugs_ATC$ATC_LEVEL2 <- strsplit(drugs_ATC$ATC_code, "\\|") %>% sapply(., function(x) unique(unlist(substr(x, 1,3))))  ## Level 2 has the letter plus two apha-numerical
drugs_ATC$ATC_LEVEL2_TERM <- sapply(drugs_ATC$ATC_LEVEL2, function(x) ATC_df$term[match(x, ATC_df$atc_codes)])


#### DO IT FOR ONLY RELEVANT 
## Tag each relevant drug with their category
drugs_ATC_rel <- data.frame(DRUG = drugbank_relevants$name, 
                            synonyms =  sapply(str_split(drugbank_relevants$name, " ") , function(x) x[1]),
                            ID = drugbank_relevants$drugbank_id, 
                            ATC_code =  drugbank_relevants$atc_codes,
                            drugbank_category1 = sapply(str_split(drugbank_relevants$categories, "\\|") , function(x) {x[1]}),
                            drugbank_category2 = sapply(str_split(drugbank_relevants$categories, "\\|") , function(x) {x[2]}),
                            drugbank_categoryLast = sapply(str_split(drugbank_relevants$categories, "\\|") , function(x) {x[length(x)]})) %>% distinct(.)

noATC <- which(drugs_ATC_rel$ATC_code == "") ## 36 Relevant drugs without ATC codes
length(noATC)

## Try to tag those that don't have an ATC code but are drugs derived from other compounds that do
drugs_ATC_rel$ATC_code[noATC] <- drugs_ATC_rel$ATC_code[match(drugs_ATC_rel$synonyms[noATC],drugs_ATC_rel$synonyms[-noATC])]
sum(is.na(drugs_ATC_rel$ATC_code)) ## 6 more drugs have an ATC code
drugs_ATC_rel$ATC_code[is.na(drugs_ATC_rel$ATC_code)] <- "--"


drugs_ATC_rel$ATC_LEVEL1 <-  strsplit(drugs_ATC_rel$ATC_code, "\\|") %>% sapply(., function(x) unique(unlist(substr(x, 1,1)))) ## Level 1 only has the letter
drugs_ATC_rel$ATC_LEVEL1_TERM <- sapply(drugs_ATC_rel$ATC_LEVEL1, function(x) ATC_df$term[match(x, ATC_df$atc_codes)])

drugs_ATC_rel$ATC_LEVEL2 <- strsplit(drugs_ATC_rel$ATC_code, "\\|") %>% sapply(., function(x) unique(unlist(substr(x, 1,3))))  ## Level 2 has the letter plus two apha-numerical
drugs_ATC_rel$ATC_LEVEL2_TERM <- sapply(drugs_ATC_rel$ATC_LEVEL2, function(x) ATC_df$term[match(x, ATC_df$atc_codes)])


### PLOT OF THE DISTRIBUTION OF DRUGS BY ATC LEVEL 1 CLASSIFICATION ###
### PLOT OF THE DISTRIBUTION OF DRUGS BY ATC LEVEL 1 CLASSIFICATION ###

all_leve1_bar <- table(unlist(drugs_ATC$ATC_LEVEL1_TERM))%>% data.frame(.) %>% .[order(.$Freq, decreasing = T),]

colnames(all_leve1_bar) <- c( "ATC", "Freq")
all_leve1_bar$ATC <- factor(all_leve1_bar$ATC, levels = all_leve1_bar$ATC[order(all_leve1_bar$Freq, decreasing = F)])

png(here("results/figures/assessed_drugs_ATC1.png"), height = 5000, width = 10000, res = 500)
ggplot(all_leve1_bar, aes(x= ATC, y = Freq)) +  
  geom_bar(stat = "identity", fill = "magenta") +
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 15, family = "candara", hjust = 1),
        axis.text.x = element_text(size = 15),
        axis.title.x = element_text(size =15))+
  coord_flip()+
  scale_y_continuous(expand = c(0, 0))
dev.off()  


leve1_bar <- table(unlist(drugs_ATC_rel$ATC_LEVEL1_TERM))%>% data.frame(.) %>% .[order(.$Freq, decreasing = T),]
colnames(leve1_bar) <- c( "ATC", "Freq")
leve1_bar$ATC <- factor(leve1_bar$ATC, levels = leve1_bar$ATC[order(leve1_bar$Freq, decreasing = F)])

png(here("results/figures/Relevant_drugs_ATC1.png"), height = 5000, width = 10000, res = 500)
ggplot(leve1_bar, aes(x= ATC, y = Freq)) +  
  geom_bar(stat = "identity", fill = "purple") +
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 15, family = "candara", hjust = 1),
        axis.text.x = element_text(size = 15),
        axis.title.x = element_text(size =15))+
  coord_flip()+
  scale_y_continuous(expand = c(0, 0))
dev.off()  


all_leve2_bar <- table(unlist(drugs_ATC$ATC_LEVEL2_TERM))%>% data.frame(.) %>% .[order(.$Freq, decreasing = T),]
colnames(all_leve2_bar) <- c( "ATC", "Freq")
all_leve2_bar$ATC <- factor(all_leve2_bar$ATC, levels = all_leve2_bar$ATC[order(all_leve2_bar$Freq, decreasing = F)])

png(here("results/figures/assessed_drugs_ATC2.png"), height = 5000, width = 10000, res = 500)
ggplot(all_leve2_bar, aes(x= ATC, y = Freq)) +  
  geom_bar(stat = "identity", fill = "magenta") +
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 8, family = "candara", hjust = 1),
        axis.text.x = element_text(size = 15),
        axis.title.x = element_text(size =15))+
  coord_flip()+
  scale_y_continuous(expand = c(0, 0))
dev.off() 

leve2_bar <- table(unlist(drugs_ATC_rel$ATC_LEVEL2_TERM))%>% data.frame(.) %>% .[order(.$Freq, decreasing = T),]
colnames(leve2_bar) <- c( "ATC", "Freq")
leve2_bar$ATC <- factor(leve2_bar$ATC, levels = leve2_bar$ATC[order(leve2_bar$Freq, decreasing = F)])

png(here("results/figures/Relevant_drugs_ATC2.png"), height = 5000, width = 10000, res = 500)
ggplot(leve2_bar, aes(x= ATC, y = Freq)) +  
  geom_bar(stat = "identity", fill = "purple") +
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 8, family = "candara", hjust = 1),
        axis.text.x = element_text(size = 15),
        axis.title.x = element_text(size =15))+
  coord_flip()+
  scale_y_continuous(expand = c(0, 0))
dev.off() 


#################################################################
#### OVER REPRESENTATION ANALYASIS OF MY DRUGS CATEGORIES 
####### with Fisher's Exact test
################################################

data_fisherTest <- tidyr::unnest(drugs_ATC[, c(1,11)], cols = c(ATC_LEVEL2_TERM))  %>% add_column("Relevant" = .$DRUG %in% drugs_ATC_rel$DRUG)

length(data_fisherTest$DRUG)
length(unique(data_fisherTest$DRUG))

fisher.test(data_fisherTest$ATC_LEVEL2_TERM, data_fisherTest$Relevant, alternative = "greater")


CrT <- gmodels:: CrossTable(data_fisherTest$ATC_LEVEL2_TERM, data_fisherTest$Relevant , resid = T, prop.chisq = TRUE, fisher = T, simulate.p.value=TRUE)


# l <- split(data_fisherTest, data_fisherTest$ATC_LEVEL2_TERM)
# lapply(l, function(x) CrossTable(x$ATC_LEVEL2_TERM, x$Relevant, 
#                                  prop.chisq = FALSE, fisher = TRUE))
