############################################
## Project: RP DRexML
## Script purpose:  Spider plot of relevan KDTs around hallmark annotations.
## Date: 25.02.2023
## Author: Marina Esteban-Medina
#######################################

library("here")
library("stringr")
library("dplyr")
library("tidyr")
library("openxlsx")
library("NMF")
library("ggnewscale")
library("ggplot2")
library("fmsb")
library("tibble")


if(!dir.exists(here("results/tables"))){
  dir.create(here("results/tables"))
}

if(!dir.exists(here("results/figures"))){
  dir.create(here("results/figures"))
}

if(!dir.exists(here("rds"))){
  dir.create(here("rds"))
}

#### 1.Read the "rds" files for de KDT_by function and DRUG_byfunction ####

relevant_KDTbyfunc <- readRDS(file = here("rds", "relevant_KDTbyfunc.rds"))

drug_byfunc<- readRDS(file = here("rds", "drug_byfunc.rds"))

## Transform count of hallmarks in percentage of the hallmark represented in a KDT or drug.
relevant_KDTbyfunc_percent <- apply(relevant_KDTbyfunc, 2, function(x) round(x/max(x)*100)) %>% data.frame(.)
drug_byfunc_percent <- apply(drug_byfunc, 2, function(x) round(x/max(x)*100)) %>% data.frame(.) 

#### 2. SPIDER PLOT  #####

# Normalize the data (Has 0-1 range) for the radar plot of whatever we would like to plot
normalize <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

drugs_fromValidated <- c("Omega-3-carboxylic acids", "Zopiclone", "Taurine")

drugs_barbiturics <- c("Phenobarbital", "Talbutal", "Butobarbital")

drugs_compare <- drug_byfunc_percent[drugs_fromValidated, ] %>% apply(., 1, normalize) %>% t(.) %>% data.frame(.)
# drugs_compare <- drug_byfunc_percent[drugs_barbiturics, ] %>% apply(., 1, normalize) %>% t(.) %>% data.frame(.) 
drugs_compare <- rbind(rep(max(drugs_compare), 9), rep(0, 9), drugs_compare)

# Define fill colors
colors_fill <- c(scales::alpha("yellow", 0.1),
                 scales::alpha("lightgreen", 0.2),
                 scales::alpha("purple", 0.2))

# Define line colors
colors_line <-  c(scales::alpha("gold", 0.9),
                  scales::alpha("darkgreen", 0.9),
                  scales::alpha("purple", 0.9))


png(here("results", "figures","spider_Omega3_Zopiclone_Taurine.png"), height = 6000, width = 9000, res = 500)
# png(here("results", "figures","spider_barbiturics.png"), height = 6000, width = 9000, res = 500)
radarchart(drugs_compare, 
           seg = 9,  # Number of axis segments
           title = "Drugs targeting validated KDTs: ELOVL4, GABRA1, GLRA2",
           pcol = colors_line,
           pfcol = colors_fill,
           plwd = 3
)

# Add a legend
legend(x=1.15, 
       y=1.20, 
       legend = rownames(drugs_compare[-c(1,2),]), 
       bty = "n", pch=20 , col = colors_line, cex = 1.5, pt.cex = 3.5,  )
dev.off()

