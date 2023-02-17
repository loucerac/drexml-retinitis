#!/usr/bin/env RScript

#########################################
### Processing GTEx V8 datasets #####
#########################################

library(here)
library(hipathia)
library(feather)
library(edgeR)
library(data.table)

dir.create(here("data", "interim"), showWarnings = FALSE, recursive=TRUE)
dir.create(here("data", "final"), showWarnings = FALSE, recursive=TRUE)


# Load  RNA-seq expression set from GTEx *.gct file.

## Read the downloaded GTEx raw counts dataset
expreset_raw <- fread(
  file = here("data", "raw", "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct"),
  header = T, sep = "\t"
) %>% as.data.frame(.)

rownames(expreset_raw) <- expreset_raw$Name
expreset_raw <- expreset_raw[, -(1:2)]
print("read...done")

# Normalization by TMM with "edgeR" package
dge <- DGEList(counts = expreset_raw)
print("dge...done")
tmm <- calcNormFactors(dge, method = "TMM")
print("tmm...done")
logcpm <- cpm(tmm, prior.count = 3, log = TRUE)
print("dge...done")

# eliminate from rownames the ".number", beacuse Hipathia does not process them well
rownames(logcpm) <- gsub("\\..*", "", rownames(logcpm))
print("normalization 1...done")


###### Hipathia Processing #####

trans_data <- translate_data(logcpm, "hsa")
exp_data <- normalize_data(trans_data)

## Loading Pathways (only physiological)

# physiological_pathways
path_list <- read.table(file = here("data", "raw", "physiological_paths.tsv"), sep = "\t")
pathways <- load_pathways("hsa", pathways_list = path_list[, 2])

## Using Hipathia to compute the signal

results <- hipathia(exp_data, pathways, decompose = FALSE, verbose = FALSE)
path_vals <- get_paths_data(results, matrix = TRUE)


saveRDS(exp_data, file = here("data", "final", "expreset_Hinorm_gtexV8.rds"))

saveRDS(path_vals, file = here("data", "final", "expreset_pathvals_gtexV8.rds"))

save_feather <- function(x, path) {
  df <- data.frame(index = row.names(x), x)

  feather::write_feather(df, path)
}


save_feather(
  exp_data,
  file.path(here("data", "final", "expreset_Hinorm_gtexV8.rds.feather"))
)


save_feather(
  path_vals,
  file.path(here("data", "final", "expreset_pathvals_gtexV8.rds.feather"))
)
