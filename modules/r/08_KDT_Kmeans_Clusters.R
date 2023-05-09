############################################
## Project: RP DRexML
## Script purpose:  Cluster Analysis between the relevant KDTs + cluster behavior along hallmarks
## Date: 25.02.2023
## Author: Marina Esteban-Medina
#######################################

if (!require(pacman, quietly = TRUE)){
  install.packages(pacman)
}else{
  library(pacman)
}

p_load("here","cluster", "factoextra","NbClust","tidyr", "tibble", "NMF")


if(!dir.exists(here("results/tables"))){
  dir.create(here("results/tables"))
}

if(!dir.exists(here("results/figures"))){
  dir.create(here("results/figures"))
}

if(!dir.exists(here("rds"))){
  dir.create(here("rds"))
}

### 1. Read files for the analysis ####

## Load the relevance scores matrix with the threshold selection matrix for filtering ###
data_folder = here("results","ml")

shap <- read.delim(file = file.path(data_folder,"shap_summary_symbol.tsv"), header = T) %>% as.data.frame()
rownames(shap)<- shap$circuit_name
shap <- shap[ ,-1]

threshold <- read.delim(file = file.path(data_folder,"shap_selection_symbol.tsv"), header = T) %>% as.data.frame() 
rownames(threshold) <- threshold$circuit_name
threshold <- threshold[ ,-1]

## Load the filtered matrix with stable circuits and relevant KDTs normalized  values to a [-1,1] range
shap_relevant_stable <- read.delim(here("results/tables/shap_relevant_stable.tsv"))
shap_entrez_relevant_stable <- read.delim( here("results/tables/shap_entrez_relevant_stable.tsv"))

## Read the "rds" files for de KDT_by function and DRUG_byfunction ####
relevant_KDTbyfunc <- readRDS(file = here("rds", "relevant_KDTbyfunc.rds"))
relevant_KDTbyfunc_percent <- apply(relevant_KDTbyfunc, 2, function(x) round(x/max(x)*100)) %>% data.frame(.) ## Transform the count of hallmarks in percentage of the hallmark represented in a KDT or drug

## Selected from the raw SHAP matrix, relevant results
df <- shap[rownames(shap_relevant_stable), colnames(shap_relevant_stable)] %>% t(.)  %>% data.frame(.)
df <-  apply(df, 2, function(x) x/max(abs(x)))  %>% data.frame(.)

### 2. Estimate optimal nº of  clusters #####
fviz_nbclust(df, kmeans, method = "gap")
ggsave(plot = last_plot(), file = here("results", "figures", "nClusters_optimal_methodGapStat.png"))

mat <- t(apply(shap, 1, function(x) x/max(abs(x)))) ## Rescale the score values so they are in a scale -1,1 on the whole matrix
mat[which(threshold == 0, arr.ind = T)] <- 0  ## set to 0 the non relevant ones
mat <- mat[rownames(shap_relevant_stable), colnames(shap_relevant_stable)]
df <- t(mat) %>% data.frame(.)

## Calculate optimal nº of cluster again
resnumclust <-NbClust(df, distance = "euclidean", min.nc=2, max.nc= 8, method = "kmeans", index = "gap", set.seed(1))
n_cluster <- fviz_nbclust(resnumclust)

### 3. K-means method for cluster analysis ####
set.seed(42)
k <- kmeans(df, centers = n_cluster["Number_clusters"], nstart = 25)

clusters <- data.frame( CLUSTER = k$cluster) %>% rownames_to_column("KDT")%>% .[order(.$CLUSTER, decreasing = F),]
write.table(clusters, here("results", "tables", "clusters_n3_relevant_stableKDTS.tsv"), row.names = F, col.names = T, sep = "\t", quote = F)

## Plot the clusters with ellipses
png(here("results", "figures","clusters3_heatmap.png"), height = 5000, width = 8000, res = 500)
fviz_cluster(k, data = df, ellipse.type = "euclid",repel = TRUE,star.plot = TRUE,  ggtheme = theme_minimal() ) #ellipse.type= "t", "norm", "euclid"
dev.off()

## Add clusters to the results with hallmarks annotations table
df_relevant_KDTbyfunc <- relevant_KDTbyfunc_percent[-1,] %>% tibble::add_column("CLUSTER" = clusters$CLUSTER[match(rownames(relevant_KDTbyfunc_percent)[-1], clusters$KDT)])
df_relevant_KDTbyfunc$CLUSTER <- factor(df_relevant_KDTbyfunc$CLUSTER)

data_long <- gather(df_relevant_KDTbyfunc, hallmark, Percentage, Fatty.acid.and.lipid.metabolism:Sensory.and.stimuli.transduction, factor_key = T )
data_long$hallmark <- gsub("\\.", " ", data_long$hallmark)

png(here("results", "figures","clusters3_heatmap_countsperHallmark.png"), height = 5000, width = 8000, res = 500)
ggplot(data_long, aes(as.factor(x = hallmark), y = Percentage, group=CLUSTER, colour = CLUSTER)) + 
  stat_summary(fun = mean, geom="pointrange", size = 1)+
  stat_summary(geom="line")+
  theme_minimal()+
  # geom_point(aes(shape=CLUSTER))+
  theme(axis.text.x = element_text(size= 20, family = "Helvetica", angle = 30, hjust = 1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size= 20, family = "Helvetica"),
        axis.title.y =  element_text(size= 20, family = "Helvetica"))+
  ylim(0,100)

dev.off()

###########################################################################
### Plot now the heatmap with the order of the clusters on the column order 
############################################################

## Load the drug_eff DF from Drugbank DB
drugbank_effects_tar<- readRDS(here("rds", "drugbank_effects_tar.rds"))
drug_ef<- read.delim(here("results/tables/drugEffects_KDT_relevant_stable_simpl.tsv"))

## Prepare heatMap function with drug effect annotations using NMF R pack
annot_testAheatmap_clust <- function(matrix, annotations, title, anot_colors, colors_main = "-RdYlBu",Colv = F, Rowv = F, hclustfun = "complete" ) {    
  NMF::aheatmap(matrix, color = colors_main, border_color = "white", annCol= annotations["Drug_effect"], Rowv = Rowv, Colv =Colv, ## If we set Rowv and Colv to false it will cluster qwith hclust funct.
                hclustfun =hclustfun, cexCol = 4, cexRow = 10, main =  title, 
                annColors = anot_colors )
}


## Select colors for drug effects
annot.color.col <- list('Drug_effect'=c('green','lightblue', "purple", "black", "grey"))

### Heatmap with colors indicating the strength of the relevance in the clusters order ####

mat <- t(apply(shap, 1, function(x) x/max(abs(x)))) ## Rescale the score values so they are in a scale -1,1 on the whole matrix
mat[which(threshold == 0, arr.ind = T)] <- 0  ## set to 0 the non relevant ones
mat <- mat[rownames(shap_relevant_stable), colnames(shap_relevant_stable)]
dim(mat) ## 207 cir, 109 KDTs

## Extract clusters from hierarchical clustering
heatmap.obj <- annot_testAheatmap_clust(matrix = mat, annotations = drug_ef[drug_ef$symbol %in% colnames(shap_relevant_stable),], anot_colors = annot.color.col,
                                        title = "RELEVANCE SCORE OF RELEVANT KDTs x RP circuits \n Score of how relevant is each KDT for the activity of each circuit with sign")

clusters_hierarch <- data.frame(CLUSTER_HIER = cutree(as.hclust(heatmap.obj$Colv), k=3)) %>% rownames_to_column("KDT") %>% merge(.,clusters ) %>% .[order(.$CLUSTER, decreasing = T),]
colnames(clusters_hierarch) <- c("KDT", "CLUSTER_HIER_COMPLETE", "CLUSTER_Kmeans" )
write.table(clusters_hierarch, here("results", "tables", "clusters_n3_kmeansVShierarch_relevant_stableKDTS.tsv"), row.names = F, col.names = T, sep = "\t", quote = F)

  
## Plot the score HeatMap with the drug effect annotations on top
png(here("results","figures","heatmap_relevancesSHAP_reorder_n3clusters.png"), height = 8000, width = 8000, res = 500)
annot_testAheatmap_clust(matrix = mat[,clusters_hierarch$KDT], annotations = drug_ef[drug_ef$symbol %in% colnames(shap_relevant_stable),], anot_colors = annot.color.col, Colv = T, Rowv = F,
                   title = "RELEVANCE SCORE OF RELEVANT KDTs x RP circuits \n Score of how relevant is each KDT for the activity of each circuit with sign")

dev.off()



