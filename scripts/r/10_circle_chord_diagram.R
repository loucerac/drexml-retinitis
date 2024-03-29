############################################
## Project: RP DRexML
## Script purpose: Chord diagrams representing the shared items among de datasets.
## Date: 05.03.2023
## Author: Marina Esteban-Medina
#######################################

library(circlize)

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

### Load de data we want to analyze ###

ha_ha <- readRDS(file.path(rds_folder, "hallmarks_byhallmarks.rds"))

# Set the colors in proper order (like in the spiderplot) from palette Set1
grid.col <- c("#FFFF33", "#E41A1C", "#984EA3","#4DAF4A", "#FF7F00", "#999999", "#A65628", "#377EB8","#F781BF") 
grid.col <- setNames(grid.col, colnames(ha_ha))

# KDts shared among Hallmarks on a chord diagram
png(file.path(figures_folder,"hallmarks_KDTshared_chord_downwards.png"), height = 8000, width = 12000, res = 500) 
chordDiagram(as.matrix(ha_ha), transparency = 0.5,
             annotationTrack = "grid", 
             preAllocateTracks = 1, 
              grid.col = grid.col)
             # directional = 1,
              #direction.type = c("diffHeight", "arrows"), 
              #link.arr.type = "big.arrow")

circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(CELL_META$xcenter ,
              ylim[1] + cm_h(20),
              sector.name,
              facing = "downward",
              niceFacing = T,
              adj = c(0.65,0.1),
              cex = 1.8,
              col= grid.col[sector.name],
              font = 2)

}, bg.border = NA)
dev.off()
