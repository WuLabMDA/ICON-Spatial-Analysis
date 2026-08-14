# delete work space
rm(list = ls(all = TRUE))
graphics.off()

library(SpatialExperiment)
library(imcRtools)
library(openxlsx)
library(pheatmap)

tempdata <- readRDS("Y:/Projects/ICON IMC/Processed data (rescale)/speCelltypes2.rds")

data <- tempdata

## cellular neighborhood analysis 
# construct interaction graph
data <- buildSpatialGraph(data, img_id = "sample_id", type = "knn", k = 20)

data <- aggregateNeighbors(data, 
                           colPairName = "knn_interaction_graph", 
                           aggregate_by = "metadata", 
                           count_by = "celltypes")

cn_1 <- kmeans(data$aggregatedNeighbors, centers = 8)
data$cn_celltypes <- as.factor(paste(rep("CN",length(cn_1$cluster)),cn_1$cluster))

# fraction of each CN made up of each cell type
for_plot <- prop.table(table(as.character(data$cn_celltypes), data$celltypes), margin = 1)

pheatmap(for_plot, 
         color = colorRampPalette(c("dark blue", "white", "dark red"))(100), 
         scale = "column", cluster_rows = FALSE)

norm_data <- data[,data@colData@listData[["location"]] == "Normal"]

nfor_plot <- prop.table(table(as.character(norm_data$cn_celltypes), norm_data$celltypes), margin = 1)

pheatmap(nfor_plot, 
         color = colorRampPalette(c("dark blue", "white", "dark red"))(100), 
         scale = "column", cluster_rows = FALSE)


tumor_data <- data[,data@colData@listData[["location"]] == "Tumor"]

tfor_plot <- prop.table(table(as.character(tumor_data$cn_celltypes), tumor_data$celltypes), margin = 1)

pheatmap(tfor_plot, 
         color = colorRampPalette(c("dark blue", "white", "dark red"))(100), 
         scale = "column", cluster_rows = FALSE)


# load data
tempMetadata <- readWorkbook(
  "Y:/Projects/ICON IMC/Clinical data/roiMetadata2.xlsx",
  sheet = "Sheet 1",
  colNames = TRUE)
tempMetadata <- tempMetadata[,-1]


roiMetadata <- tempMetadata
rownames(roiMetadata) <- roiMetadata$ROI.NO

cn_celltypes <- as.character(unique(data@colData@listData[["cn_celltypes"]]))

freqTable <- as.data.frame(matrix(0, nrow = nrow(roiMetadata), ncol = length(cn_celltypes)))
colnames(freqTable) <- cn_celltypes
rownames(freqTable) <- roiID <- roiMetadata$ROI.NO

propTable <- as.data.frame(matrix(0, nrow = nrow(roiMetadata), ncol = length(cn_celltypes)))
colnames(propTable) <- cn_celltypes
rownames(propTable) <- roiID 

densityTable <- as.data.frame(matrix(0, nrow = nrow(roiMetadata), ncol = length(cn_celltypes)))
colnames(densityTable) <- cn_celltypes
rownames(densityTable) <- roiID 

for (ii in 1:length(roiID)){
  message("Processing: ",ii," of ",length(roiID))
  tidx <- roiMetadata[roiMetadata$ROI.NO==roiID[ii],]
  fdata <- matrix(0, nrow = 1, ncol = length(cn_celltypes))
  colnames(fdata) <- cn_celltypes
  tdata <- data[,grepl(tidx$ROI.NO, data@colData@listData[["sample_id"]])]
  df <- as.numeric()
  
  # Cell density
  for (k in 1:length(cn_celltypes)){
    df[k] <- sum(tdata@colData@listData[["cn_celltypes"]]==cn_celltypes[k])
  }
  fdata <- fdata + ((df*1e6)/tidx$Area)
  densityTable[roiID[ii],cn_celltypes] <- fdata
  
  # Cell frequency
  cellFreqs <- as.data.frame(table(tdata@colData@listData[["cn_celltypes"]]))
  freqTable[roiID[ii],as.character(cellFreqs$Var1)] <- cellFreqs$Freq
  
  # Cell proportiion
  cellProp <- as.data.frame(prop.table(table(tdata@colData@listData[["cn_celltypes"]])))
  cellProp$Freq <- round(cellProp$Freq, digits = 2)
  propTable[roiID[ii],as.character(cellProp$Var1)] <- cellProp$Freq
}


# save.image("Y:/Projects/ICON IMC/Final/S3 cellular neighborhood/Data/CN_workspace.rds")

# write.xlsx(densityTable, file = "Y:/Projects/ICON IMC/Experimental results/Latest section 1/Data/densityCN.xlsx",
#            colNames = TRUE, rowNames = TRUE)
# write.xlsx(freqTable, file = "Y:/Projects/ICON IMC/Experimental results/Latest section 1/Data/frequencyCN.xlsx",
#            colNames = TRUE, rowNames = TRUE)
# write.xlsx(propTable, file = "Y:/Projects/ICON IMC/Experimental results/Latest section 1/Data/proportionCN.xlsx",
#            colNames = TRUE, rowNames = TRUE)

