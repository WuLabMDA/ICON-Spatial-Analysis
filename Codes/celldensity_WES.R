# delete work space
rm(list = ls(all = TRUE))
graphics.off()

library(openxlsx)
library(ggplot2)
library(RColorBrewer)

tempdata <- readRDS("Y:/Projects/ICON IMC/Processed data (rescale)/speCelltypes2.rds")
tempdata <- tempdata[, tempdata$celltypes != "Undefined"]

ptID <- unique(tempdata@colData@listData[["patient_id"]])

# call either normal or tumor
roiLocation <- "Tumor"
data <- tempdata[,colData(tempdata)$location==roiLocation]


# load data
mGenes <- readWorkbook(
  "Y:/Projects/ICON IMC/Experimental results/WES section/coMut plot/mutatedGenes.xlsx",
  sheet = "Sheet 1",
  colNames = TRUE,
  rowNames = TRUE)
mGenes <- mGenes[,colSums(mGenes !=0)>0]


tempMetadata <- readWorkbook(
  "Y:/Projects/ICON IMC/Experimental results/Latest section 1/Data/roiMetaData_LindaRevision.xlsx",
  sheet = "Sheet 1",
  colNames = TRUE)

roiMetadata <- tempMetadata
roiMetadata <- roiMetadata[roiMetadata$NSCLC=="YES",]
roiMetadata <- roiMetadata[roiMetadata$Location=="Tumor",]
rownames(roiMetadata) <- roiMetadata$ROI.NO


# unq1 <- unique(roiMetadata$SampleID)

celltypes <- as.character(unique(data@colData@listData[["celltypes"]]))

freqTable <- as.data.frame(matrix(0, nrow = nrow(roiMetadata), ncol = length(celltypes)))
colnames(freqTable) <- celltypes
rownames(freqTable) <- roiID <- roiMetadata$ROI.NO

propTable <- as.data.frame(matrix(0, nrow = nrow(roiMetadata), ncol = length(celltypes)))
colnames(propTable) <- celltypes
rownames(propTable) <- roiID 

densityTable <- as.data.frame(matrix(0, nrow = nrow(roiMetadata), ncol = length(celltypes)))
colnames(densityTable) <- celltypes
rownames(densityTable) <- roiID 

for (ii in 1:length(roiID)){
  message("Processing: ",ii," of ",length(roiID))
  tidx <- roiMetadata[roiMetadata$ROI.NO==roiID[ii],]
  fdata <- matrix(0, nrow = 1, ncol = length(celltypes))
  colnames(fdata) <- celltypes
  tdata <- data[,grepl(tidx$ROI.NO, data@colData@listData[["sample_id"]])]
  df <- as.numeric()
  
  # Cell density
  for (k in 1:length(celltypes)){
    df[k] <- sum(tdata@colData@listData[["celltypes"]]==celltypes[k])
  }
  fdata <- fdata + ((df*1e6)/tidx$Area)
  densityTable[roiID[ii],celltypes] <- fdata
  
  # Cell frequency
  cellFreqs <- as.data.frame(table(tdata@colData@listData[["celltypes"]]))
  freqTable[roiID[ii],as.character(cellFreqs$Var1)] <- cellFreqs$Freq
  
  # Cell proportiion
  cellProp <- as.data.frame(prop.table(table(tdata@colData@listData[["celltypes"]])))
  cellProp$Freq <- round(cellProp$Freq, digits = 2)
  propTable[roiID[ii],as.character(cellProp$Var1)] <- cellProp$Freq
}


int1 <- intersect(rownames(densityTable),rownames(roiMetadata))
# densityTable <- densityTable[int1,]
densityTable[int1,"TID"] <- roiMetadata[int1,"TID"]

int2 <- intersect(rownames(mGenes),densityTable$TID)

fdensityTable <- densityTable[densityTable$TID %in% int2,]

prop <- as.data.frame(matrix(0, nrow = length(int2), ncol = ncol(fdensityTable)-1))
colnames(prop) <- colnames(fdensityTable)[1:ncol(fdensityTable)-1]
rownames(prop) <- int2 

for (ii in 1:length(int2)){
  tidx <- roiMetadata[roiMetadata$TID==int2[ii],]
  numROI <- nrow(tidx)
  
  tdata <- fdensityTable[fdensityTable$TID %in% tidx$TID,]
  prop[int2[ii],colnames(tdata[,1:ncol(fdensityTable)-1])] <- colSums(tdata[,1:ncol(fdensityTable)-1])/numROI
}
prop <- prop[,colSums(prop !=0)>0]

fGenes <- mGenes[int2,]

cdata <- cbind(prop,fGenes)
cdata[,15:ncol(cdata)] <- lapply(cdata[,15:ncol(cdata)], factor)
features <- colnames(cdata)

corr <- pvalues <- as.data.frame(matrix(0, nrow = ncol(prop), ncol = ncol(fGenes)))
rownames(corr) <- rownames(pvalues) <- colnames(prop)
colnames(corr) <- colnames(pvalues) <- colnames(fGenes)

for (i in 1:ncol(fGenes)){
  for (j in 1:ncol(prop)){
    model <- lm(cdata[[j]] ~ cdata[[i+14]], data = cdata)
    corr[j,i] <- summary(model)$r.squared
    pvalues[j,i] <- summary(model)$coefficients[,4][2]
  }
}

corr$celltypes <- rownames(corr)
pvalues$celltypes <- rownames(pvalues)

corr <- gather(corr, genes, Rsquared, -celltypes) 
pvalues <- gather(pvalues, genes, Pvalues, -celltypes)

fdata <- corr
fdata$Pvalues <- pvalues$Pvalues
fdata[is.na(fdata$Pvalues),"Pvalues"] <- 1


colGEX = c("grey85", brewer.pal(7, "Reds"))
ggplot(fdata, aes(genes, celltypes, size = Pvalues, color = Rsquared)) + 
  geom_point() + ggtitle("Major celltypes") +
  scale_size_continuous(limits = c(min(fdata$Pvalues), max(fdata$Pvalues)), range = c(10,1), breaks = c(0.5,0.05,0.005)) + 
  theme_linedraw(base_size = 18) + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) + 
  scale_color_gradientn(colors = colGEX, limits = c(0,0.7), na.value = colGEX[8])




selected_genes <- c(
  "KRAS",
  "EGFR",
  "TP53",
  "STK11",
  "KEAP1",
  "SMARCA4",
  "RB1",
  "PIK3CA",
  "CDKN2A",
  "BRAF",
  "MET",
  "NFE2L2",
  "ARID1A",
  "PDGFRA",
  "GNAS",
  "MAP2K4",
  "MAP3K1",
  "PTEN",
  "NOTCH1",
  "RET",
  "FGFR1",
  "FBXW7"
)

fdata_sub <- fdata %>%
  filter(genes %in% selected_genes) %>%
  mutate(
    genes = factor(genes, levels = selected_genes),
    logP = -log10(Pvalues)
  )



colGEX <- c("grey90", brewer.pal(7, "Reds"))

ggplot(fdata_sub, aes(genes, celltypes)) +
  
  # Non-significant (no border)
  geom_point(
    data = subset(fdata_sub, Pvalues >= 0.05),
    aes(size = logP, fill = Rsquared),
    shape = 21,
    color = "lightgrey",
    alpha = 1
  ) +
  
  # Significant (black border)
  geom_point(
    data = subset(fdata_sub, Pvalues < 0.05),
    aes(size = logP, fill = Rsquared),
    shape = 21,
    color = "black",
    stroke = 1
  ) +
  
  scale_size_continuous(
    range = c(0, 9),
    name = expression(-log[10](italic(p)))
  ) +
  
  scale_fill_gradientn(
    colors = colGEX,
    limits = c(0, 0.7),
    name = expression(R^2)
  ) +
  
  labs(
    title = "Major cell types",
    x = "Genes",
    y = "Cell types"
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    panel.grid.major = element_line(color = "grey90", size = 0.3),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )


