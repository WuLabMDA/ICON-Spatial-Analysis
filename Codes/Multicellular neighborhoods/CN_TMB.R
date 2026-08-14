# delete work space
rm(list = ls(all = TRUE))
graphics.off()

library(ggplot2)
library(dplyr)
library(openxlsx)
library(pheatmap)
library(tidyverse)
library(tidytext)

load("Y:/Projects/ICON IMC/Final/S3 cellular neighborhood/Data/CN_workspace.rds")

unq <- rownames(densityTable)

tempMetadata <- readWorkbook(
  "Y:/Projects/ICON IMC/Experimental results/Latest section 1/Data/roiMetadata.xlsx",
  sheet = "Sheet 1",
  colNames = TRUE)

roiMetadata <- tempMetadata
#roiMetadata <- roiMetadata[roiMetadata$NSCLC=="YES",]
rownames(roiMetadata) <- roiMetadata$ROI.NO

for (i in 1:length(unq)){
  densityTable[unq[i], "Histology"] <- roiMetadata[unq[i],"LUAD_LUSC"]
  freqTable[unq[i], "Histology"] <- roiMetadata[unq[i],"LUAD_LUSC"]
  propTable[unq[i], "Histology"] <- roiMetadata[unq[i],"LUAD_LUSC"]
}

density_LUAD <- densityTable[densityTable$Histology == "LUAD",c(1:8)]

density_LUSC <- densityTable[densityTable$Histology == "LUSC",c(1:8)]

wesdata <- readWorkbook(
  "Y:/Projects/ICON IMC/Experimental results/WES section/Data/ICON_WES_Muhammad/SUMMARY_WES_for_Muhammad.xlsx",
  sheet = "Purity_Driver_mutations_TMB_108",
  colNames = TRUE)
rownames(wesdata) <- wesdata$TID

# load data
metaData <- readWorkbook(
  "Y:/Projects/ICON IMC/Clinical data/sortedClinical_LindaRevision.xlsx",
  sheet = "Bubble metaData",
  colNames = TRUE)
rownames(metaData) <- ptID <- metaData$SampleID


roiMetadata <- roiMetadata[roiMetadata$NSCLC=="YES",]
roiMetadata <- roiMetadata[roiMetadata$Location=="Tumor",]
# roiMetadata <- roiMetadata[complete.cases(roiMetadata),]
# rownames(roiMetadata) <- roiID <- roiMetadata$ROI.NO

roiMetadata_LUAD <- roiMetadata[roiMetadata$LUAD_LUSC=="LUAD",]
roiMetadata_LUSC <- roiMetadata[roiMetadata$LUAD_LUSC=="LUSC",]

# LUAD
dens_LUAD <- as.data.frame(matrix(0, nrow = length(ptID), ncol = ncol(density_LUAD)))
colnames(dens_LUAD) <- colnames(density_LUAD)
rownames(dens_LUAD) <- ptID

for (ii in 1:length(ptID)){
  tidx_LUAD <- roiMetadata_LUAD[roiMetadata_LUAD$SampleID==ptID[ii],]
  numROI <- nrow(tidx_LUAD)
  
  tdata_LUAD <- density_LUAD[rownames(density_LUAD) %in% rownames(tidx_LUAD),]
  dens_LUAD[ptID[ii],colnames(tdata_LUAD)] <- colSums(tdata_LUAD)/numROI
}
dens_LUAD <- dens_LUAD[complete.cases(dens_LUAD),]
dens_LUAD$TID <- metaData[rownames(dens_LUAD),"TID"]
rownames(dens_LUAD) <- dens_LUAD$TID

dens_LUAD$TMB <- wesdata[rownames(dens_LUAD),"TMB"]
dens_LUAD$TNB <- wesdata[rownames(dens_LUAD),"Neoantig.Burden"]
dens_LUAD$cnvGain <- wesdata[rownames(dens_LUAD),"CNV.gain"]
dens_LUAD$cnvLoss <- wesdata[rownames(dens_LUAD),"CNV.loss"]
dens_LUAD <- dens_LUAD[complete.cases(dens_LUAD),]

corr_LUAD <- as.data.frame(matrix(0, nrow = 8, ncol = 4))
colnames(corr_LUAD) <- c("TMB","NeoAntigen burden","CNV gain","CNV loss")
rownames(corr_LUAD) <- colnames(dens_LUAD)[1:8]

cormat_LUAD <- round(cor(dens_LUAD[, !colnames(dens_LUAD) %in% c("TID")]),2)
cormat_LUAD <- cormat_LUAD[c(9:12),c(1:8)]

pheatmap(cormat_LUAD, 
         color = colorRampPalette(c("dark blue", "white", "dark red"))(100), 
         scale = "none")



corr_df_luad <- as.data.frame(cormat_LUAD) %>%
  rownames_to_column("GenomicFeature") %>%
  pivot_longer(
    cols = -GenomicFeature,
    names_to = "CN",
    values_to = "Correlation"
  )

ggplot(corr_df_luad, aes(CN, GenomicFeature)) +
  
  geom_point(
    aes(size = abs(Correlation),
        fill = Correlation),
    shape = 21,
    color = "black",
    stroke = 0.4
  ) +
  
  scale_fill_gradient2(
    low = "#2C7BB6",
    mid = "white",
    high = "#D7191C",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Correlation"
  ) +
  
  scale_size(range = c(3, 12), name = "|Correlation|") +
  
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.position = "right"
  ) +
  
  labs(x = "", y = "")



cnColors <- c("#17BECF",  
              "#33A02C",  
              "#FF7F00",  
              "#e0ab09",  
              "#B15928",
              "#8d67f5",
              "#E31A1C",
              "#1F78B4")

cncelltypes <- setNames(cnColors, c("CN 1", "CN 2", "CN 3", "CN 4", "CN 5", "CN 6", "CN 7", "CN 8"))


corr_df_luad <- corr_df_luad %>%
  mutate(
    CN_reordered = reorder_within(CN, Correlation, GenomicFeature)
  )



ggplot(corr_df_luad,
       aes(Correlation, CN_reordered, color = CN)) +
  
  geom_point(size = 4) +
  
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "grey40") +
  
  facet_grid(GenomicFeature ~ ., scales = "free_y") +
  
  scale_y_reordered() +
  
  scale_color_manual(values = cncelltypes) +
  
  theme_classic(base_size = 14) +
  theme(
    axis.text.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    strip.text  = element_text(face = "bold", size = 14),
    legend.position = "right"
  ) +
  
  labs(title = "LUAD", x = "Correlation", y = "", color = "CN")





# LUSC
dens_LUSC <- as.data.frame(matrix(0, nrow = length(ptID), ncol = ncol(density_LUSC)))
colnames(dens_LUSC) <- colnames(density_LUSC)
rownames(dens_LUSC) <- ptID

for (ii in 1:length(ptID)){
  tidx_LUSC <- roiMetadata_LUSC[roiMetadata_LUSC$SampleID==ptID[ii],]
  numROI <- nrow(tidx_LUSC)
  
  tdata_LUSC <- density_LUSC[rownames(density_LUSC) %in% rownames(tidx_LUSC),]
  dens_LUSC[ptID[ii],colnames(tdata_LUSC)] <- colSums(tdata_LUSC)/numROI
}
dens_LUSC <- dens_LUSC[complete.cases(dens_LUSC),]
dens_LUSC$TID <- metaData[rownames(dens_LUSC),"TID"]
rownames(dens_LUSC) <- dens_LUSC$TID

dens_LUSC$TMB <- wesdata[rownames(dens_LUSC),"TMB"]
dens_LUSC$TNB <- wesdata[rownames(dens_LUSC),"Neoantig.Burden"]
dens_LUSC$cnvGain <- wesdata[rownames(dens_LUSC),"CNV.gain"]
dens_LUSC$cnvLoss <- wesdata[rownames(dens_LUSC),"CNV.loss"]
dens_LUSC <- dens_LUSC[complete.cases(dens_LUSC),]

corr_LUSC <- as.data.frame(matrix(0, nrow = 8, ncol = 4))
colnames(corr_LUSC) <- c("TMB","NeoAntigen burden","CNV gain","CNV loss")
rownames(corr_LUSC) <- colnames(dens_LUSC)[1:8]

cormat_LUSC <- round(cor(dens_LUSC[, !colnames(dens_LUSC) %in% c("TID")]),2)
cormat_LUSC <- cormat_LUSC[c(9:12),c(1:8)]

pheatmap(cormat_LUSC, 
         color = colorRampPalette(c("dark blue", "white", "dark red"))(100), 
         scale = "none")


corr_df_lusc <- as.data.frame(cormat_LUSC) %>%
  rownames_to_column("GenomicFeature") %>%
  pivot_longer(
    cols = -GenomicFeature,
    names_to = "CN",
    values_to = "Correlation"
  )

corr_df_lusc <- corr_df_lusc %>%
  mutate(
    CN_reordered = reorder_within(CN, Correlation, GenomicFeature)
  )


ggplot(corr_df_lusc,
       aes(Correlation, CN_reordered, color = CN)) +
  
  geom_point(size = 4) +
  
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "grey40") +
  
  facet_grid(GenomicFeature ~ ., scales = "free_y") +
  
  scale_y_reordered() +
  
  scale_color_manual(values = cncelltypes) +
  
  theme_classic(base_size = 14) +
  theme(
    axis.text.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    strip.text  = element_text(face = "bold", size = 14),
    legend.position = "right"
  ) +
  
  labs(title = "LUSC", x = "Correlation", y = "", color = "CN")


