# delete work space
rm(list = ls(all = TRUE))
graphics.off()

library(ggplot2)
library(cowplot)
library(patchwork)
library(dplyr)
library(openxlsx)
library(Seurat)
library(viridis)
library(ggpubr)
library(reshape2)

tempdata <- readRDS("Y:/Projects/ICON IMC/Processed data (rescale)/speCelltypes2.rds")

tempdata <- tempdata[, tempdata$celltypes != "Undefined"]

data <- tempdata

# load data
metaData <- readWorkbook(
  "Y:/Projects/ICON IMC/Clinical data/sortedClinical_012224.xlsx",
  sheet = "Bubble metaData",
  colNames = TRUE)
rownames(metaData) <- metaData$SampleID

# load data
tempMetadata <- readWorkbook(
  "Y:/Projects/ICON IMC/Clinical data/roiMetadata2.xlsx",
  sheet = "Sheet 1",
  colNames = TRUE)
tempMetadata <- tempMetadata[,-1]

roiMetadata <- tempMetadata
rownames(roiMetadata) <- roiMetadata$ROI.NO

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

densityTable$Location <- roiMetadata[roiID,"Location"]
freqTable$Location <- roiMetadata[roiID,"Location"]
propTable$Location <- roiMetadata[roiID,"Location"]


################
# load data
tempTCR <- readWorkbook(
  "Y:/Projects/ICON IMC/Experimental results/TCR section/Data/20240319_ICON_all.xlsx",
  sheet = "Normal_TCR",
  colNames = TRUE)

TCR <- tempTCR[,c(1,6,7)]
colnames(TCR) <- c("SampleID","Richness","Clonality")

roi_location <- "Normal"
# roiMetadata <- roiMetadata[roiMetadata$NSCLC=="YES",]
roiMetadata <- roiMetadata[roiMetadata$Location==roi_location,]
rownames(roiMetadata) <- roiID <- roiMetadata$ROI.NO

densityTable <- densityTable[densityTable$Location == roi_location,c(1:14)]

freqTable <- freqTable[freqTable$Location == roi_location,c(1:14)]

propTable <- propTable[propTable$Location == roi_location,c(1:14)]


pTID <- roiMetadata[,c(1,2,3)]

int <- intersect(TCR$SampleID,pTID$TID)

fTCR <- TCR[TCR$SampleID %in% int,]

prop <- as.data.frame(matrix(0, nrow = length(int), ncol = ncol(propTable)))
colnames(prop) <- colnames(propTable)
rownames(prop) <- int 

ffTCR <- as.data.frame(matrix(0, nrow = length(int), ncol = ncol(fTCR)-1))
colnames(ffTCR) <- colnames(fTCR)[c(2:3)]
rownames(ffTCR) <- int

for (ii in 1:length(int)){
  tidx <- roiMetadata[roiMetadata$TID==int[ii],]
  numROI <- nrow(tidx)
  
  tdata <- propTable[rownames(propTable) %in% rownames(tidx),]
  prop[int[ii],colnames(tdata)] <- colSums(tdata)/numROI
  
  tid <- fTCR[fTCR$SampleID==int[ii],c(2:3)]
  tid$Richness <- as.numeric(tid$Richness)
  tid$Clonality <- as.numeric(tid$Clonality)
  numS <- nrow(tid)
  ffTCR[int[ii],] <- colSums(tid)/numS
}

# prop <- prop[complete.cases(prop),]
# int1 <- intersect(rownames(prop),rownames(metaData))
fprop <- cbind(prop[int,],ffTCR[int,])

# library(car)
# fprop[,c(1:16)] <- 1*logit(fprop[,c(1:16)], percents = FALSE, adjust = 0.25)


cor_test <- cor.test(
  fprop$`CD8 T cells`,
  fprop$Clonality,
  method = "spearman"
)

rho <- round(cor_test$estimate, 2)
pval <- signif(cor_test$p.value, 2)


ggplot(
    fprop,
    aes(
      x = `CD8 T cells`,
      y = Clonality
    )
  ) +
  
geom_smooth(
  method = "lm",
  colour = "#D73027",
  fill = "#D9D9D9",
  linewidth = 1.2
) +

geom_point(
  shape = 21,
  size = 3.6,
  stroke = 0.6,
  fill = "#2C7FB8",
  colour = "black"
) +

annotate(
  "text",
  x = min(fprop$`CD8 T cells`),
  y = max(fprop$Clonality),
  hjust = 0,
  vjust = 1,
  size = 5,
  fontface = "bold",
  label = paste0(
    "Spearman ",
    "\u03C1 = ",
    rho,
    "\nP = ",
    format.pval(pval, digits = 2)
  )
) +
  
labs(
  title = "CD8 T-cell abundance is associated with TCR clonality",
  x = "CD8 T-cell proportion",
  y = "TCR clonality"
) +
  
theme_classic(base_size = 16) +
  
  theme(
    
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 18
    ),
    
    axis.title = element_text(
      face = "bold",
      size = 15
    ),
    
    axis.text = element_text(
      face = "bold",
      colour = "black",
      size = 13
    )
    
  )



cor_test <- cor.test(
  fprop$`CD8 T cells`,
  fprop$Richness,
  method = "spearman"
)

rho <- round(cor_test$estimate, 2)
pval <- signif(cor_test$p.value, 2)

ggplot(
    fprop,
    aes(
      x = `CD8 T cells`,
      y = Richness
    )
  ) +

geom_smooth(
  method = "lm",
  colour = "#D73027",
  fill = "#D9D9D9",
  linewidth = 1.2
) +

geom_point(
  shape = 21,
  size = 3.6,
  stroke = 0.6,
  fill = "#2C7FB8",
  colour = "black"
) +

annotate(
  "text",
  x = min(fprop$`CD8 T cells`),
  y = max(fprop$Richness),
  hjust = 0,
  vjust = 1,
  size = 5,
  fontface = "bold",
  label = paste0(
    "Spearman ",
    "\u03C1 = ",
    rho,
    "\nP = ",
    format.pval(pval, digits = 2)
  )
) +
  
labs(
  title = "CD8 T-cell abundance is associated with TCR richness",
  x = "CD8 T-cell proportion",
  y = "TCR richness"
) +
  
theme_classic(base_size = 16) +
  
  theme(
    
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 18
    ),
    
    axis.title = element_text(
      face = "bold",
      size = 15
    ),
    
    axis.text = element_text(
      face = "bold",
      colour = "black",
      size = 13
    )
    
  )



CN_long <- propTable %>%
  tibble::rownames_to_column("ROI") %>%
  pivot_longer(cols = -ROI,
               names_to = "Celltypes",
               values_to = "Proportion")

CN_mean <- colMeans(propTable, na.rm = TRUE)

CN_mean

barplot(sort(CN_mean, decreasing = TRUE),
        col = "steelblue",
        las = 2,
        ylab = "Mean Celltypes proportion")

library(dplyr)

CNs <- colnames(propTable)

res <- lapply(CNs, function(cn){
  
  cor_clon <- cor.test(fprop[[cn]], fprop$Clonality, method="spearman")
  cor_rich <- cor.test(fprop[[cn]], fprop$Richness, method="spearman")
  
  data.frame(
    CN = cn,
    clonality_rho = cor_clon$estimate,
    clonality_p = cor_clon$p.value,
    richness_rho = cor_rich$estimate,
    richness_p = cor_rich$p.value
  )
})

CN_cor <- bind_rows(res)

CN_cor

fprop <- fprop[complete.cases(fprop), ]
fprop$CD8_group <- ifelse(fprop$`CD8 T cells` > median(fprop$`CD8 T cells`),
                          "High","Low")

fprop$CD4_group <- ifelse(fprop$`CD4 T cells` > median(fprop$`CD4 T cells`),
                          "High","Low")


fprop$CD8_group <- factor(fprop$CD8_group, levels = c("Low","High"))
fprop$CD4_group <- factor(fprop$CD4_group, levels = c("Low","High"))


fprop$CD8_group <- factor(
  fprop$CD8_group,
  levels = c("Low", "High")
)

pval <- wilcox.test(
  Clonality ~ CD8_group,
  data = fprop
)$p.value

label <- paste0(
  "Wilcoxon P = ",
  format.pval(
    pval,
    digits = 2,
    eps = 1e-4
  )
)

ggplot(
    
    fprop,
    
    aes(
      x = CD8_group,
      y = Clonality,
      fill = CD8_group
    )
    
  ) +

geom_violin(
  
  trim = FALSE,
  
  alpha = 0.30,
  
  linewidth = 0.4,
  
  colour = NA
  
) +

geom_boxplot(
  
  width = 0.18,
  
  outlier.shape = NA,
  
  linewidth = 0.5,
  
  colour = "black"
  
) +

geom_jitter(
  
  aes(colour = CD8_group),
  
  width = 0.08,
  
  size = 2.3,
  
  alpha = 0.8,
  
  show.legend = FALSE
  
) +

stat_summary(
  
  fun = median,
  
  geom = "point",
  
  shape = 23,
  
  size = 3,
  
  fill = "black",
  
  colour = "black"
  
) +
  
annotate(
  
  "text",
  
  x = 1.5,
  
  y = max(fprop$Clonality) * 1.08,
  
  label = label,
  
  fontface = "bold",
  
  size = 5
  
) +
  
scale_fill_manual(
  
  values = c(
    
    Low = "#56B4E9",
    
    High = "#D55E00"
    
  )
  
) +
  
  scale_colour_manual(
    
    values = c(
      
      Low = "#2C7FB8",
      
      High = "#B2182B"
      
    )
    
  ) +
  
labs(
  
  title = "High CD8 T-cell abundance is associated with increased TCR clonality",
  
  x = "CD8 T-cell abundance",
  
  y = "TCR clonality"
  
) +
  
theme_classic(base_size = 16) +
  
  theme(
    
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 18
    ),
    
    axis.title = element_text(
      face = "bold",
      size = 15
    ),
    
    axis.text = element_text(
      face = "bold",
      colour = "black",
      size = 13
    ),
    
    legend.position = "none"
    
  )


