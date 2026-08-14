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

tempdata <- readRDS("Y:/Projects/ICON IMC/Experimental results/Section 2/Data/speCelltypes2.rds")

CD4T <- tempdata[,colData(tempdata)$celltypes=="CD4 T cells"]
CD4T <- CD4T@assays@data@listData[["exprs"]]

features <- c("GranzymeB","TIM3","ICOS","TIGIT","PD-L1","LAG3","PD-1","CTLA-4","Ki67","CD45RO")

CD4T <- CreateSeuratObject(counts = CD4T,
                           project = "ICON IMC",
                           assay = "CD4T")

CD4T <- NormalizeData(CD4T, verbose = FALSE)
CD4T <- FindVariableFeatures(CD4T, selection.method = "vst", nfeatures = 45,
                             verbose = FALSE)
CD4T <- ScaleData(CD4T, verbose = FALSE)
CD4T <- RunPCA(CD4T, npcs = 35, verbose = FALSE, features = rownames(CD4T))
CD4T <- RunUMAP(CD4T, reduction = "pca", dims = 1:35, verbose = FALSE)

# Cluster endothelail data
CD4T <- FindNeighbors(CD4T, dims = 1:35)
CD4T <- FindClusters(CD4T, resolution = 0.04, algorithm = 2)

DimPlot(CD4T, reduction = "umap", label = TRUE, repel = TRUE) 

DotPlot(CD4T, features=features, dot.scale=10, cols="RdBu")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

cellNames <- c("Prolif CD4 T cells","Prolif CD4 T cells","CD4 T cells","CD4 T cells","CD4 T cells","Exhausted CD4 T cells",
               "Memory CD4 T cells","CD4 T cells")

names(cellNames) <- levels(CD4T)
CD4T <- RenameIdents(object = CD4T, cellNames)
DimPlot(CD4T, label = FALSE, reduction = 'umap')
DotPlot(CD4T, features=features, dot.scale=10, cols="RdBu")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

FeaturePlot(CD4T, features = features,cols = c("gray","red"))

CD4T@meta.data[["celltypes"]] <- CD4T@active.ident

CD4clusters <- as.data.frame(as.character(CD4T$celltypes))
rownames(CD4clusters) <- colnames(CD4T)
colnames(CD4clusters) <- c("celltypes")

newClusters <- as.data.frame(as.character(colData(tempdata)$celltypes))
rownames(newClusters) <- colnames(tempdata)
colnames(newClusters) <- c("celltypes")

int <- intersect(rownames(newClusters),rownames(CD4clusters))

newClusters[int,c("celltypes")] <- CD4clusters[int,c("celltypes")]
all(colnames(tempdata)==rownames(newClusters))

colData(tempdata)$cellsubtypes <- as.factor(newClusters$celltypes)

# associate CD4T subtypes with recurrence/survival
# call either normal or tumor
roiLocation <- "Tumor"
data <- tempdata[,colData(tempdata)$location==roiLocation]
# data <- tempdata

ptID <- unique(tempdata@colData@listData[["patient_id"]])

# load data
metaData <- readWorkbook(
  "Y:/Projects/ICON IMC/Clinical data/sortedClinical_LindaRevision.xlsx",
  sheet = "OS_RFS",
  colNames = TRUE)
rownames(metaData) <- metaData$SampleID

# load data
tempMetadata <- readWorkbook(
  "Y:/Projects/ICON IMC/Clinical data/roiMetadata2.xlsx",
  sheet = "Sheet 1",
  colNames = TRUE)
tempMetadata <- tempMetadata[,-1]

roiMetadata <- tempMetadata[tempMetadata$Location==roiLocation,]
# roiMetadata <- tempMetadata
unq1 <- unique(roiMetadata$SampleID)

celltypes <- as.character(unique(data@colData@listData[["cellsubtypes"]]))

densityTable <- as.data.frame(matrix(0, nrow = length(ptID), ncol = length(celltypes)))
colnames(densityTable) <- celltypes
rownames(densityTable) <- ptID 
# unq <- unique(data@colData@listData[["sample_id"]])

for (ii in 1:length(unq1)){
  tidx <- roiMetadata[roiMetadata$SampleID==unq1[ii],]
  numROI <- nrow(tidx)
  fdata <- matrix(0, nrow = 1, ncol = length(celltypes))
  colnames(fdata) <- celltypes
  for (j in 1:numROI){
    tdata <- data[,grepl(tidx$ROI.NO[j], data@colData@listData[["sample_id"]])]
    df <- as.numeric()
    for (k in 1:length(celltypes)){
      df[k] <- sum(tdata@colData@listData[["cellsubtypes"]]==celltypes[k])
    }
    fdata <- fdata + ((df*1e6)/tidx$Area[j])
  }
  densityTable[unq1[ii],celltypes] <- fdata/numROI
}

int <- intersect(rownames(densityTable),rownames(metaData))
metaData <- metaData[int,]
metaData$prolifCD4 <- densityTable[int,c("Prolif CD4 T cells")]
metaData$memCD4 <- densityTable[int,c("Memory CD4 T cells")]
metaData$exhCD4 <- densityTable[int,c("Exhausted CD4 T cells")]

metaData$prolif[densityTable$`Prolif CD4 T cells`<= mean(densityTable$`Prolif CD4 T cells`)] <- "Low density"
metaData$prolif[densityTable$`Prolif CD4 T cells`> mean(densityTable$`Prolif CD4 T cells`)] <- "High density"

metaData$mem[densityTable$`Memory CD4 T cells`<= mean(densityTable$`Memory CD4 T cells`)] <- "Low density"
metaData$mem[densityTable$`Memory CD4 T cells`> mean(densityTable$`Memory CD4 T cells`)] <- "High density"

metaData$exh[densityTable$`Exhausted CD4 T cells`<= mean(densityTable$`Exhausted CD4 T cells`)] <- "Low density"
metaData$exh[densityTable$`Exhausted CD4 T cells`> mean(densityTable$`Exhausted CD4 T cells`)] <- "High density"

custom_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(size = 14, color = "black",hjust=0.5,face = "bold"),
      axis.text.x = element_text(size = 14, color = "black", face = "bold"),
      legend.text = element_text(size = 14, color = "black", face = "bold"),
      legend.title = element_text(size = 14, color = "black", face = "bold"),
      axis.text.y = element_text(size = 14, color = "black", face = "bold"),
      axis.title.x = element_text(size = 14, color = "black", face = "bold"),
      axis.title.y = element_text(size = 14, color = "black", face = "bold", angle = 90) , #angle=(90))
    )
}

# fit <- coxph(Surv(RFS_month,RFS_status)
#                      ~ endothelial, data = metaData)
fit <- survfit(Surv(OS_month/12,OS_status)
               ~ mem, data = metaData)

summary(coxph(Surv(OS_month/12,OS_status)  ~ mem, data = metaData))

ggsurvplot(fit, data = metaData,title = "",ggtheme=custom_theme(),
           conf.int = FALSE,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           xlab = "Time (Years)",
           ylab = "Overall Survival (%)",
           xlim = c(0, 8),
           risk.table.fontsize =5,
           size = 2,
           linetype = "solid",
           palette = c("#00468BFF","#be0000"),
           
           risk.table.col = "strata",
           #legend = "bottom",
           legend.title = "")

fit2 <- survfit(Surv(RFS_month/12,RFS_status)
                ~ mem, data = metaData)

summary(coxph(Surv(RFS_month/12,RFS_status)  ~ mem, data = metaData))

ggsurvplot(fit2, data = metaData,title = "",ggtheme=custom_theme(),
           conf.int = FALSE,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           xlab = "Time (Years)",
           ylab = "Recurrence Free Survival (%)",
           xlim = c(0, 8),
           risk.table.fontsize =5,
           size = 2,
           linetype = "solid",
           palette = c("#00468BFF","#be0000"),
           
           risk.table.col = "strata",
           #legend = "bottom",
           legend.title = "")
