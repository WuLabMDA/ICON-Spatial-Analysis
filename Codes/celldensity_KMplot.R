rm(list = ls(all = TRUE))
graphics.off()

library(ggplot2)
library(openxlsx)
library(survival)
library(survminer)
library(SummarizedExperiment)


tempdata <- readRDS("Y:/Projects/ICON IMC/Processed data (rescale)/speCelltypes2.rds")

ptID <- unique(tempdata@colData@listData[["patient_id"]])

# call either normal or tumor
roiLocation <- "Normal"
data <- tempdata[,colData(tempdata)$location==roiLocation]
# data <- tempdata

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

celltypes <- as.character(unique(data@colData@listData[["celltypes"]]))

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
      df[k] <- sum(tdata@colData@listData[["celltypes"]]==celltypes[k])
    }
    fdata <- fdata + ((df*1e6)/tidx$Area[j])
  }
  densityTable[unq1[ii],celltypes] <- fdata/numROI
}

int <- intersect(rownames(densityTable),rownames(metaData))
metaData <- metaData[int,]
metaData$endo <- densityTable[int,c("Endothelial cells")]
metaData$macro <- densityTable[int,c("Macrophages")]
metaData$cd4 <- densityTable[int,c("CD4 T cells")]
metaData$cd8 <- densityTable[int,c("CD8 T cells")]
metaData$treg <- densityTable[int,c("Treg")]
metaData$bcells <- densityTable[int,c("B cells")]
metaData$epi <- densityTable[int,c("Epithelial cells")]
metaData$fibro <- densityTable[int,c("Fibroblast")]
metaData$other <- densityTable[int,c("Other immune")]
metaData$mdsc <- densityTable[int,c("MDSC")]
metaData$mono <- densityTable[int,c("Monocytes")]
metaData$dendri <- densityTable[int,c("Dendritic")]
metaData$neutro <- densityTable[int,c("Neutrophils")]
metaData$nk <- densityTable[int,c("NK cells")]


metaData$Endothelial[densityTable$`Endothelial cells`<= mean(densityTable$`Endothelial cells`)] <- "Low density"
metaData$Endothelial[densityTable$`Endothelial cells`> mean(densityTable$`Endothelial cells`)] <- "High density"

metaData$Macrophages[densityTable$Macrophages <= mean(densityTable$Macrophages)] <- "Low density"
metaData$Macrophages[densityTable$Macrophages> mean(densityTable$Macrophages)] <- "High density"

metaData$CD4Tcells[densityTable$`CD4 T cells` <= mean(densityTable$`CD4 T cells`)] <- "Low density"
metaData$CD4Tcells[densityTable$`CD4 T cells`> mean(densityTable$`CD4 T cells`)] <- "High density"

metaData$CD8Tcells[densityTable$`CD8 T cells` <= mean(densityTable$`CD8 T cells`)] <- "Low density"
metaData$CD8Tcells[densityTable$`CD8 T cells`> mean(densityTable$`CD8 T cells`)] <- "High density"

metaData$Treg[densityTable$Treg <= mean(densityTable$Treg)] <- "Low density"
metaData$Treg[densityTable$Treg> mean(densityTable$Treg)] <- "High density"

metaData$Bcell[densityTable$`B cells` <= mean(densityTable$`B cells`)] <- "Low density"
metaData$Bcell[densityTable$`B cells`> mean(densityTable$`B cells`)] <- "High density"

metaData$Epithelial[densityTable$`Epithelial cells` <= mean(densityTable$`Epithelial cells`)] <- "Low density"
metaData$Epithelial[densityTable$`Epithelial cells`> mean(densityTable$`Epithelial cells`)] <- "High density"

metaData$Fibroblast[densityTable$Fibroblast <= mean(densityTable$Fibroblast)] <- "Low density"
metaData$Fibroblast[densityTable$Fibroblast> mean(densityTable$Fibroblast)] <- "High density"

metaData$OtherImmune[densityTable$`Other immune` <= mean(densityTable$`Other immune`)] <- "Low density"
metaData$OtherImmune[densityTable$`Other immune`> mean(densityTable$`Other immune`)] <- "High density"

metaData$MDSC[densityTable$MDSC <= mean(densityTable$MDSC)] <- "Low density"
metaData$MDSC[densityTable$MDSC> mean(densityTable$MDSC)] <- "High density"

metaData$Monocytes[densityTable$Monocytes <= mean(densityTable$Monocytes)] <- "Low density"
metaData$Monocytes[densityTable$Monocytes> mean(densityTable$Monocytes)] <- "High density"

metaData$Dendritic[densityTable$Dendritic <= mean(densityTable$Dendritic)] <- "Low density"
metaData$Dendritic[densityTable$Dendritic> mean(densityTable$Dendritic)] <- "High density"

metaData$Neutrophils[densityTable$Neutrophils <= mean(densityTable$Neutrophils)] <- "Low density"
metaData$Neutrophils[densityTable$Neutrophils> mean(densityTable$Neutrophils)] <- "High density"

metaData$NKcells[densityTable$`NK cells` <= mean(densityTable$`NK cells`)] <- "Low density"
metaData$NKcells[densityTable$`NK cells`> mean(densityTable$`NK cells`)] <- "High density"


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
               ~ CD4Tcells, data = metaData)

summary(coxph(Surv(OS_month/12,OS_status)  ~ CD4Tcells, data = metaData))

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
                ~ CD4Tcells, data = metaData)

summary(coxph(Surv(RFS_month/12,RFS_status)  ~ CD4Tcells, data = metaData))

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

