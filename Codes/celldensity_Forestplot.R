# delete work space
rm(list = ls(all = TRUE))
graphics.off()

library(openxlsx)
library(survival)
library(dplyr)
library(forestplot)
library(grid)
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


celltype_vars <- c(
  "Endothelial cells",
  "Macrophages",
  "CD4 T cells",
  "CD8 T cells",
  "Treg",
  "B cells",
  "Epithelial cells",
  "Fibroblast",
  "Other immune",
  "MDSC",
  "Monocytes",
  "Dendritic",
  "Neutrophils",
  "NK cells"
)

run_cox <- function(time, status) {
  
  results <- data.frame()
  
  for (ct in celltype_vars) {
    
    varname <- make.names(ct)
    metaData[[varname]] <- densityTable[int, ct]
    
    cox_model <- coxph(
      as.formula(paste0("Surv(", time, ", ", status, ") ~ scale(", varname, ")")),
      data = metaData
    )
    
    s <- summary(cox_model)
    
    results <- rbind(
      results,
      data.frame(
        Celltype = ct,
        HR = s$coefficients[1,2],
        lower = s$conf.int[1,3],
        upper = s$conf.int[1,4],
        pvalue = s$coefficients[1,5]
      )
    )
  }
  
  results$FDR <- p.adjust(results$pvalue, method = "fdr")
  
  results
}


results_OS  <- run_cox("OS_month/12",  "OS_status")
results_RFS <- run_cox("RFS_month/12", "RFS_status")


prepare_fp <- function(df) {
  
  df <- df %>%
    mutate(
      HR_text = sprintf("%.2f (%.2f–%.2f)", HR, lower, upper),
      p_text  = ifelse(pvalue < 0.001,
                       "<0.001",
                       sprintf("%.3f", pvalue))
    ) %>%
    arrange(desc(HR))
  
  df
}

plot_fp_style <- function(df, title_text) {
  
  df <- prepare_fp(df)
  
  res <- data.frame(
    names   = df$Celltype,
    HR      = df$HR_text,
    p.value = df$p_text,
    mean    = df$HR,
    lower   = df$lower,
    upper   = df$upper
  )
  
  res |>
    forestplot(
      labeltext = c(names, HR, p.value),
      mean  = mean,
      lower = lower,
      upper = upper,
      clip = c(min(res$lower)*0.8, max(res$upper)*1.2),
      boxsize = .5,
      fn.ci_norm = "fpDrawCircleCI",
      xlog = TRUE,
      zero = 1,
      title = title_text,
      col = fpColors(box = "#660000",
                     line = "black",
                     summary = "black")
    ) |>
    
    fp_set_style(
      box = "#660000",
      line = "black",
      txt_gp = fpTxtGp(
        label = gpar(fontsize = 12),
        ticks = gpar(cex = 1),
        xlab  = gpar(cex = 1.3),
        title = gpar(fontsize = 14, fontface = "bold")
      )
    )
}

plot_fp_style(
  results_OS,
  "Overall Survival"
)

plot_fp_style(
  results_RFS,
  "Recurrence-Free Survival"
)

