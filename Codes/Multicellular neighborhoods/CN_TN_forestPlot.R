rm(list = ls(all = TRUE))
graphics.off()

library(ggplot2)
library(dplyr)
library(openxlsx)
library(survival)
library(survminer)
library(forestplot)
library(grid)
library(forcats)


load("Y:/Projects/ICON IMC/Final/S3 cellular neighborhood/Data/CN_workspace.rds")

ptID <- unique(tempdata@colData@listData[["patient_id"]])

# call either normal or tumor
roiLocation <- "Normal"
tdata <- data[,colData(data)$location==roiLocation]
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

celltypes <- as.character(unique(data@colData@listData[["cn_celltypes"]]))

densityTable <- as.data.frame(matrix(0, nrow = length(unq1), ncol = length(celltypes)))
colnames(densityTable) <- celltypes
rownames(densityTable) <- unq1 
# unq <- unique(data@colData@listData[["sample_id"]])

for (ii in 1:length(unq1)){
  message("Processing: ", ii, " of ", length(unq1))
  tidx <- roiMetadata[roiMetadata$SampleID==unq1[ii],]
  numROI <- nrow(tidx)
  fdata <- matrix(0, nrow = 1, ncol = length(celltypes))
  colnames(fdata) <- celltypes
  for (j in 1:numROI){
    tdata <- data[,grepl(tidx$ROI.NO[j], data@colData@listData[["sample_id"]])]
    df <- as.numeric()
    for (k in 1:length(celltypes)){
      df[k] <- sum(tdata@colData@listData[["cn_celltypes"]]==celltypes[k])
    }
    fdata <- fdata + ((df*1e6)/tidx$Area[j])
  }
  densityTable[unq1[ii],celltypes] <- fdata/numROI
}

int <- intersect(rownames(densityTable),rownames(metaData))
metaData <- metaData[int,]


celltype_vars <- c(
  "CN 1",
  "CN 2",
  "CN 3",
  "CN 4",
  "CN 5",
  "CN 6",
  "CN 7",
  "CN 8"
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
      boxsize = .3,
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


##############

metaData$CN3_high <- ifelse(
  metaData$CN3 > median(metaData$CN3, na.rm=TRUE),
  "High", "Low"
)

fit <- survfit(Surv(RFS_month/12, RFS_status) ~ CN3_high, data=metaData)

ggsurvplot(fit, data=metaData, pval = TRUE)

############
int <- intersect(rownames(densityTable),rownames(metaData))
metaData <- metaData[int,]
metaData$CN1 <- densityTable[int,c("CN 1")]
metaData$CN2 <- densityTable[int,c("CN 2")]
metaData$CN3 <- densityTable[int,c("CN 3")]
metaData$CN4 <- densityTable[int,c("CN 4")]
metaData$CN5 <- densityTable[int,c("CN 5")]
metaData$CN6 <- densityTable[int,c("CN 6")]
metaData$CN7 <- densityTable[int,c("CN 7")]
metaData$CN8 <- densityTable[int,c("CN 8")]



metaData$CN1[densityTable$`CN 1`<= mean(densityTable$`CN 1`)] <- "Low density"
metaData$CN1[densityTable$`CN 1`> mean(densityTable$`CN 1`)] <- "High density"

metaData$CN2[densityTable$`CN 2` <= mean(densityTable$`CN 2`)] <- "Low density"
metaData$CN2[densityTable$`CN 2`> mean(densityTable$`CN 2`)] <- "High density"

metaData$CN3[densityTable$`CN 3` <= mean(densityTable$`CN 3`)] <- "Low density"
metaData$CN3[densityTable$`CN 3`> mean(densityTable$`CN 3`)] <- "High density"

metaData$CN4[densityTable$`CN 4` <= mean(densityTable$`CN 4`)] <- "Low density"
metaData$CN4[densityTable$`CN 4`> mean(densityTable$`CN 4`)] <- "High density"

metaData$CN5[densityTable$`CN 5` <= mean(densityTable$`CN 5`)] <- "Low density"
metaData$CN5[densityTable$`CN 5`> mean(densityTable$`CN 5`)] <- "High density"

metaData$CN6[densityTable$`CN 6` <= mean(densityTable$`CN 6`)] <- "Low density"
metaData$CN6[densityTable$`CN 6`> mean(densityTable$`CN 6`)] <- "High density"

metaData$CN7[densityTable$`CN 7` <= mean(densityTable$`CN 7`)] <- "Low density"
metaData$CN7[densityTable$`CN 7`> mean(densityTable$`CN 7`)] <- "High density"

metaData$CN8[densityTable$`CN 8` <= mean(densityTable$`CN 8`)] <- "Low density"
metaData$CN8[densityTable$`CN 8`> mean(densityTable$`CN 8`)] <- "High density"


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
               ~ CN1, data = metaData)

summary(coxph(Surv(OS_month/12,OS_status)  ~ CN1, data = metaData))

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
                ~ CN8, data = metaData)

summary(coxph(Surv(RFS_month/12,RFS_status)  ~ CN8, data = metaData))

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

