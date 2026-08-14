# delete work space
rm(list = ls(all = TRUE))
graphics.off()

library(car)
library(dplyr)
library(ggplot2)

# load data
tempMetadata <- readWorkbook(
  "Y:/Projects/ICON IMC brief/Clinical data/roiMetaData_LindaRevision.xlsx",
  sheet = "revised data",
  colNames = TRUE)

roiMetadata <- tempMetadata
# roiMetadata <- roiMetadata[roiMetadata$NSCLC=="YES",]
# roiMetadata <- roiMetadata[roiMetadata$LUAD_LUSC=="LUAD",]
rownames(roiMetadata) <- roiMetadata$ROI.NO
ptID <- unique(roiMetadata$SampleID)

# Load your bulk RNA-seq data
tBulkRNA <- read.table(file = "Y:/Projects/ICON IMC/Experimental results/Integrative analysis/Data/all_Rna_samples_TMM_normalized_log_cpm.txt",
                       header = TRUE)
rownames(tBulkRNA) <- tBulkRNA[,1]
tBulkRNA <- tBulkRNA[,-c(1)]

# Bulk RNA 
tBulkRNA <- as.data.frame(t(tBulkRNA))
rownames(tBulkRNA) <- gsub("X", "", rownames(tBulkRNA))

troi <- roiMetadata[roiMetadata$TID %in% rownames(tBulkRNA),c("TID","SampleID")]
troi <- troi[!duplicated(troi$TID),]
rownames(troi) <- troi$TID

intB <- intersect(rownames(tBulkRNA),rownames(troi))
tBulkRNA <- tBulkRNA[intB,]
troi <- troi[intB,]
rownames(tBulkRNA) <- troi$SampleID


propTable <- readWorkbook(
  "Y:/Projects/ICON IMC/Experimental results/Latest section 1/Majorcelltypes differential abundance/Data/proportionMajorcelltypes.xlsx",
  sheet = "Sheet 1",
  colNames = TRUE,
  rowNames = TRUE)
propTable <- propTable[,-c(1)]


int <- intersect(rownames(tBulkRNA),ptID)

frna <- tBulkRNA[int,]

prop <- as.data.frame(matrix(0, nrow = length(ptID), ncol = ncol(propTable)))
colnames(prop) <- colnames(propTable)
rownames(prop) <- ptID 

for (ii in 1:length(ptID)){
  tidx <- roiMetadata[roiMetadata$SampleID==ptID[ii],]
  numROI <- nrow(tidx)
  
  tdata <- propTable[rownames(propTable) %in% rownames(tidx),]
  prop[ptID[ii],colnames(tdata)] <- colSums(tdata)/numROI
}

prop <- prop[complete.cases(prop),]
fprop <- prop[int,]

fprop <- 1 * logit(fprop, percents = FALSE, adjust = 0.025)


genes <- colnames(frna)

cor_summary <- data.frame(
  celltype = colnames(fprop),
  n_pos = NA,
  n_neg = NA
)

for (k in 1:ncol(fprop)) {
  
  cat("Processing", colnames(fprop)[k], "\n")
  
  cor_vals <- numeric(length(genes))
  p_vals <- numeric(length(genes))
  
  for (kk in seq_along(genes)) {
    res <- cor.test(frna[, genes[kk]], fprop[, k])
    cor_vals[kk] <- res$estimate
    p_vals[kk] <- res$p.value
  }
  
  cor_summary$n_pos[k] <- sum(cor_vals > 0 & p_vals < 0.05, na.rm = TRUE)
  cor_summary$n_neg[k] <- sum(cor_vals < 0 & p_vals < 0.05, na.rm = TRUE)
}


ggplot(cor_summary, aes(x = n_neg, y = n_pos, label = celltype)) +
  
  geom_point(size = 4, color = "#7A0177") +
  
  geom_text(nudge_y = 50, size = 4, fontface = "bold") +
  
  theme_classic(base_size = 14) +
  
  labs(
    x = "Number of negatively correlated genes",
    y = "Number of positively correlated genes"
  )


# Compute correlation between n_neg and n_pos
cor_test <- cor.test(cor_summary$n_neg, cor_summary$n_pos)

r_value <- round(cor_test$estimate, 3)
p_value <- signif(cor_test$p.value, 3)

ggplot(cor_summary, aes(x = n_neg, y = n_pos)) +
  
  # Points
  geom_point(size = 4, color = "#7A0177") +
  
  # Cell type labels
  geom_text(aes(label = celltype),
            nudge_y = max(cor_summary$n_pos) * 0.03,
            size = 4,
            fontface = "bold") +
  
  # Annotate correlation
  annotate("text",
           x = max(cor_summary$n_neg) * 0.6,
           y = max(cor_summary$n_pos) * 0.9,
           label = paste0("R = ", r_value,
                          "\nP = ", p_value),
           size = 5,
           fontface = "bold") +
  
  theme_classic(base_size = 14) +
  
  labs(
    x = "Number of negatively correlated genes",
    y = "Number of positively correlated genes"
  )

