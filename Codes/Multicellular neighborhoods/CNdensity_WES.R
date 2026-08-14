rm(list = ls(all = TRUE))
graphics.off()

library(openxlsx)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(caret)
library(pROC)


load("Y:/Projects/ICON IMC/Final/S3 cellular neighborhood/Data/CN_workspace.rds")

# call either normal or tumor
unq <- rownames(densityTable)
for (i in 1:length(unq)){
  densityTable[unq[i], "Location"] <- roiMetadata[unq[i],"Location"]
  freqTable[unq[i], "Location"] <- roiMetadata[unq[i],"Location"]
  propTable[unq[i], "Location"] <- roiMetadata[unq[i],"Location"]
}

densityTable <- densityTable[densityTable$Location == "Tumor",c(1:8)]

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
cdata[,9:ncol(cdata)] <- lapply(cdata[,9:ncol(cdata)], factor)
features <- colnames(cdata)

corr <- pvalues <- as.data.frame(matrix(0, nrow = ncol(prop), ncol = ncol(fGenes)))
rownames(corr) <- rownames(pvalues) <- colnames(prop)
colnames(corr) <- colnames(pvalues) <- colnames(fGenes)

for (i in 1:ncol(fGenes)){
  for (j in 1:ncol(prop)){
    model <- lm(cdata[[j]] ~ cdata[[i+8]], data = cdata)
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
  "CDKN2B",
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
  "FBXW7",
  "AR",
  "DDR2",
  "BRCA1",
  "ERBB2",
  "NF1",
  "APC"
)

fdata_sub <- fdata %>%
  filter(genes %in% selected_genes) %>%
  mutate(
    genes = factor(genes, levels = selected_genes),
    logP = -log10(Pvalues)
  )


colGEX <- c("grey90", brewer.pal(7, "Reds"))

ggplot(fdata_sub, aes(genes, celltypes)) +
  
  geom_point(
    data = subset(fdata_sub, Pvalues >= 0.05),
    aes(size = logP, fill = Rsquared),
    shape = 21,
    color = "lightgrey",
    alpha = 1
  ) +
  
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
    title = "",
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


#############
# classification
library(caret)
library(pROC)
library(dplyr)

cn_vars <- c("CN 1","CN 2","CN 3","CN 4",
             "CN 5","CN 6","CN 7","CN 8")

mutation_vars <- colnames(cdata)[9:49]


cdata_scaled <- cdata
cdata_scaled[, cn_vars] <- scale(cdata_scaled[, cn_vars])

cdata_scaled[, mutation_vars] <- 
  lapply(cdata_scaled[, mutation_vars], function(x)
    factor(x, levels = c(0,1), labels = c("WT","Mut")))

ctrl <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 10,
  classProbs = TRUE,
  summaryFunction = twoClassSummary
)

auc_results <- data.frame()

for (cn in cn_vars) {
  
  for (mut in mutation_vars) {
    message("Using: ",cn, " 2predict ",mut)
    
    formula_text <- paste0(
      mut, " ~ `", cn, "`"
    )
    
    model <- train(
      as.formula(formula_text),
      data = cdata_scaled,
      method = "glm",
      family = binomial,
      metric = "ROC",
      trControl = ctrl
    )
    
    auc_results <- rbind(
      auc_results,
      data.frame(
        CN = cn,
        Mutation = mut,
        AUC = model$results$ROC
      )
    )
  }
}

library(tidyr)
library(ggplot2)

auc_mat <- auc_results %>%
  pivot_wider(names_from = Mutation, values_from = AUC)

auc_long <- auc_results

ggplot(auc_long,
       aes(Mutation, CN, fill = AUC)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "white",
    mid = "#FDB863",
    high = "#B2182B",
    midpoint = 0.6,
    limits = c(0.5, max(auc_long$AUC))
  ) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "AUC")


cnColors <- c("#17BECF",  
              "#33A02C",  
              "#FF7F00",  
              "#e0ab09",  
              "#B15928",
              "#8d67f5",
              "#E31A1C",
              "#1F78B4")

cncelltypes <- setNames(cnColors, c("CN 1", "CN 2", "CN 3", "CN 4", "CN 5", "CN 6", "CN 7", "CN 8"))


ggplot(auc_long,
       aes(x = AUC,
           y = reorder(Mutation, AUC),
           color = CN)) +
  
  geom_point(size = 3.2, alpha = 0.95) +
  
  geom_vline(xintercept = 0.5,
             linetype = "dashed",
             color = "grey50",
             linewidth = 0.6) +
  
  scale_color_manual(values = cncelltypes) +
  
  theme_classic(base_size = 13) +
  
  theme(
    axis.text.y = element_text(size = 9, face = "bold"),
    axis.text.x = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  ) +
  
  labs(
    x = "Cross-validated AUC",
    y = "Mutation",
    color = "Cellular Neighborhood"
  )


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
  "NFE2L2",
  "FGFR1",
  "RET",
  "GNAS",
  "FBXW7"
)


auc_selected <- auc_long %>%
  filter(Mutation %in% selected_genes)

auc_selected$Mutation <- factor(
  auc_selected$Mutation,
  levels = selected_genes
)

auc_selected <- auc_selected %>%
  group_by(Mutation) %>%
  mutate(maxAUC = max(AUC)) %>%
  ungroup()

auc_selected$Mutation <- reorder(
  auc_selected$Mutation,
  auc_selected$maxAUC
)

ggplot(auc_selected,
       aes(x = AUC,
           y = Mutation,
           color = CN)) +
  
  geom_point(size = 3, alpha = 0.9) +
  
  geom_vline(xintercept = 0.5,
             linetype = "dashed",
             color = "grey50") +
  
  scale_color_manual(values = cncelltypes) +
  
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.2),
    limits = c(0, 1)
  ) +
  
  theme_classic(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 10, face = "bold"),
    legend.position = "right"
  ) +
  
  labs(x = "AUC",
       y = "Mutation",
       color = "CN")

