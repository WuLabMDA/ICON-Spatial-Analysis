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
library(stringr)

load("Y:/Projects/ICON IMC/Final/S3 cellular neighborhood/Data/CN_workspace.rds")

# call either normal or tumor
unq <- rownames(densityTable)
for (i in 1:length(unq)){
  densityTable[unq[i], "Location"] <- roiMetadata[unq[i],"Location"]
  freqTable[unq[i], "Location"] <- roiMetadata[unq[i],"Location"]
  propTable[unq[i], "Location"] <- roiMetadata[unq[i],"Location"]
}


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


rna <- readWorkbook(
  "Y:/Projects/ICON IMC brief/Data/GeoMx_Expr_NormalTumor.xlsx",
  sheet = "Sheet 1",
  colNames = TRUE,
  rowNames = TRUE)


tempTCR <- readWorkbook(
  "Y:/Projects/ICON IMC/Experimental results/TCR section/Data/20240319_ICON_all.xlsx",
  sheet = "Normal_TCR",
  colNames = TRUE)

TCR <- tempTCR[,c(1,6,7)]
colnames(TCR) <- c("SampleID","Richness","Clonality")


roi_location <- "Tumor"
# roiMetadata <- roiMetadata[roiMetadata$NSCLC=="YES",]
roiMetadata <- roiMetadata[roiMetadata$Location==roi_location,]
rownames(roiMetadata) <- roiID <- roiMetadata$ROI.NO

densityTable <- densityTable[densityTable$Location == roi_location,c(1:8)]

freqTable <- freqTable[freqTable$Location == roi_location,c(1:8)]

propTable <- propTable[propTable$Location == roi_location,c(1:8)]


pTID <- roiMetadata[,c(1,2,6)]


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

rna$TID <- pTID$TID[match(rownames(rna), pTID$SampleID)]

frna <- rna[rna$TID %in% int,]


CN_long <- propTable %>%
  tibble::rownames_to_column("ROI") %>%
  pivot_longer(cols = -ROI,
               names_to = "CN",
               values_to = "Proportion")

CN_mean <- colMeans(propTable, na.rm = TRUE)


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

fprop$CN1_group <- ifelse(fprop$`CN 3` > median(fprop$`CN 3`),
                          "High","Low")

fprop$CN1_group <- factor(fprop$CN1_group, levels = c("Low","High"))


fprop$CN1_group <- factor(
  fprop$CN1_group,
  levels = c("Low", "High")
)

pval <- wilcox.test(
  Clonality ~ CN1_group,
  data = fprop
)$p.value

label <- paste0(
  "Wilcoxon P = ",
  format.pval(pval, digits = 2, eps = 1e-4)
)

ggplot(
    fprop,
    aes(
      x = CN1_group,
      y = Clonality,
      fill = CN1_group
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
  aes(colour = CN1_group),
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
  title = "TLS-like neighborhoods are associated with increased T-cell clonal expansion",
  x = "TLS-like CN abundance",
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


frna$CN1_group <- fprop$CN1_group[match(frna$TID,rownames(fprop))]
frna$Clonality <- fprop$Clonality[match(frna$TID,rownames(fprop))]
frna$Richness <- fprop$Richness[match(frna$TID,rownames(fprop))]

TLS_genes <- c(
  "CXCL13",
  "CXCR5",
  "CCL19",
  "CCL21",
  "LTA",
  "LTB",
  "CD79A",
  "CD79B",
  "MS4A1",
  "BANK1",
  "TNFRSF13C",
  "ICOS"
)

TLS_genes <- TLS_genes[TLS_genes %in% colnames(frna)]

frna$TLS_score <- rowMeans(frna[, TLS_genes], na.rm = TRUE)


frna$CN1_group <- factor(
  frna$CN1_group,
  levels = c("Low", "High")
)

pval <- wilcox.test(
  TLS_score ~ CN1_group,
  data = frna
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
    frna,
    aes(
      x = CN1_group,
      y = TLS_score,
      fill = CN1_group
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
  aes(colour = CN1_group),
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
  y = max(frna$TLS_score, na.rm = TRUE) * 1.08,
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
  title = "TLS transcriptional activity is enriched in TLS-like CN-high tumors",
  x = "IMC TLS-like CN abundance",
  y = "GeoMx TLS signature score"
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


genes <- setdiff(colnames(frna),
                 c("TID","CN1_group","TLS_score","Clonality","Richness"))

DE_list <- lapply(genes, function(g){
  
  test <- wilcox.test(frna[[g]] ~ frna$CN1_group)
  
  data.frame(
    gene = g,
    pvalue = test$p.value,
    logFC = mean(frna[[g]][frna$CN1_group=="High"], na.rm=TRUE) -
      mean(frna[[g]][frna$CN1_group=="Low"], na.rm=TRUE)
  )
})

DE <- bind_rows(DE_list)
DE$FDR <- p.adjust(DE$pvalue, method="fdr")
DE <- DE[order(DE$FDR),]

high_genes <- DE %>%
  filter(pvalue < 0.05 & logFC > 0) %>%
  pull(gene)

low_genes <- DE %>%
  filter(pvalue < 0.05 & logFC < 0) %>%
  pull(gene)

# high_genes <- DE %>%
#   filter(logFC > 0) %>%
#   arrange(pvalue) %>%
#   dplyr::slice(1:100) %>%
#   pull(gene)
# 
# low_genes <- DE %>%
#   filter(logFC < 0) %>%
#   arrange(pvalue) %>%
#   dplyr::slice(1:100) %>%
#   pull(gene)


#library(gprofiler2)

# gost_high <- gost(
#   query = high_genes,
#   organism = "hsapiens",
#   sources = c("GO:BP","KEGG","REAC")
# )
# 
# gost_low <- gost(
#   query = low_genes,
#   organism = "hsapiens",
#   sources = c("GO:BP","KEGG","REAC")
# )

gost_high <- gost(
  query = high_genes,
  organism = "hsapiens",
  sources = c(
    "GO:BP",
    "GO:MF",
    "GO:CC",
    "KEGG",
    "REAC"
  )
)

gost_low <- gost(
  query = low_genes,
  organism = "hsapiens",
  sources = c(
    "GO:BP",
    "GO:MF",
    "GO:CC",
    "KEGG",
    "REAC"
  )
)


pathways_interest <- c(
  "B cell activation",
  "B cell differentiation",
  "B cell receptor signaling pathway",
  "T cell activation",
  "alpha-beta T cell activation",
  "T cell receptor signaling pathway",
  "T cell costimulation",
  "T cell proliferation",
  "Cytokine Signaling in Immune system",
  "Signaling by Interleukins",
  "NF-kappa B signaling pathway",
  "Interleukin-2 signaling",
  "Interleukin-15-mediated signaling pathway",
  "immune response",
  "immune effector process",
  "antigen receptor-mediated signaling pathway"
)

plot_df <- gost_high$result[
  gost_high$result$term_name %in% pathways_interest,]

plot_df <- plot_df[-c(12),]



ggplot(plot_df,
       aes(x=reorder(term_name,-log10(p_value)),
           y=-log10(p_value),
           size=intersection_size)) +
  
  geom_point(color="#d95f02") +
  
  coord_flip() +
  
  theme_classic(base_size=14) +
  
  labs(
    x="Pathway",
    y="-log10(p-value)",
    size="Gene count",
    title="Adaptive immune signaling enriched only in TLS-high tumors"
  )


celltypes <- colnames(prop)

res_list <- lapply(celltypes, function(ct){
  
  df <- fprop[,c(ct,"Clonality","Richness")]
  
  cor_clon <- cor.test(df[[ct]], df$Clonality, method="spearman")
  cor_rich <- cor.test(df[[ct]], df$Richness, method="spearman")
  
  data.frame(
    Celltype = ct,
    metric = c("Clonality","Richness"),
    rho = c(cor_clon$estimate, cor_rich$estimate),
    pval = c(cor_clon$p.value, cor_rich$p.value)
  )
})


celltypes <- colnames(fprop)[!colnames(fprop) %in% c("Clonality","Richness")]

form <- as.formula(
  paste("Clonality ~", paste(sprintf("`%s`", celltypes), collapse = " + "))
)

fit <- lm(form, data = fprop)

summary(fit)


genes <- setdiff(
  colnames(frna),
  c("TID","CN1_group","TLS_score","Clonality","Richness")
)

rank.df <- lapply(genes, function(g){
  
  test <- wilcox.test(
    frna[[g]] ~ frna$CN1_group
  )
  
  data.frame(
    gene = g,
    logFC =
      mean(frna[[g]][frna$CN1_group=="High"], na.rm=TRUE) -
      mean(frna[[g]][frna$CN1_group=="Low"], na.rm=TRUE),
    pvalue = test$p.value
    
  )
  
})

rank.df <- bind_rows(rank.df)

rank.df$stat <-
  sign(rank.df$logFC) *
  -log10(rank.df$pvalue)

gene.ranks <- rank.df$stat

names(gene.ranks) <- rank.df$gene

gene.ranks <- sort(
  gene.ranks,
  decreasing = TRUE
)

library(msigdbr)

hallmark <- msigdbr(
  species = "Homo sapiens",
  category = "C2"
)

go.bp <- msigdbr(
  species="Homo sapiens",
  category="C5",
  subcategory="GO:BP"
)

immune.pathways <-
  go.bp %>%
  filter(
    grepl(
      "T_CELL|T_CELL_RECEPTOR|B_CELL|CYTOTOXIC|INTERLEUKIN_2|ADAPTIVE",
      gs_name,
      ignore.case=TRUE))

pathways <- split(
    immune.pathways$gene_symbol,
    immune.pathways$gs_name)

library(fgsea)
        
fg <- fgsea(
  pathways = pathways,
  stats = gene.ranks,
  minSize = 10,
  maxSize = 500)

fg <- fg %>% arrange(desc(NES))

keep <- c(
  "T_CELL_RECEPTOR_SIGNALING_PATHWAY",
  "T_CELL_ACTIVATION",
  "ALPHA_BETA_T_CELL_ACTIVATION",
  "T_CELL_PROLIFERATION",
  "POSITIVE_REGULATION_OF_T_CELL_ACTIVATION",
  "B_CELL_ACTIVATION",
  "B_CELL_RECEPTOR_SIGNALING_PATHWAY",
  "B_CELL_DIFFERENTIATION",
  "INTERLEUKIN_2_MEDIATED_SIGNALING_PATHWAY",
  "CYTOTOXIC_T_CELL_DIFFERENTIATION",
  "IMMUNE_EFFECTOR_PROCESS")


plot.df2 <- fg %>%
  mutate(
    Program = case_when(
      str_detect(pathway, regex("T_CELL_RECEPTOR", ignore_case = TRUE)) ~
        "TCR signaling",
      
      str_detect(pathway, regex("INTERLEUKIN_2|IL2", ignore_case = TRUE)) ~
        "IL-2 signaling",
      
      str_detect(pathway, regex("T_CELL_ACTIVATION", ignore_case = TRUE)) ~
        "T-cell activation",
      
      str_detect(pathway, regex("T_CELL_PROLIFERATION", ignore_case = TRUE)) ~
        "T-cell proliferation",
      
      str_detect(pathway, regex("B_CELL_RECEPTOR", ignore_case = TRUE)) ~
        "B-cell receptor signaling",
      
      str_detect(pathway, regex("B_CELL_ACTIVATION", ignore_case = TRUE)) ~
        "B-cell activation",
      
      str_detect(pathway, regex("B_CELL_DIFFERENTIATION", ignore_case = TRUE)) ~
        "B-cell differentiation",
      
      str_detect(pathway, regex("IMMUNE_EFFECTOR", ignore_case = TRUE)) ~
        "Cytotoxic immunity",
      
      TRUE ~ NA_character_
    )
    
  ) %>% filter(!is.na(Program))

plot.df2 <- plot.df2 %>%
  group_by(Program) %>%
  slice_max(
    NES,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

plot.df2$Program <- factor(
    plot.df2$Program,
    levels =
      plot.df2$Program[
        order(plot.df2$NES)])


# plot.df2 <- plot.df2 %>%
#   mutate(
#     Program = recode(
#       Program,
#       "T-cell receptor signaling" = "TCR signaling",
#       "B-cell receptor signaling" = "BCR signaling"
#     ))

plot.df2$Program <- factor(
  plot.df2$Program,
  levels = c(
    "TCR signaling",
    "T-cell activation",
    "T-cell proliferation",
    "IL-2 signaling",
    "B-cell activation",
    "B-cell differentiation",
    "BCR signaling",
    "B-cell receptor signaling"
  ))

plot.df2$Program <- factor(
  plot.df2$Program,
  levels = plot.df2$Program[order(plot.df2$NES, decreasing = FALSE)]
)

plot.df2$stars <- case_when(
  plot.df2$padj < 0.001 ~ "***",
  plot.df2$padj < 0.01  ~ "**",
  plot.df2$padj < 0.05  ~ "*",
  TRUE                  ~ ""
)

ggplot(
  plot.df2,
  aes(
    x = NES,
    y = Program
  )
) +
  
geom_vline(
  xintercept = 0,
  linetype = 2,
  colour = "grey85"
) +
geom_segment(
  aes(
    x = 0,
    xend = NES,
    y = Program,
    yend = Program
  ),
  colour = "grey88",
  linewidth = 0.8
) +
  
geom_point(
  aes(
    size = size,
    colour = -log10(padj)
  ),
  alpha = 0.95
) +
geom_text(
  aes(
    x = NES + 0.015,
    label = stars
  ),
  fontface = "bold",
  size = 6
) +
scale_colour_gradient2(
  low = "#2C7BB6",
  mid = "#F7F7F7",
  high = "#D7191C",
  midpoint = median(-log10(plot.df2$padj)),
  name = expression(-log[10](FDR))
) +
scale_size_continuous(
  range = c(4,8),
  name = "Gene count"
) +
coord_cartesian(
  xlim = c(1.55, 2.1)
) +
labs(
  title = "Adaptive immune programs enriched in TLS-high tumors",
  x = "Normalized enrichment score (NES)",
  y = NULL
) +
theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 18
    ),
    axis.text.y = element_text(
      face = "bold",
      size = 13,
      colour = "black"
    ),
    
    axis.text.x = element_text(
      face = "bold",
      size = 12,
      colour = "black"
    ),
    
    axis.title.x = element_text(
      face = "bold",
      size = 14
    ),
    
    legend.title = element_text(
      face = "bold",
      size = 13
    ),
    
    legend.text = element_text(
      face = "bold",
      size = 12
    )
    
  )



