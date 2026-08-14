rm(list = ls(all = TRUE))
graphics.off()

library(openxlsx)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(msigdbr)
library(fgsea)
library(scales)
library(qvalue)
library(car)

tempMetadata <- readWorkbook(
  "Y:/Projects/ICON IMC brief/Clinical data/roiMetaData_LindaRevision.xlsx",
  sheet = "revised data",
  colNames = TRUE)

roiMetadata <- tempMetadata

rownames(roiMetadata) <- roiMetadata$ROI.NO

ptID <- unique(roiMetadata$SampleID)

tBulkRNA <- read.table(file = "Y:/Projects/ICON IMC/Experimental results/Integrative analysis/Data/all_Rna_samples_TMM_normalized_log_cpm.txt",
                       header = TRUE)

rownames(tBulkRNA) <- tBulkRNA[,1]

tBulkRNA <- tBulkRNA[,-1]

tBulkRNA <- as.data.frame(t(tBulkRNA))

rownames(tBulkRNA) <- gsub("X", "", rownames(tBulkRNA))

troi <- roiMetadata[roiMetadata$TID %in% rownames(tBulkRNA), c("TID", "SampleID")]

troi <- troi[!duplicated(troi$TID), ]

rownames(troi) <- troi$TID

intB <- intersect(rownames(tBulkRNA), rownames(troi))

tBulkRNA <- tBulkRNA[intB,]

troi <- troi[intB,]

rownames(tBulkRNA) <- troi$SampleID

propTable <- readWorkbook(
  "Y:/Projects/ICON IMC/Experimental results/Latest section 1/Majorcelltypes differential abundance/Data/proportionMajorcelltypes.xlsx",
  sheet = "Sheet 1",
  colNames = TRUE,
  rowNames = TRUE)

propTable <- propTable[,-1]

int <- intersect(rownames(tBulkRNA), ptID)

frna <- tBulkRNA[int,]

prop <- as.data.frame(matrix(0, nrow = length(ptID), ncol = ncol(propTable)))

colnames(prop) <- colnames(propTable)

rownames(prop) <- ptID

for(ii in seq_along(ptID)){
  tidx <- roiMetadata[roiMetadata$SampleID == ptID[ii], ]
  numROI <- nrow(tidx)
  tdata <- propTable[rownames(propTable) %in% rownames(tidx), ]
  prop[ptID[ii], colnames(tdata)] <- colSums(tdata) / numROI
}

prop <- prop[complete.cases(prop), ]

fprop <- prop[int,]

fprop <- logit(fprop, percents = FALSE, adjust = 0.025)

genes <- colnames(frna)

cor_summary <- data.frame(celltype = colnames(fprop), n_pos = NA, n_neg = NA)

all_cor <- list()

for(k in seq_len(ncol(fprop))){
  cat("Processing:", colnames(fprop)[k], "\n")
  cor_vals <- numeric(length(genes))
  p_vals <- numeric(length(genes))
  
  for(kk in seq_along(genes)){
    res <- cor.test(frna[, genes[kk]], fprop[, k], method = "pearson")
    cor_vals[kk] <- res$estimate
    p_vals[kk] <- res$p.value
  }
  tmp <- data.frame(celltype = colnames(fprop)[k], gene = genes, correlation = cor_vals, pvalue = p_vals)
  all_cor[[k]] <- tmp
  cor_summary$n_pos[k] <- sum(cor_vals > 0 & p_vals < 0.05, na.rm = TRUE)
  cor_summary$n_neg[k] <- sum(cor_vals < 0 & p_vals < 0.05, na.rm = TRUE)
}

corr.df <- bind_rows(all_cor)

cor_test <- cor.test(cor_summary$n_neg, cor_summary$n_pos)

ggplot(
  cor_summary,
  aes(x = n_neg, y = n_pos)
) +
  geom_point(
    size = 4,
    color = "#7A0177"
  ) +
  geom_text(
    aes(label = celltype),
    nudge_y = max(cor_summary$n_pos)*0.03,
    size = 4,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = max(cor_summary$n_neg)*0.6,
    y = max(cor_summary$n_pos)*0.9,
    label = paste0("R = ", round(cor_test$estimate, 3), "\nP = ", signif(cor_test$p.value, 3)),
    size = 5,
    fontface = "bold"
  ) +
  theme_classic(base_size = 14)


top.genes <- corr.df %>%
  filter(correlation > 0, pvalue < 0.05) %>%
  group_by(celltype) %>%
  slice_max(correlation, n = 3) %>%
  ungroup()

cell_order <- c(
  "B.cells",
  "CD4.T.cells",
  "CD8.T.cells",
  "Dendritic",
  "Endothelial.cells",
  "Epithelial.cells",
  "Fibroblast",
  "MDSC",
  "Macrophages",
  "Monocytes",
  "NK.cells",
  "Neutrophils",
  "Other.immune",
  "Treg"
)

top.genes$celltype <- factor(top.genes$celltype, levels = rev(cell_order))

gene_order <- top.genes %>%
  group_by(gene) %>%
  summarise(max_cor = max(correlation), .groups = "drop") %>%
  arrange(max_cor) %>%
  pull(gene)

top.genes$gene <- factor(top.genes$gene, levels = gene_order)

ggplot(
  top.genes,
  aes(x = gene, y = celltype)
) +
  geom_point(
    aes(size = correlation, color = correlation)
  ) +
  scale_color_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0
  ) +
  scale_size_continuous(range = c(3, 10)) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title = element_blank(),
    panel.grid.major.y = element_line(color = "grey90")
  ) +
  labs(
    color = "Pearson r",
    size = "Pearson r"
  )


marker_df <- tribble(
  ~celltype,            ~gene,
  "Macrophages",        "C1QA",
  "Macrophages",        "C1QB",
  "Macrophages",        "APOE",
  "Macrophages",        "FCER1G",
  "Macrophages",        "TYROBP",
  "Macrophages",        "AIF1",
  "CD4.T.cells",        "CD3D",
  "CD4.T.cells",        "CD3E",
  "CD4.T.cells",        "IL7R",
  "CD4.T.cells",        "TRAC",
  "CD4.T.cells",        "LTB",
  "CD8.T.cells",        "CD8A",
  "CD8.T.cells",        "CD8B",
  "CD8.T.cells",        "TRBC1",
  "CD8.T.cells",        "GZMB",
  "CD8.T.cells",        "CCL5",
  "CD8.T.cells",        "CTSW",
  "B.cells",            "MS4A1",
  "B.cells",            "CD79A",
  "B.cells",            "CD74",
  "B.cells",            "BANK1",
  "Fibroblast",         "COL1A1",
  "Fibroblast",         "COL1A2",
  "Fibroblast",         "COL3A1",
  "Fibroblast",         "DCN",
  "Fibroblast",         "LUM",
  "Fibroblast",         "C7",
  "Fibroblast",         "CFD",
  "Monocytes",          "FCN1",
  "Monocytes",          "VCAN",
  "Monocytes",          "S100A12",
  "Monocytes",          "CTSS",
  "Monocytes",          "SAT1",
  "Monocytes",          "TYMP",
  "Dendritic",          "FCER1A",
  "Dendritic",          "HLA-DRA",
  "Dendritic",          "CLEC10A",
  "Dendritic",          "CD74",
  "Dendritic",          "CLEC9A",
  "Treg",               "FOXP3",
  "Treg",               "IL2RA",
  "Treg",               "TIGIT",
  "Treg",               "CTLA4",
  "NK.cells",           "NKG7",
  "NK.cells",           "GNLY",
  "NK.cells",           "FCGR3A",
  "NK.cells",           "KLRD1",
  "NK.cells",           "TYROBP",
  "NK.cells",           "CTSW",
  "Endothelial.cells",  "KDR",
  "Endothelial.cells",  "PECAM1",
  "Endothelial.cells",  "VWF",
  "Endothelial.cells",  "EMCN",
  "Endothelial.cells",  "CDH5",
  "Epithelial.cells",   "EPCAM",
  "Epithelial.cells",   "KRT8",
  "Epithelial.cells",   "KRT18",
  "Epithelial.cells",   "KRT19",
  "Epithelial.cells",   "TACSTD2",
  "Neutrophils",        "S100A8",
  "Neutrophils",        "S100A9",
  "Neutrophils",        "FCGR3B",
  "Neutrophils",        "CXCR2",
  "Neutrophils",        "CSF3R"
)

plot.df <- corr.df %>%
  inner_join(marker_df, by = c("celltype", "gene")) %>%
  filter(!is.na(correlation))

plot.df <- plot.df %>%
  mutate(neglog10P = -log10(pvalue))

gene_order <- c(
  "C1QA","C1QB","APOE","FCER1G","TYROBP","AIF1",
  "CD3D","CD3E","IL7R","TRAC","LTB",
  "CD8A","CD8B","TRBC1","GZMB","CCL5","CTSW",
  "MS4A1","CD79A","CD74","BANK1",
  "COL1A1","COL1A2","COL3A1","DCN","LUM","C7","CFD",
  "FCN1","VCAN","S100A12","CTSS","SAT1","TYMP",
  "FCER1A","HLA-DRA","CLEC10A","CD74","CLEC9A",
  "FOXP3","IL2RA","TIGIT","CTLA4",
  "NKG7","GNLY","FCGR3A","KLRD1","TYROBP","CTSW",
  "KDR","PECAM1","VWF","EMCN","CDH5",
  "EPCAM","KRT8","KRT18","KRT19","TACSTD2",
  "S100A8","S100A9","FCGR3B","CXCR2","CSF3R"
)

gene_order <- gene_order[gene_order %in% plot.df$gene]

gene_order <- unique(gene_order)

plot.df$gene <- factor(plot.df$gene, levels = gene_order)

cell_order <- c(
  "B.cells",
  "CD4.T.cells",
  "CD8.T.cells",
  "Dendritic",
  "Endothelial.cells",
  "Epithelial.cells",
  "Fibroblast",
  "Macrophages",
  "Monocytes",
  "NK.cells",
  "Neutrophils",
  "Treg"
)

plot.df$celltype <- factor(plot.df$celltype, levels = rev(cell_order))

ggplot(plot.df, aes(x = gene, y = celltype)) +
  geom_point(
    aes(color = correlation, size = neglog10P),
    alpha = 0.9
  ) +
  scale_color_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = "Pearson r"
  ) +
  scale_size_continuous(range = c(2, 10), name = expression(-log[10](P))) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold")
  ) +
  labs(
    title = "Canonical lineage markers validate IMC-derived cellular phenotypes"
  )


cell_order <- c(
  "Macrophages",
  "CD4.T.cells",
  "CD8.T.cells",
  "B.cells",
  "Fibroblast",
  "Monocytes",
  "Dendritic",
  "Treg",
  "NK.cells",
  "Endothelial.cells",
  "Epithelial.cells",
  "Neutrophils"
)

plot.df$celltype_tmp <- factor(plot.df$celltype, levels = cell_order)

gene_order <- plot.df %>%
  group_by(celltype_tmp) %>%
  arrange(desc(abs(correlation)), .by_group = TRUE) %>%
  ungroup() %>%
  pull(gene) %>%
  unique()

plot.df$gene <- factor(plot.df$gene, levels = gene_order)

plot.df$celltype <- factor(plot.df$celltype, levels = cell_order)

ggplot(plot.df, aes(x = gene, y = celltype)) +
  geom_point(
    aes(color = correlation, size = neglog10P),
    alpha = 0.9
  ) +
  scale_color_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = "Pearson r"
  ) +
  scale_size_continuous(range = c(2, 10), name = expression(-log[10](P))) +
  scale_y_discrete(limits = rev) +
  theme_classic(
    base_size = 16
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold")
  ) +
  labs(
    title =
      "Canonical lineage markers validate IMC-derived cellular phenotypes"
  )



signature_list <- list(
  Macrophages = c("C1QA", "C1QB", "APOE", "FCER1G", "AIF1"),
  CD4.T.cells = c("CD3D", "CD3E", "IL7R", "TRAC", "LTB"),
  CD8.T.cells = c("CD8A", "CD8B", "TRBC1", "GZMB", "CCL5"),
  B.cells = c("MS4A1", "CD79A", "BANK1"),
  Fibroblast = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "C7"),
  Monocytes = c("FCN1", "VCAN", "S100A12", "CTSS", "SAT1"),
  Dendritic = c("FCER1A", "HLA-DRA", "CLEC10A", "CLEC9A"),
  Treg = c("FOXP3", "IL2RA", "TIGIT", "CTLA4"),
  NK.cells = c("NKG7", "GNLY", "FCGR3A", "KLRD1"),
  Endothelial.cells = c("KDR", "PECAM1", "VWF", "EMCN", "CDH5"),
  Epithelial.cells = c("EPCAM", "KRT8", "KRT18", "KRT19"),
  Neutrophils = c("S100A8", "S100A9", "FCGR3B", "CXCR2")
)

signature_scores <- data.frame(row.names = rownames(frna))

for(ct in names(signature_list)){
  genes_use <- intersect(signature_list[[ct]], colnames(frna))
  cat(ct, ":", length(genes_use), "genes found\n")
  if(length(genes_use) < 2){
    cat("Skipping:", ct, "\n")
    next
  }
  signature_scores[, ct] <- rowMeans(frna[, genes_use, drop=FALSE])
}

sig_cor <- data.frame()

for(ct in intersect(colnames(signature_scores), colnames(fprop))){
  cat("Correlating:", ct, "\n")
  res <- cor.test(signature_scores[, ct], fprop[, ct], method = "pearson")
  sig_cor <- rbind(sig_cor, data.frame(Celltype = ct, Correlation = as.numeric(res$estimate), Pvalue = res$p.value))
}

sig_cor$FDR <- p.adjust(sig_cor$Pvalue, method = "BH")

sig_cor$neglog10P <- -log10(sig_cor$Pvalue)

sig_cor$Celltype <- factor(sig_cor$Celltype, levels = sig_cor$Celltype[order(sig_cor$Correlation)])
sig_cor %>%
  arrange(desc(Correlation))

sig_cor <- sig_cor %>%
  arrange(Correlation)

sig_cor$Celltype <- factor(sig_cor$Celltype, levels = sig_cor$Celltype)

sig_cor$label <- ifelse(
  sig_cor$Pvalue < 0.001,
  "***",
  ifelse(sig_cor$Pvalue < 0.01, "**", ifelse(sig_cor$Pvalue < 0.05, "*", ""))
)

ggplot(
  sig_cor,
  aes(y = Celltype, x = Correlation)
) +
  geom_segment(
    aes(x = 0, xend = Correlation, y = Celltype, yend = Celltype),
    linewidth = 1.2,
    color = "grey70"
  ) +
  geom_point(
    aes(size = -log10(Pvalue)),
    color = "#2166AC"
  ) +
  geom_text(
    aes(label = label),
    nudge_x = 0.025,
    size = 6,
    fontface = "bold"
  ) +
  geom_vline(xintercept = 0, linetype = 2, color = "black") +
  scale_size_continuous(range = c(4, 10)) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.y = element_text(face = "bold"),
    legend.position = "right"
  ) +
  labs(
    x = "Pearson correlation",
    y = NULL,
    color = "Pearson r",
    size = expression(-log[10](P)),
    title = "Bulk RNA lineage signatures validate IMC phenotypes"
  )
