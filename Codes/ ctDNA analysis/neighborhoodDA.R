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
library(tidyr)
library(MASS)

load("Y:/Projects/ICON IMC/Final/S3 cellular neighborhood/Data/CN_workspace.rds")


# load data
tempMetadata <- readWorkbook(
  "Y:/Projects/ICON IMC/Clinical data/roiMetadata2.xlsx",
  sheet = "Sheet 1",
  colNames = TRUE)
tempMetadata <- tempMetadata[,-1]

roiMetadata <- tempMetadata
rownames(roiMetadata) <- roiMetadata$ROI.NO

location <- "Normal"
roiMetadata <- roiMetadata[roiMetadata$Location == location,]


# load data
ctdnaData <- readWorkbook(
  "Y:/Projects/ICON IMC/Experimental results/ctDNA section/Data/ctDNA_Sorted.xlsx",
  sheet = "combined",
  colNames = TRUE,
  rowNames = TRUE)

ptID <- intersect(rownames(ctdnaData),unique(roiMetadata$SampleID))

dens <- as.data.frame(matrix(0, nrow = length(ptID), ncol = ncol(densityTable)))
colnames(dens) <- colnames(densityTable)
rownames(dens) <- ptID 

for (ii in 1:length(ptID)){
  tidx <- roiMetadata[roiMetadata$SampleID==ptID[ii],]
  numROI <- nrow(tidx)
  
  tdata <- densityTable[rownames(densityTable) %in% rownames(tidx),]
  dens[ptID[ii],colnames(tdata)] <- colSums(tdata)/numROI
}

dens <- dens[complete.cases(dens),]
dens$ctdna <- ctdnaData[rownames(dens),"ctDNA_status"]

densityTable_avg <- dens %>%
  mutate(
    ctdna = case_when(
      ctdna == "0" ~ "Not detected",
      ctdna == "1" ~ "Detected",
      TRUE ~ NA_character_
    )
  )


da_long <- densityTable_avg %>%
  pivot_longer(
    cols = 1:8,                
    names_to = "CN",
    values_to = "Density"
  ) %>%
  mutate(
    ctdna = factor(ctdna, levels = c("Detected","Not detected"))
  )


nb_results <- lapply(split(da_long, da_long$CN), function(df) {
  
  fit <- glm.nb(
    Density ~ ctdna,
    data = df
  )
  
  coef_summary <- summary(fit)$coefficients
  
  data.frame(
    CN = unique(df$CN),
    logFC    = coef_summary["ctdnaNot detected", "Estimate"],
    p_value  = coef_summary["ctdnaNot detected", "Pr(>|z|)"]
  )
}) %>% bind_rows()



plot_df <- nb_results %>%
  mutate(
    pSignif = p_value < 0.05,
    CN = factor(CN)
  )

plot_df <- plot_df %>%
  arrange(logFC) %>%
  mutate(CN = factor(CN, levels = CN))


ggbarplot(
  plot_df,
  x = "CN",
  y = "logFC",
  width = 1,
  fill = "pSignif",
  color = "white",
  palette = c("ivory4", "aquamarine4"),
  sort.val = "asc",
  sort.by.groups = FALSE,
  x.text.angle = 90,
  ylab = "log2FC",
  xlab = FALSE,
  legend.title = "p < 0.05"
) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1),
    axis.text   = element_text(face = "bold"),
    axis.title  = element_text(face = "bold")
  ) +
  coord_flip()


plot_df <- nb_results %>%
  mutate(
    Direction = case_when(
      logFC > 0  ~ "Not detected",
      logFC < 0  ~ "Detected",
      TRUE       ~ "No change"
    ),
    signif_label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ ""
    )
  ) %>%
  arrange(logFC) %>%
  mutate(
    CN = factor(CN, levels = CN),
    star_pos = ifelse(logFC > 0, logFC + 0.05, logFC - 0.05)
  )


ggplot(plot_df,
       aes(x = CN, y = logFC, fill = Direction)) +
  
  geom_col(width = 0.9) +
  
  geom_text(
    aes(y = star_pos, label = signif_label),
    size = 5,
    fontface = "bold",
    hjust = ifelse(plot_df$logFC > 0, 0, 1)
  ) +
  
  scale_fill_manual(values = c(
    "Detected" = "#2C7BB6",   # calm blue
    "Not detected" = "#D95F02",
    "No change"      = "grey70"
  )) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  theme_classic() +
  theme(
    axis.text  = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) +
  labs(
    x = "",
    y = "Log fold-change (Not detected / Detected)",
    fill = "Higher abundance in",
    title = "CN differential abundance (NB model)"
  )


plot_df$logFC_plot <- -plot_df$logFC

ggplot(plot_df,
       aes(x = CN, y = logFC_plot, fill = Direction)) +
  
  geom_col(width = 0.9) +
  
  geom_text(
    aes(y = ifelse(logFC_plot > 0,
                   logFC_plot + 0.05,
                   logFC_plot - 0.05),
        label = signif_label),
    size = 5,
    fontface = "bold",
    hjust = ifelse(plot_df$logFC_plot > 0, 0, 1)
  ) +
  
  scale_fill_manual(values = c(
    "Detected" = "#2C7BB6",   # calm blue
    "Not detected" = "#D95F02",
    "No change"      = "grey70"
  )) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  theme_classic() +
  theme(
    axis.text  = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) +
  labs(
    x = "",
    y = "Log fold-change (Detected / Not detected)",
    fill = "Higher abundance in",
    title = "CN differential abundance (NB model)"
  )


plot_df2 <- plot_df %>%
  arrange(logFC_plot) %>%
  mutate(CN = factor(CN, levels = CN))


ggplot(plot_df2,
       aes(x = CN, y = logFC_plot, fill = Direction)) +
  
  geom_col(width = 0.9) +
  
  geom_text(
    aes(y = ifelse(logFC_plot > 0,
                   logFC_plot + 0.05,
                   logFC_plot - 0.05),
        label = signif_label),
    size = 5,
    fontface = "bold",
    hjust = ifelse(plot_df$logFC_plot > 0, 0, 1)
  ) +
  
  scale_fill_manual(values = c(
    "Detected" = "#2C7BB6",   # calm blue
    "Not detected" = "#D95F02",
    "No change"      = "grey70"
  )) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  theme_classic() +
  theme(
    axis.text  = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) +
  labs(
    x = "",
    y = "-log(FC)",
    fill = "Higher abundance in",
    title = "CN differential abundance by ctdna"
  )


nb_results <- lapply(split(da_long, da_long$CN), function(df){
  
  fit <- glm.nb(Density ~ ctdna, data = df)
  
  coef <- summary(fit)$coefficients
  ci <- suppressMessages(confint.default(fit))
  
  data.frame(
    CN = unique(df$CN),
    logFC = -coef["ctdnaNot detected","Estimate"],
    lower = -ci["ctdnaNot detected",1],
    upper = -ci["ctdnaNot detected",2],
    p_value = coef["ctdnaNot detected","Pr(>|z|)"]
  )
  
}) %>% bind_rows()

plot_df <- nb_results %>%
  mutate(
    Direction = case_when(
      logFC > 0 ~ "Detected",
      logFC < 0 ~ "Not detected",
      TRUE ~ "No change"
    ),
    FDR = p.adjust(p_value, method = "BH"),
    neglogFDR = -log10(FDR),
    stars = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ ""
    )
  ) %>%
  arrange(logFC) %>%
  mutate(
    CN = factor(CN, levels = CN),
    star_x = ifelse(logFC > 0, upper + 0.08, lower - 0.08)
  )

ggplot(plot_df, aes(logFC, CN)) +
  
  geom_vline(xintercept = 0, linetype = 2,
             colour = "grey70", linewidth = 0.6) +
  
  geom_errorbarh(aes(xmin = lower, xmax = upper),
                 height = 0.15,
                 colour = "grey60",
                 linewidth = 0.9) +
  
  geom_segment(aes(x = 0, xend = logFC,
                   y = CN, yend = CN),
               colour = "grey88",
               linewidth = 0.8) +
  
  geom_point(aes(colour = Direction,
                 size = neglogFDR),
             alpha = 0.95) +
  
  geom_text(aes(x = star_x,
                label = stars),
            fontface = "bold",
            size = 5) +
  
  scale_colour_manual(values = c(
    "Detected" = "#2C7BB6",   # calm blue
    "Not detected" = "#D95F02",
    "No change"      = "grey70"
  )) +
  
  scale_size(range = c(3, 9),
             name = expression(-log[10](p_value))) +
  
  labs(
    title = "Differential cell abundance by ctDNA status",
    x = expression(log[2]*" fold-change (Detected vs Not detected)"),
    y = NULL,
    colour = "Higher abundance in"
  ) +
  
  theme_classic(base_size = 16) +
  
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(face = "bold", colour = "black"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold")
  )


