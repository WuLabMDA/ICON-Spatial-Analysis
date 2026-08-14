rm(list = ls(all = TRUE))
graphics.off()

library(SpatialExperiment)
library(ggplot2)
library(dplyr)
library(openxlsx)
library(reshape2)
library(tibble)
library(forcats)
library(ggnewscale)
library(patchwork)

tempdata <- readRDS("Y:/Projects/ICON IMC/Processed data (rescale)/speCelltypes2.rds")
data <- tempdata
data <- data[, data$celltypes != "Undefined"]

clinData <- readWorkbook(
  "Y:/Projects/ICON IMC/Clinical data/sortedClinical_LindaRevision.xlsx",
  sheet = "Sheet 1",
  colNames = TRUE)
rownames(clinData) <- clinData$SampleID


tempMetadata <- readWorkbook(
  "Y:/Projects/ICON IMC/Clinical data/roiMetadata2.xlsx",
  sheet = "Sheet 1",
  colNames = TRUE)
tempMetadata <- tempMetadata[,-1]

roiMetadata <- tempMetadata[tempMetadata$ROI.NO %in% unique(data$sample_id),]
rownames(roiMetadata) <- roiMetadata$ROI.NO

unq <- unique(roiMetadata$ROI.NO)

celltypes <- as.character(unique(data@colData@listData[["celltypes"]]))

freqTable <- as.data.frame(matrix(0, nrow = length(unq), ncol = length(celltypes)))
colnames(freqTable) <- celltypes
rownames(freqTable) <- roiID <- unq

propTable <- as.data.frame(matrix(0, nrow = length(unq), ncol = length(celltypes)))
colnames(propTable) <- celltypes
rownames(propTable) <- roiID

densityTable <- as.data.frame(matrix(0, nrow = length(unq), ncol = length(celltypes)))
colnames(densityTable) <- celltypes
rownames(densityTable) <- roiID

for (ii in 1:length(roiID)){
  message("Processing: ",ii," of ",length(roiID))
  tidx <- roiMetadata[roiMetadata$ROI.NO==roiID[ii],]
  fdata <- matrix(0, nrow = 1, ncol = length(celltypes))
  colnames(fdata) <- celltypes
  tdata <- data[,grepl(unq[ii], data@colData@listData[["sample_id"]])]
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


for (i in 1:length(unq)){
  densityTable[unq[i], "Location"] <- roiMetadata[unq[i],"Location"]
  freqTable[unq[i], "Location"] <- roiMetadata[unq[i],"Location"]
  propTable[unq[i], "Location"] <- roiMetadata[unq[i],"Location"]
}

densityTable <- na.omit(densityTable)
freqTable <- na.omit(freqTable)
propTable <- na.omit(propTable)

tpropTable <- propTable


tprop_long <- tpropTable %>%
  rownames_to_column(var = "Sample") %>%
  melt(id.vars = c("Sample", "Location"), variable.name = "CellType", value.name = "Proportion")
colnames(tprop_long) <- c("Sample", "Location","CellType","Proportion")

epithelial_order <- tprop_long %>%
  filter(CellType == "Epithelial cells") %>%
  arrange(Proportion)  

tprop_long$Sample <- factor(tprop_long$Sample, levels = epithelial_order$Sample)

celltype_colors <- data@metadata[["color_vectors"]][["celltype"]]

tprop_long$CellType <- factor(tprop_long$CellType, levels = rev(names(celltype_colors)))


epi_rank <- tprop_long %>% 
  filter(CellType == "Epithelial cells") %>% 
  group_by(Location) %>%          # note the extra Location level
  arrange(Proportion, .by_group = TRUE) %>% 
  mutate(rank = row_number()) %>% 
  ungroup() %>% 
  dplyr::select(Sample, Location, rank)

tprop_long <- tprop_long %>% 
  left_join(epi_rank, by = c("Sample", "Location")) %>% 
  arrange(Location, rank) %>% 
  mutate(Sample = factor(Sample, levels = unique(Sample)))

tprop_long <- tprop_long %>%
  mutate(
    CellType = forcats::fct_rev(
      forcats::fct_relevel(CellType, "Epithelial cells")
    )
  )




loc_cols <- c("Tumor" = "firebrick", "Normal" = "steelblue")

make_panel <- function(df, location_label) {
  
  sample_strip <- distinct(df[, c("Sample", "Location")])
  
  ggplot() +
    
    geom_tile(data = sample_strip,
              aes(x = Sample, y = 1.05, fill = Location),
              height = 0.05) +
    scale_fill_manual(values = loc_cols, name = NULL) +   # name=NULL => no legend title
    ggnewscale::new_scale_fill() +
    

    geom_bar(data = df,
             aes(x = Sample, y = Proportion, fill = CellType),
             stat = "identity", position = "fill") +
    scale_fill_manual(values = celltype_colors, name = "Cell Type") +
    
    coord_cartesian(ylim = c(0, 1.1)) +
    
    labs(title = paste(location_label, "samples"),
         x = NULL, y = "Proportion") +
    
    theme_minimal() +
    theme(
      strip.text.x  = element_text(face = "bold"),
      axis.text.x   = element_blank(),
      axis.ticks.x  = element_blank(),
      legend.position = "right"
    )
}

plot_pre <- make_panel(filter(tprop_long, Location == "Normal"), "Normal")
plot_post     <- make_panel(filter(tprop_long, Location == "Tumor"), "Tumor")

(plot_pre / plot_post) + plot_layout(guides = "collect")






