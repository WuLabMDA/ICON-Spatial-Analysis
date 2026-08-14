rm(list = ls(all = TRUE))
graphics.off()

library(ggplot2)
library(dplyr)
library(imcRtools)
library(rstatix)
library(ggforce)
library(ggnewscale)
library(scales)


load("Y:/Projects/ICON IMC/Final/S1 landscape/Data/CCI_workspace_TN.rds")

rec_data <- tempdata[,tempdata@colData@listData[["recurrence"]] %in% c("YES","NO")]
rec_y <- rec_data[,rec_data@colData@listData[["recurrence"]] == "YES"]
rec_y <- unique(rec_y@colData@listData[["sample_id"]])

rec_n <- rec_data[,rec_data@colData@listData[["recurrence"]] == "NO"]
rec_n <- unique(rec_n@colData@listData[["sample_id"]])

pre_out_recy <- pre_out[pre_out@listData[["group_by"]] %in% rec_y,]
pre_out_recn <- pre_out[pre_out@listData[["group_by"]] %in% rec_n,]


pre_recy_df <- pre_out_recy %>%
  as_tibble() %>%
  group_by(from_label, to_label) %>%
  summarize(
    mean_sigval = mean(sigval, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  mutate(context = "pre_recy")

pre_recn_df <- pre_out_recn %>%
  as_tibble() %>%
  group_by(from_label, to_label) %>%
  summarize(
    mean_sigval = mean(sigval, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  mutate(context = "pre_recn")

comb_df <- bind_rows(pre_recy_df, pre_recn_df)

comb_df <- comb_df %>%
  mutate(
    x = as.numeric(factor(from_label)),
    y = as.numeric(factor(to_label)),
    start = ifelse(context == "pre_recy",  pi, 0),
    end   = ifelse(context == "pre_recy", 2*pi, pi)
  )

global_lim <- quantile(abs(comb_df$mean_sigval), 0.98, na.rm = TRUE)

comb_df <- comb_df %>%
  mutate(
    r = scales::rescale(
      sqrt(abs(mean_sigval)),
      to = c(0.08, 0.45)
    )
  )



stat_df <- bind_rows(
  pre_out_recy  %>% as_tibble() %>% mutate(context = "pre_recy"),
  pre_out_recn %>% as_tibble() %>% mutate(context = "pre_recn")
) %>%
  group_by(from_label, to_label) %>%
  wilcox_test(sigval ~ context) %>%
  adjust_pvalue(method = "fdr") %>%
  mutate(
    signif = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE          ~ ""
    )
  ) %>%
  ungroup()

stat_df <- stat_df %>%
  mutate(
    x = as.numeric(factor(from_label, levels = levels(factor(comb_df$from_label)))),
    y = as.numeric(factor(to_label,   levels = levels(factor(comb_df$to_label)))),
    y_star = y + 0.35
  )



ggplot() +
  
  ## ---- PRE ----
geom_arc_bar(
  data = comb_df %>% filter(context == "pre_recy"),
  aes(
    x0 = x, y0 = y,
    r0 = 0,
    r = r,
    start = start,
    end = end,
    fill = mean_sigval
  ),
  color = "grey30",
  linewidth = 0.1
) +
  scale_fill_gradient2(
    low = "#3B4CC0",
    mid = "white",
    high = "#B40426",
    midpoint = 0,
    limits = c(-global_lim, global_lim),
    oob = scales::squish,
    name = "Interaction strength\n(mean per sample)"
  ) +
  
  ggnewscale::new_scale_fill() +
  
  ## ---- POST ----
geom_arc_bar(
  data = comb_df %>% filter(context == "pre_recn"),
  aes(
    x0 = x, y0 = y,
    r0 = 0,
    r = r,
    start = start,
    end = end,
    fill = mean_sigval
  ),
  color = "grey30",
  linewidth = 0.1
) +
  scale_fill_gradient2(
    low = "#3B4CC0",
    mid = "white",
    high = "#B40426",
    midpoint = 0,
    limits = c(-global_lim, global_lim),
    oob = scales::squish,
    name = "Interaction strength\n(mean per sample)"
  ) +
  
  ## ---- STARS ----
geom_text(
  data = stat_df %>% filter(signif != ""),
  aes(x = x, y = y_star, label = signif),
  inherit.aes = FALSE,
  size = 5,
  fontface = "bold",
  color = "black"
) +
  
  ## ---- AXES ----
scale_x_continuous(
  breaks = seq_along(levels(factor(comb_df$from_label))),
  labels = levels(factor(comb_df$from_label))
) +
  scale_y_continuous(
    breaks = seq_along(levels(factor(comb_df$to_label))),
    labels = levels(factor(comb_df$to_label))
  ) +
  coord_fixed() +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12)
  ) +
  labs(
    x = "From cell type",
    y = "To cell type",
    title = "Cell–Cell Interaction Strength Normal (Recurrence)",
    subtitle = "Left half = Recurred, Right half = Not recurred\nBubble size = |mean interaction strength per sample|"
  )




post_out_recy <- post_out[post_out@listData[["group_by"]] %in% rec_y,]
post_out_recn <- post_out[post_out@listData[["group_by"]] %in% rec_n,]


post_recy_df <- post_out_recy %>%
  as_tibble() %>%
  group_by(from_label, to_label) %>%
  summarize(
    mean_sigval = mean(sigval, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  mutate(context = "post_recy")

post_recn_df <- post_out_recn %>%
  as_tibble() %>%
  group_by(from_label, to_label) %>%
  summarize(
    mean_sigval = mean(sigval, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  mutate(context = "post_recn")

comb_df <- bind_rows(post_recy_df, post_recn_df)

comb_df <- comb_df %>%
  mutate(
    x = as.numeric(factor(from_label)),
    y = as.numeric(factor(to_label)),
    start = ifelse(context == "post_recy",  pi, 0),
    end   = ifelse(context == "post_recy", 2*pi, pi)
  )

global_lim <- quantile(abs(comb_df$mean_sigval), 0.98, na.rm = TRUE)

comb_df <- comb_df %>%
  mutate(
    r = scales::rescale(
      sqrt(abs(mean_sigval)),
      to = c(0.08, 0.45)
    )
  )



stat_df <- bind_rows(
  post_out_recy  %>% as_tibble() %>% mutate(context = "post_recy"),
  post_out_recn %>% as_tibble() %>% mutate(context = "post_recn")
) %>%
  group_by(from_label, to_label) %>%
  wilcox_test(sigval ~ context) %>%
  adjust_pvalue(method = "fdr") %>%
  mutate(
    signif = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE          ~ ""
    )
  ) %>%
  ungroup()

stat_df <- stat_df %>%
  mutate(
    x = as.numeric(factor(from_label, levels = levels(factor(comb_df$from_label)))),
    y = as.numeric(factor(to_label,   levels = levels(factor(comb_df$to_label)))),
    y_star = y + 0.35
  )



ggplot() +
  
  ## ---- PRE ----
geom_arc_bar(
  data = comb_df %>% filter(context == "post_recy"),
  aes(
    x0 = x, y0 = y,
    r0 = 0,
    r = r,
    start = start,
    end = end,
    fill = mean_sigval
  ),
  color = "grey30",
  linewidth = 0.1
) +
  scale_fill_gradient2(
    low = "#3B4CC0",
    mid = "white",
    high = "#B40426",
    midpoint = 0,
    limits = c(-global_lim, global_lim),
    oob = scales::squish,
    name = "Interaction strength\n(mean per sample)"
  ) +
  
  ggnewscale::new_scale_fill() +
  
  ## ---- POST ----
geom_arc_bar(
  data = comb_df %>% filter(context == "post_recn"),
  aes(
    x0 = x, y0 = y,
    r0 = 0,
    r = r,
    start = start,
    end = end,
    fill = mean_sigval
  ),
  color = "grey30",
  linewidth = 0.1
) +
  scale_fill_gradient2(
    low = "#3B4CC0",
    mid = "white",
    high = "#B40426",
    midpoint = 0,
    limits = c(-global_lim, global_lim),
    oob = scales::squish,
    name = "Interaction strength\n(mean per sample)"
  ) +
  
  ## ---- STARS ----
geom_text(
  data = stat_df %>% filter(signif != ""),
  aes(x = x, y = y_star, label = signif),
  inherit.aes = FALSE,
  size = 5,
  fontface = "bold",
  color = "black"
) +
  
  ## ---- AXES ----
scale_x_continuous(
  breaks = seq_along(levels(factor(comb_df$from_label))),
  labels = levels(factor(comb_df$from_label))
) +
  scale_y_continuous(
    breaks = seq_along(levels(factor(comb_df$to_label))),
    labels = levels(factor(comb_df$to_label))
  ) +
  coord_fixed() +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12)
  ) +
  labs(
    x = "From cell type",
    y = "To cell type",
    title = "Cell–Cell Interaction Strength Tumor (Recurrence)",
    subtitle = "Left half = Recurred, Right half = Not recurred\nBubble size = |mean interaction strength per sample|"
  )


