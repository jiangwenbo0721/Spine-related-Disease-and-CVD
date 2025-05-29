library(tidyverse)
library(cowplot)
library(readxl)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(patchwork)
library(dplyr)
library(tidyr)
library(readxl)
library(writexl)
getwd()
setwd("") 


# 1. Read and preprocess data
aa <- read_excel("brainMRI_data_FDR.xlsx", sheet = 1)
aa$outcome <- make.names(aa$outcome)

# Data Conversion (Retain Region Information)
transformed_data <- aa %>%
  select(term, standardized_estimate, FDR, outcome, Region) %>%
  pivot_wider(
    names_from = term,
    values_from = c(standardized_estimate, FDR),
    names_glue = "{term}_{.value}"
  ) %>%
  rename(Metric = outcome) %>%
  select(
    Metric,
    Region,
    M47_status_standardized_estimate,
    M48_status_standardized_estimate,
    M47_status_FDR,
    M48_status_FDR
  ) %>%
  rename(
    M47_status = M47_status_standardized_estimate,
    M48_status = M48_status_standardized_estimate,
    M47_FDR = M47_status_FDR,
    M48_FDR = M48_status_FDR
  )

aa <- transformed_data

# 2. Clean up the indicator names (remove the left and right marks)

aa <- aa %>%
  mutate(
    Metric = gsub("^IDP\\.T1\\.FAST\\.ROIs\\.", "", Metric),
    Metric = gsub("^Weighted\\.mean\\.", "", Metric),
    Metric = gsub("\\.orientation\\.dispersion\\.index\\.\\.in\\.tract\\.", "", Metric),
    Metric = gsub("\\.diffusion\\.tensor\\.mode\\.\\.in\\.tract\\.", "", Metric),
    Metric = gsub("\\.mean\\.diffusivity\\.\\.in\\.tract\\.", "", Metric),
    Metric = gsub("\\.isotropic\\.or\\.free\\.water\\.volume\\.fraction\\.\\.in\\.tract\\.", "", Metric),
    Metric = gsub("\\.intra\\.cellular\\.volume\\.fraction\\.\\.in\\.tract\\.", "", Metric),
    Metric = gsub("\\.fractional\\.anisotropy\\.\\.in\\.tract\\.", "", Metric),
    Metric = gsub("\\.\\.from\\.dMRI\\.data\\.$", "", Metric),
    # Remove the left and right marks (retain the names of the main parts)
    Metric_modified = str_remove(Metric, "^(FA|MD|MO|OD|ICVF|ISOVF)\\.(left|right)\\."),
    Metric_modified = str_remove(Metric_modified, "^(L|R)\\.")
  )

# 3. Add the Category category
aa <- aa %>%
  mutate(
    Category = case_when(
      grepl("^L|^R", Metric) ~ "GMV",
      grepl("^FA", Metric) ~ "FA",
      grepl("^MD", Metric) ~ "MD",
      grepl("^MO", Metric) ~ "MO",
      grepl("^OD", Metric) ~ "OD",
      grepl("^ICVF", Metric) ~ "ICVF",
      grepl("^ISOVF", Metric) ~ "ISOVF",
      TRUE ~ "Other"
    )
  )

# 4.-log10(FDR)
aa <- aa %>%
  mutate(
    M47_log_FDR = -log10(M47_FDR),
    M48_log_FDR = -log10(M48_FDR)
  )

# 5. Sorting (by Category and maximum FDR value, GMV is ranked first)
aa <- aa %>% 
  arrange(factor(Category, levels = c("GMV", "FA", "MD", "MO", "OD", "ICVF", "ISOVF")), 
          desc(pmax(M47_log_FDR, M48_log_FDR))) %>%
  mutate(Metric = factor(Metric, levels = unique(Metric)))

# 6. Create plotting data (Solve NA problems)
plot_data <- bind_rows(
  # M47L
  aa %>% filter(Region == "L") %>% 
    select(Metric, Metric_modified, Category, status = M47_status, FDR = M47_FDR) %>%
    mutate(model_region = "M47_L"),
  # M47R
  aa %>% filter(Region == "R") %>% 
    select(Metric, Metric_modified, Category, status = M47_status, FDR = M47_FDR) %>%
    mutate(model_region = "M47_R"),
  # M48L
  aa %>% filter(Region == "L") %>% 
    select(Metric, Metric_modified, Category, status = M48_status, FDR = M48_FDR) %>%
    mutate(model_region = "M48_L"),
  # M48R
  aa %>% filter(Region == "R") %>% 
    select(Metric, Metric_modified, Category, status = M48_status, FDR = M48_FDR) %>%
    mutate(model_region = "M48_R")
) %>%
  mutate(
    model_region = factor(model_region,
                          levels = c("M47_L", "M47_R", "M48_L", "M48_R"),
                          labels = c("M47 (L)", "M47 (R)", "M48 (L)", "M48 (R)"))
  )

# 7. Define color mapping
category_colors <- c(
  GMV = "#458B74",
  FA = "#4682B4",
  MD = "#1D1D8F",
  MO = "#DAA520",
  OD = "#CD5B45",
  ICVF = "#C71585",
  ISOVF = "#91D1C2",
  Other = "gray"
)

# 8. Create a bubble chart (no legend displayed, to be added uniformly later)
bubble_plot <- ggplot(plot_data, aes(x = model_region, y = Metric_modified)) +
  geom_point(aes(size = abs(status)), color = "black", shape = 1, stroke = 0.8) +
  geom_point(aes(size = abs(status), color = status), alpha = 0.8) +
  scale_size_continuous(range = c(3, 10), name = "Effect size") +
  scale_color_gradient2(
    low = "#08519C",
    mid = "white",
    high = "#A50F15",
    midpoint = 0,
    name = "Effect direction"
  ) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",  # Do not show the legend for now
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 0, face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 10),
    strip.text = element_text(face = "bold", size = 12),
    panel.spacing = unit(0.5, "lines"),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5)
  ) +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y")

# 9. Create an FDR bar chart (distinguishing left from right)

fdr_data <- plot_data %>%
  mutate(
    model = str_extract(model_region, "M47|M48"),
    region = str_extract(model_region, "\\(L\\)|\\(R\\)") %>% 
      str_remove_all("[()]"),
    log_FDR = -log10(FDR)
  ) %>%
  unite("model_region", model, region, sep = "_") %>%
  mutate(model_region = factor(model_region, 
                               levels = c("M47_L", "M47_R", "M48_L", "M48_R"),
                               labels = c("M47 (L)", "M47 (R)", "M48 (L)", "M48 (R)")))

fdr_plot <- ggplot(fdr_data, aes(y = Metric_modified, x = log_FDR, 
                                 fill = Category, alpha = model_region)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, 
           color = "black", linewidth = 0.3) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
  scale_fill_manual(values = category_colors, name = "Category") +
  scale_alpha_manual(
    values = c("M47 (L)" = 0.1, "M47 (R)" = 0.4, "M48 (L)" = 0.7, "M48 (R)" = 1),
    name = "Model & Region"
  ) +
  labs(x = "-log10(FDR)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(0.2, "cm"),
    legend.key.size = unit(0.8, "cm"),
    plot.margin = margin(5, 15, 5, 5),
    strip.text = element_blank()
  ) +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y")

# 10. Combined graphics
# Extract Legend
legend_bubble <- cowplot::get_legend(
  bubble_plot + 
    guides(size = guide_legend(order = 1), color = guide_legend(order = 2)) +
    theme(legend.position = "right",
          legend.box = "vertical",
          legend.spacing.y = unit(0.5, "cm"))
)

legend_fdr <- cowplot::get_legend(
  fdr_plot + 
    guides(fill = guide_legend(order = 1), alpha = guide_legend(order = 2)) +
    theme(legend.position = "right",
          legend.box = "vertical",
          legend.spacing.y = unit(0.5, "cm"))
)

# Combined main Image (without Legend)
main_plot <- bubble_plot + fdr_plot + 
  plot_layout(widths = c(4, 1.5)) & 
  theme(legend.position = "none")

# Add legend to the main image
final_plot <- main_plot +
  inset_element(legend_bubble, left = 0.9, bottom = 0.7, right = 1, top = 1) +
  inset_element(legend_fdr, left = 0.9, bottom = 0.4, right = 1, top = 0.65) +
  plot_annotation(theme = theme(plot.margin = margin(1, 5, 1, 1, "cm")))
# 1. Create a bubble chart (with legend)
bubble_plot_with_legend <- bubble_plot +
  theme(legend.position = "right",
        legend.box = "vertical",
        legend.spacing.y = unit(0.5, "cm"))

# 2. Create an FDR bar chart (with legend)
fdr_plot_with_legend <- fdr_plot +
  theme(legend.position = "right",
        legend.box = "vertical",
        legend.spacing.y = unit(0.5, "cm"))

# 3. Extract the legend
legend_bubble <- get_legend(bubble_plot_with_legend)
legend_fdr <- get_legend(fdr_plot_with_legend)

# 4. Create a main image without legends
bubble_plot_no_legend <- bubble_plot + theme(legend.position = "none")
fdr_plot_no_legend <- fdr_plot + theme(legend.position = "none")

# 5. Use plot_grid to combine graphs
combined_plot <- plot_grid(
  bubble_plot_no_legend, 
  fdr_plot_no_legend,
  nrow = 1,
  align = "h",
  axis = "tb",
  rel_widths = c(4, 1.5)
)

# 6. Add legends

final_plot <- plot_grid(
  combined_plot,
  plot_grid(legend_bubble, legend_fdr, ncol = 1),
  nrow = 1,
  rel_widths = c(8, 2)
)



# 7. Save the graphics
ggsave("brian_MRI.pdf",
       final_plot, 
       width = 8,
       height = 4 + length(unique(plot_data$Metric))*0.15,
       device = cairo_pdf,
       limitsize = FALSE)

# Display Graphics
print(final_plot)



