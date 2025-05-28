install.packages("patchwork")
install.packages("gridExtra")
library(tidyverse)
library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyr)
library(data.table)
library(readxl)
library(gridExtra)

# Spine-heart heatmap
# Read CSV file with proper encoding for Chinese
xinzang <- read.csv("/code/Barplot/New_Cardiac_Phenotype_Results_Standardized_Coefficients_FDR_Adjusted.csv", 
                    fileEncoding = "GBK")

# Check data
head(xinzang)

# Clean "...Instance.2" from outcome column if it exists
if ("outcome" %in% colnames(xinzang)) {
  xinzang$outcome <- gsub("\\.\\.\\.Instance\\.2$", "", xinzang$outcome)
  print("Cleaned outcome column")
} else {
  warning("No 'outcome' column found in the dataframe!")
}

# View results
head(xinzang$outcome)

# Create new dataframe retu_xinzang2 with renamed disease names
retu_xinzang <- xinzang %>%
  mutate(
    outcome = fct_recode(
      factor(outcome),  # Ensure outcome is a factor
      # ========== Name mapping ==========
      "Asc.aorta.D" = "Ascending.aorta.distensibility",
      "Asc.aorta.max" = "Ascending.aorta.maximum.area",
      "Asc.aorta.min" = "Ascending.aorta.minimum.area",
      "Desc.aorta.D" = "Descending.aorta.distensibility",
      "Desc.aorta.max" = "Descending.aorta.maximum.area",
      "Desc.aorta.min" = "Descending.aorta.minimum.area",
      "LAEF" = "LA.ejection.fraction",
      "LA.vol.max" = "LA.maximum.volume",
      "LA.vol.min" = "LA.minimum.volume",
      "LA.SV" = "LA.stroke.volume",
      "LV.CO" = "LV.cardiac.output",
      "LV.CS.A1" = "LV.circumferential.strain.AHA.1",
      "LV.CS.A2" = "LV.circumferential.strain.AHA.2",
      "LV.CS.A3" = "LV.circumferential.strain.AHA.3",
      "LV.CS.A4" = "LV.circumferential.strain.AHA.4",
      "LV.CS.A5" = "LV.circumferential.strain.AHA.5",
      "LV.CS.A6" = "LV.circumferential.strain.AHA.6",
      "LV.CS.A7" = "LV.circumferential.strain.AHA.7",
      "LV.CS.A8" = "LV.circumferential.strain.AHA.8",
      "LV.CS.A9" = "LV.circumferential.strain.AHA.9",
      "LV.CS.A10" = "LV.circumferential.strain.AHA.10",
      "LV.CS.A11" = "LV.circumferential.strain.AHA.11",
      "LV.CS.A12" = "LV.circumferential.strain.AHA.12",
      "LV.CS.A13" = "LV.circumferential.strain.AHA.13",
      "LV.CS.A14" = "LV.circumferential.strain.AHA.14",
      "LV.CS.A15" = "LV.circumferential.strain.AHA.15",
      "LV.CS.A16" = "LV.circumferential.strain.AHA.16",
      "LV.CS.global" = "LV.circumferential.strain.global",
      "LVEF" = "LV.ejection.fraction",
      "LV.EDV" = "LV.end.diastolic.volume",
      "LV.ESV" = "LV.end.systolic.volume",
      "LV.LS.S1" = "LV.longitudinal.strain.Segment.1",
      "LV.LS.S2" = "LV.longitudinal.strain.Segment.2",
      "LV.LS.S3" = "LV.longitudinal.strain.Segment.3",
      "LV.LS.S4" = "LV.longitudinal.strain.Segment.4",
      "LV.LS.S5" = "LV.longitudinal.strain.Segment.5",
      "LV.LS.S6" = "LV.longitudinal.strain.Segment.6",
      "LV.LS.global" = "LV.longitudinal.strain.global",
      "LV.WT.A1" = "LV.mean.myocardial.wall.thickness.AHA.1",
      "LV.WT.A2" = "LV.mean.myocardial.wall.thickness.AHA.2",
      "LV.WT.A3" = "LV.mean.myocardial.wall.thickness.AHA.3",
      "LV.WT.A4" = "LV.mean.myocardial.wall.thickness.AHA.4",
      "LV.WT.A5" = "LV.mean.myocardial.wall.thickness.AHA.5",
      "LV.WT.A6" = "LV.mean.myocardial.wall.thickness.AHA.6",
      "LV.WT.A7" = "LV.mean.myocardial.wall.thickness.AHA.7",
      "LV.WT.A8" = "LV.mean.myocardial.wall.thickness.AHA.8",
      "LV.WT.A9" = "LV.mean.myocardial.wall.thickness.AHA.9",
      "LV.WT.A10" = "LV.mean.myocardial.wall.thickness.AHA.10",
      "LV.WT.A11" = "LV.mean.myocardial.wall.thickness.AHA.11",
      "LV.WT.A12" = "LV.mean.myocardial.wall.thickness.AHA.12",
      "LV.WT.A13" = "LV.mean.myocardial.wall.thickness.AHA.13",
      "LV.WT.A14" = "LV.mean.myocardial.wall.thickness.AHA.14",
      "LV.WT.A15" = "LV.mean.myocardial.wall.thickness.AHA.15",
      "LV.WT.A16" = "LV.mean.myocardial.wall.thickness.AHA.16",
      "LV.WT.global" = "LV.mean.myocardial.wall.thickness.global",
      "LV.Mass" = "LV.myocardial.mass",
      "LV.RS.A1" = "LV.radial.strain.AHA.1",
      "LV.RS.A2" = "LV.radial.strain.AHA.2",
      "LV.RS.A3" = "LV.radial.strain.AHA.3",
      "LV.RS.A4" = "LV.radial.strain.AHA.4",
      "LV.RS.A5" = "LV.radial.strain.AHA.5",
      "LV.RS.A6" = "LV.radial.strain.AHA.6",
      "LV.RS.A7" = "LV.radial.strain.AHA.7",
      "LV.RS.A8" = "LV.radial.strain.AHA.8",
      "LV.RS.A9" = "LV.radial.strain.AHA.9",
      "LV.RS.A10" = "LV.radial.strain.AHA.10",
      "LV.RS.A11" = "LV.radial.strain.AHA.11",
      "LV.RS.A12" = "LV.radial.strain.AHA.12",
      "LV.RS.A13" = "LV.radial.strain.AHA.13",
      "LV.RS.A14" = "LV.radial.strain.AHA.14",
      "LV.RS.A15" = "LV.radial.strain.AHA.15",
      "LV.RS.A16" = "LV.radial.strain.AHA.16",
      "LV.RS.global" = "LV.radial.strain.global",
      "LV.SV" = "LV.stroke.volume",
      "RAEF" = "RA.ejection.fraction",
      "RA.vol.max" = "RA.maximum.volume",
      "RA.vol.min" = "RA.minimum.volume",
      "RA.SV" = "RA.stroke.volume",
      "RVEF" = "RV.ejection.fraction",
      "RV.EDV" = "RV.end.diastolic.volume",
      "RV.ESV" = "RV.end.systolic.volume",
      "RV.SV" = "RV.stroke.volume"
    )
  )

# Check result
head(retu_xinzang$disease)
levels(retu_xinzang$disease)

# -------------------- Data preprocessing --------------------
# Sort by Category_self
plot_data <- retu_xinzang %>%
  mutate(
    standardized_estimate = as.numeric(standardized_estimate),
    term = factor(term, levels = unique(term)),  # Preserve original order of term
    Category_self = factor(
      Category_self,
      levels = sort(unique(Category_self), decreasing = TRUE)
    )
  ) %>%
  arrange(Category_self, desc(p.value)) %>%
  mutate(outcome = factor(outcome, levels = unique(outcome))) %>%
  group_by(outcome) %>%
  mutate(y_pos = cur_group_id()) %>%
  ungroup()

# -------------------- Create background strip --------------------
term_count <- length(unique(plot_data$term))  # Dynamically count term groups
background_data <- plot_data %>%
  group_by(Category_self) %>%
  summarise(
    ymin = min(y_pos) - 0.5,
    ymax = max(y_pos) + 0.5,
    .groups = "drop"
  ) %>%
  mutate(
    xmin = 0.5,
    xmax = term_count + 0.5
  )

# -------------------- Heatmap plot --------------------
p <- ggplot() +
  geom_rect(
    data = background_data,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "grey", alpha = 0.15
  ) +
  geom_point(
    data = plot_data,
    aes(
      x = term,
      y = y_pos,
      size = abs(standardized_estimate),
      color = standardized_estimate
    ),
    shape = 19, alpha = 0.8
  ) +
  scale_color_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#B2182B",
    midpoint = 0,
    name = "Standardized\nEstimate",
    guide = guide_colorbar(order = 1)
  ) +
  scale_size_continuous(
    range = c(6, 14),
    name = "|Standardized\nEstimate|",
    guide = guide_legend(order = 2)
  ) +
  scale_x_discrete(
    name = "Term Groups",
    expand = expansion(add = 0.2)
  ) +
  scale_y_continuous(
    breaks = unique(plot_data$y_pos),
    labels = levels(plot_data$outcome),
    expand = expansion(add = 0.5)
  ) +
  labs(title = "Association Analysis by Term Groups") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.spacing.y = unit(0.5, "cm"),
    legend.key.size = unit(1.2, "cm")
  )

# Export PDF
ggsave("/results/Heart_Spine_Heatmap.pdf", 
       plot = p,
       width = 8,
       height = length(unique(plot_data$outcome)) * 0.4 + 5,
       device = "pdf")

# -------------------- FDR bar plot --------------------
# Create term color mapping
term_colors <- c(
  "M45_status" = "#073E7F",
  "M46_status" = "#83A0BE", 
  "M47_status" = "#866AA3",
  "M48_status" = "#C1B7CF",
  "M50_status" = "#BE0E23",
  "M51_status" = "#E18791"
)

plot_data <- plot_data %>%
  arrange(Category_self, outcome) %>%
  mutate(
    outcome = factor(outcome, levels = unique(outcome)),
    term = factor(term, levels = names(term_colors))
  ) %>%
  group_by(outcome) %>%
  mutate(
    log10_p = -log10(FDR),
    sig_label = ifelse(FDR < 0.05, "*", ""),
    cumulative_p = cumsum(log10_p),
    label_y = cumulative_p - log10_p/2
  ) %>%
  ungroup()

# Plotting FDR bar chart
p <- ggplot(plot_data, aes(x = outcome, y = log10_p, fill = term)) +
  geom_col(
    width = 0.6,
    position = position_stack(reverse = TRUE)
  ) +
  scale_fill_manual(
    values = term_colors,
    breaks = names(term_colors),
    name = "Term Groups"
  ) +
  geom_text(
    aes(y = label_y, label = sig_label),
    color = "white", 
    size = 5,
    fontface = "bold",
    hjust = 0.5
  ) +
  scale_x_discrete(limits = levels(plot_data$outcome)) +
  coord_flip() +
  labs(
    title = "Association FDR by Term Groups",
    y = "-log10(FDR)",
    x = "Outcome"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10, margin = margin(r = 15)),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    panel.grid.major.y = element_blank(),
    plot.title = element_text(hjust = 0.5, margin = margin(b = 15))
  )

# Dynamically calculate height
chart_height <- length(unique(plot_data$outcome)) * 0.4 + 3

ggsave("/results/Heart_Spine_FDR_Barplot.pdf", p,
       width = 10, 
       height = chart_height,
       limitsize = FALSE)
