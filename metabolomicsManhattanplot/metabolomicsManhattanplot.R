# Clear workspace
rm(list = ls())
install.packages("ggtext")
install.packages("ggrepel")
# Load required libraries
library(ggrepel)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggtext)
library(readxl)

# Read data from Excel file
extreme_cox_padj_bubble_plot1 <- read_xlsx("/code/metabolomicsManhattanplot/Manhattanplotsixspinaldisorders.xlsx")

# Convert Category to factor and set order
extreme_cox_padj_bubble_plot1 <- extreme_cox_padj_bubble_plot1 %>%
  mutate(Category = factor(Category, levels = c("Bone and joint", "Renal/Liver function",
                                                "Inflammation","Endocrine",
                                                "Amino acids","Glycolysis related metabolites", 
                                                "Ketone bodies","Fatty acids",
                                                "Triglycerides","Free Cholesterol",
                                                "Cholesteryl esters","Cholesterol",
                                                "Phospholipids", "Other lipids",
                                                "Total lipids",
                                                "Lipoprotein particle concentrations", 
                                                "Lipoprotein particle sizes",
                                                "Lipoproteins and Related Proteins", 
                                                "Red blood cell","Platelet","White blood cell"))) 

# Calculate chromosome size and cumulative position
extreme_cox_padj_bubble_plot1 <- extreme_cox_padj_bubble_plot1 %>%
  group_by(Category) %>%
  summarise(chr_len = max(BP)) %>%
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%
  left_join(extreme_cox_padj_bubble_plot1, ., by = "Category") %>%
  arrange(Category, BP) %>%
  mutate(BPcum = BP + tot)

# Calculate chromosome center positions
axisdf <- extreme_cox_padj_bubble_plot1 %>%
  group_by(Category) %>%
  summarize(center = (max(BPcum) + min(BPcum)) / 2)

# Create color variable based on logFC
extreme_cox_padj_bubble_plot1$color_var <- ifelse(
  extreme_cox_padj_bubble_plot1$logFC > 0, "#E53935",  # Dark red if logFC > 0
  "#7094d0"   # Dark blue otherwise
)

# Set point shapes based on disease type
extreme_cox_padj_bubble_plot1 <- extreme_cox_padj_bubble_plot1 %>%
  mutate(shape_var = case_when(
    disease == "M45_status" ~ 8,  # shape 8
    disease == "M46_status" ~ 9,  # shape 9
    disease == "M47_status" ~ 10, # shape 10
    disease == "M48_status" ~ 19, # shape 19
    disease == "M50_status" ~ 17, # shape 17
    disease == "M51_status" ~ 5,  # shape 5
    TRUE ~ 16                   # Default shape (filled circle)
  ))
extreme_cox_padj_bubble_plot1$shape_var <- as.factor(extreme_cox_padj_bubble_plot1$shape_var)

# Filter for significant metabolites (adjusted_p < 0.05)
extreme_cox_padj_bubble_plot <- extreme_cox_padj_bubble_plot1 %>%
  group_by(Category) %>%
  filter(adjusted_p < 0.05) %>%
  mutate(
    n_in_category = n(),  # Count of significant metabolites in this category
    rank_size_logFC = if_else(n_in_category >= 1, rank(-size_logFC), NA_real_),
    anno = case_when(
      rank_size_logFC == 1 ~ Abbreviation,  # Only show label for top-ranked metabolite
      TRUE ~ NA_character_  # No label for others
    )
  ) %>%
  ungroup()

# Create ggplot object
plot_p_value_filtered <- ggplot() +
  
  # Plot significant points (adjusted_p < 0.05)
  geom_point(
    data = filter(extreme_cox_padj_bubble_plot, adjusted_p < 0.05),
    aes(x = BPcum, y = logFC, color = color_var, size = size_logFC1, shape = shape_var),
    position = position_jitter(width = 0.5)
  ) +
  scale_size_continuous(range = c(0, 30)) +  # Set size range
  
  # Plot non-significant points (adjusted_p > 0.05) in gray
  geom_jitter(
    data = filter(extreme_cox_padj_bubble_plot1, adjusted_p > 0.05),
    aes(x = BPcum, y = logFC),
    color = "grey", size = 2, position = position_jitter(width = 0.5)
  ) +
  
  # Add labels for top metabolites
  geom_text_repel(
    data = filter(extreme_cox_padj_bubble_plot, adjusted_p < 0.05), 
    aes(x = BPcum, y = logFC, label = anno, color = color_var),
    size = 9, na.rm = TRUE,
    max.overlaps = 20,  # Increase maximum allowed overlaps
    min.segment.length = 0.15,  # Set minimum segment length
    nudge_x = 0.05,  # Adjust label position on x-axis
    nudge_y = ifelse(filter(extreme_cox_padj_bubble_plot, adjusted_p < 0.05)$logFC < 0, -0.08, 0.05)
  ) +
  
  # Set color scheme
  scale_color_manual(values = c("#E53935" = "#E53935", "#7094d0" = "#7094d0", "grey" = "grey")) + 
  
  # Set point shapes with legend
  scale_shape_manual(
    values = c(8, 9, 10, 19, 17, 5, 16),  # Shape values
    breaks = c("8", "9", "10", "19", "17", "5", "16"),  # Shapes to show in legend
    labels = c("M45_status", "M46_status", "M47_status", "M48_status", "M50_status", "M51_status", "Other")  # Legend labels
  ) +
  
  # Add horizontal reference line
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  
  # Configure x-axis
  scale_x_continuous(label = axisdf$Category, breaks = axisdf$center) +
  
  # Configure y-axis
  scale_y_continuous(
    limits = c(-0.2, 0.40),
    breaks = c(-0.2, -0.15, -0.1, -0.05, 0,
               0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4),
    expand = c(0, 0)
  ) +
  
  # Set axis labels
  labs(
    x = "",
    y = "logFC"
  ) +
  
  # Apply theme settings
  theme_classic() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    axis.text.x = element_text(angle = 20, hjust = 1, size = 16),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_markdown(size = 20)
  )

# Save plot as PDF
ggsave("/results/manhattan.pdf",   height = 1200,
  width = 3000,
  units = "px",
  device = "pdf", dpi = 72
)
