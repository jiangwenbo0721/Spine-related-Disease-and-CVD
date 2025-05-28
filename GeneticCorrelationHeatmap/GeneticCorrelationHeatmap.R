install.packages("ggnewscale")
library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)
library(ggnewscale)
library(readxl)
# Read data from Excel
data <- read_excel("/code/GeneticCorrelationHeatmap/table.xlsx")

# Create custom color scale function for LDSC
ldsc_fill <- scale_fill_gradientn(
  colours = c("lightblue", "#7094d0"),
  limits = range(data$value)
)

# Create custom color scale function for HDL
hdl_fill <- scale_fill_gradientn(
  colours = c("white", "#E53935"),
  limits = range(data$value)
)

# Modified plot code (using circles and removing grid lines)
ggplot() +
  # Remove background strips (by removing geom_rect part)
  
  # Add LDSC points (using new fill scale)
  ggnewscale::new_scale_fill() + 
  geom_point(data = subset(data, method == "LDSC"),
             aes(x = value, y = Category, fill = value),
             shape = 21, size = 9, color = "black", stroke = 0.5) + 
  scale_fill_gradientn(
    colours = c("lightblue", "#7094d0"),
    limits = range(data$value)
  ) +
  
  # Add LDSC p-value labels
  geom_text(data = subset(data, method == "LDSC"),
            aes(x = value, y = Category, label = sprintf("%.2g", p_value)),
            hjust = -0.1,
            size = 3,
            color = "black") +
  
  # Add HDL points (using another new fill scale)
  ggnewscale::new_scale_fill() +
  geom_point(data = subset(data, method == "HDL"),
             aes(x = value, y = Category, fill = value),
             shape = 21, size = 9, color = "black", stroke = 0.5) +
  scale_fill_gradientn(
    colours = c("white", "#E53935"),
    limits = range(data$value)
  ) +
  
  # Add HDL p-value labels
  geom_text(data = subset(data, method == "HDL"),
            aes(x = value, y = Category, label = sprintf("%.2g", p_value)),
            hjust = -0.1,
            size = 3,
            color = "black") +
  
  # Add plot title and axis labels
  labs(title = "LDSC and HDL Values by Category1",
       x = NULL,  # Remove x-axis label
       y = "Category") +
  
  # Apply minimal theme
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.title.x = element_blank(),  # Hide x-axis title
    axis.text.x = element_blank(),   # Hide x-axis ticks
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.5, "lines"),
    plot.margin = margin(5, 5, 5, 5, "mm"),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black")  # Add black border lines for the plot
  ) +
  
  # Facet the plot by method with separate columns
  facet_wrap(~method, ncol = 2, scales = "free_y") +  # Display each method in a separate column
  scale_y_discrete(expand = expansion(add = 0.2)) +  # Add space to y-axis to avoid overlap
  coord_cartesian(clip = "off")  # Prevent clipping of background strips
ggsave("/results/genetic correlation.png", width = 6, height = 4, dpi = 300)
