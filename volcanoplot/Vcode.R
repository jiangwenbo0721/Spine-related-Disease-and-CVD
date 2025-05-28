# Load required R packages
install.packages("ggrepel")
library(tidyverse)            # For data manipulation and visualization
library(ggrepel)              # For repelling overlapping text labels in plots

load('/code/volcanoplot/Vdata.Rdata') # Load the data file (Note: The comment indicates a request to change labels, etc.)

# Define color palette
Hcolor = c('#B9C19A', '#BAADC6', '#5B8D71', '#A92424', '#6376A0',  '#DDA7A7')

# Plot adjustments
p02 = ggplot() + 
  # Bars and names
  geom_col(data = df_summary, aes(x = cluster, y = 1.1*Max), width = 0.9, alpha = 0.05 ) +  # Upper limit bars
  geom_col(data = df_summary, aes(x = cluster, y = 1.1*Min), width = 0.9, alpha = 0.1) +  # Lower limit bars
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster), 
            data = df_nummary, color = "black", linewidth = .5, show.legend = FALSE) +  # Rectangles for clusters
  scale_fill_manual(values = Hcolor) +  # Custom colors as requested
  geom_text(data = df_summary, aes(label = scluster), x=df_summary$cluster,
            y = 0, vjust = 0.5, hjust = "center", size = 7, color = 'white') +  # Cluster labels
  labs(y = "log2FC", col = '-log10FDR', size = '-log10FDR') +  # Axis and legend labels
  # Points and theme
  geom_point(data = cdata, aes(x = jcluster, y = logFC, col = logPC, size = logPC),
             alpha = 0.9) +  # Points with jittering
  scale_color_gradientn(colors = c("#699ECA", "#FFCB5B", "#CC5B45", 'darkred'),
                        values = c(0, 0.25, 0.5, 1)) +  # Gradient color scale for points
  theme(
    panel.background = element_blank(), # Remove panel background
    axis.ticks.x = element_blank(), # Remove x-axis tick marks
    axis.text.x = element_blank(), # Remove x-axis text
    axis.title.x = element_blank(), # Remove x-axis title
    axis.line.y = element_line(color = "black"), # Keep y-axis line
    axis.title.y = element_text(size = 16),  # Adjust y-axis title font size
    axis.text.y = element_text(size = 12),   # Adjust y-axis tick label font size
    legend.title = element_text(size = 14),  # Adjust legend title font size
    legend.text = element_text(size = 12)    # Adjust legend item font size
  ) +
  # Labels
  ggrepel::geom_text_repel(   # Use geom_text_repel to avoid overlapping labels
    data = clabel,            # Data for labels
    max.overlaps = 20,        # Maximum allowed overlaps
    aes(x = jcluster, y = logFC, # fill = cluster, 
        label = Protein),        # Label mapping
    size = 5,                 # Label size
    alpha = 0.7,              # Label transparency
    fontface = "bold",        # Bold text
    box.padding = unit(0.3, "lines"), # Padding around labels
    min.segment.length = 0.2, # Minimum segment length for lines
    point.padding = 0.1,      # Distance between labels and points
    segment.color = 'black',  # Color of connecting lines
    show.legend = FALSE       # Do not show legend for this layer
  ) 

p02

# Save the plot as a PDF file
ggsave(filename = paste0('/results/Vpic.pdf') , 
       width = 12, height = 8) # Save the plot with specified dimensions

