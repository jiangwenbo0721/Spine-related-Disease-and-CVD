library(ggplot2)
library(readxl)
library(dplyr)

# Read the complete dataset
data_cox <- read_excel("/code/ProteinanddiseaseForest/plot_Protein_and_disease_cox_result_data.xlsx", sheet = 1)

# Calculate the height of each rectangle as half of the confidence interval
data_cox <- data_cox %>%
  mutate(tile_height = (conf.high - conf.low) / 2)

# Custom color palette
my_colors <- c('#E64B35FF', "#4DBBD5FF", "#00A087FF", 
               "#F39B7FFF", "#91D1C2FF", "#8491B4FF", 
               "#7E6148FF", "#e47dc6")

# Get all protein names
protein_list <- unique(data_cox$protein)

# Loop through each protein
for (prot in protein_list) {
  data_subset <- filter(data_cox, protein == prot)
  
  p <- ggplot(data_subset) + 
    geom_errorbar(aes(x = outcome, ymin = conf.low, ymax = conf.high, color = outcome), 
                  width = 0.1, size = 0.2) + 
    geom_tile(aes(x = outcome, y = estimate, 
                  width = 0.3, height = tile_height, fill = outcome), 
              color = NA) + 
    geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 1) +
    scale_fill_manual(values = my_colors) +
    scale_color_manual(values = my_colors) +
    scale_y_continuous() +
    scale_x_discrete(expand = c(0.1, 0)) +
    xlab("") + 
    ggtitle(paste("Forest Plot -", prot)) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "bold", color = "black"),
      axis.text.y = element_text(color = "black", size = 12),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  # Save the plot as a PDF file
  ggsave(filename = paste0("/results/forest_", prot, ".pdf"), plot = p, width = 5, height = 4)
}

