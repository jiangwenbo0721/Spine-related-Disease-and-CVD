# Load required R packages
install.packages("funkyheatmap")
install.packages("Cairo")

library(tidyverse)            # For data manipulation and visualization
library(funkyheatmap)         # For creating heatmap with additional features

# Load the data file
load('/code/funkyheatmap/Hdata-N.Rdata') # Load the data file containing necessary data and variables

# Create a funky heatmap
p01 = funky_heatmap(
  data          = j01_CSEA,           # Data frame to be used for the heatmap
  column_info   = j01_column_info,    # Information for column annotations
  palettes      = j01_palettes,       # Color palettes for the heatmap
  position_args = position_arguments( # Adjustments for column annotations
    col_annot_offset = 13,            # Offset for column annotations (adjust if column names are long)
    col_annot_angle = 50,             # Angle for column annotations (adjust if column names are long)
  ), 
  scale_column = F 
) # This creates the heatmap with legends

# Save the heatmap as a PDF file
Cairo::CairoPDF( file = "/results/Hpic-N.pdf",         # Use Cairo to create a high-quality PDF
                 width = unit(15,'cm'),     # Set the width of the PDF
                 height = unit(15,'cm')     # Set the height of the PDF
)
p01 ; dev.off() 
