install.packages("forestplot")
# Load data
df <- read.csv("/code/forest/Figure3_E_forest.csv")
head(df)
# Load required packages
library(ggplot2)
library(forestplot)
library(grid)
# Split CI into low and high values
df$CI_low <- as.numeric(sub("-.*", "", df$CI))
df$CI_high <- as.numeric(sub(".*-", "", df$CI))
# Prepare label text: Outcome | OR (95% CI) | P-value
labeltext <- cbind(
  c("Outcome", as.character(df$Outcome)),
  c("OR (95% CI)", sprintf("%.3f (%.3f-%.3f)", df$OR, df$CI_low, df$CI_high)),
  c("P-value", sprintf("%.3f", df$P_value))
)
# Create forest plot
p=forestplot(
  labeltext = labeltext,
  mean = c(NA, df$OR),         # OR values
  lower = c(NA, df$CI_low),    # CI lower bounds
  upper = c(NA, df$CI_high),   # CI upper bounds
  zero = 1,                    # Reference line at OR=1
  graph.pos = 3,               # Forest plot position (3rd column)
  boxsize = 0.2,               # Size of OR points
  lineheight = "auto",         # Automatic line height
  colgap = unit(5, "mm"),      # Column spacing
  lwd.zero = 2,                # Reference line width
  lwd.ci = 2,                  # CI line width
  col = fpColors(box = "royalblue", line = "darkblue"),  # Colors
  
  # Graphic parameters for points (circular shape)
  gp = gpar(fill = "royalblue", col = "darkblue", lwd = 1.5),
  
  # Text formatting
  txt_gp = fpTxtGp(
    label = gpar(cex = 0.9),   # Label text size
    ticks = gpar(cex = 0.8),   # Axis tick size
    xlab = gpar(cex = 1)       # X-axis label size
  )
)
p
png("/results/forest.png", width = 6, height = 4, units = "in", res = 300)
print(p)
dev.off()
