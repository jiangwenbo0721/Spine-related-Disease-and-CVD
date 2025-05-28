# Load required libraries
install.packages("reshape2")
library(ggplot2)
library(reshape2)  # For data reshaping

# Create sample fold change data (replace with your actual data)
foldchange_values <- c(-0.0274, -0.0434, -0.0387, -0.0398, -0.0019, -0.0504, 
                       -0.0586, 0.0848, 0.1790, -0.0121,
                       0.0238, 0.0112, 0.0278,
                       0.0659, 0.1168, 0.1105, 0.1159, 0.0885, 0.0980, 0.1182, 
                       0.0792, 0.0954, 0.0793)

# Create data frame with fold change values
df <- data.frame(FoldChange = foldchange_values)

# Gene names corresponding to the fold change values
gene_names <- c("X23577.0.0", "X23437.0.0", "X23434.0.0", "X23430.0.0", "X30060.0.0", 
                "X30120.0.0", "X23444.0.0", "X30130.0.0", "X30820.0.0", "X23469.0.0",
                "X23460.0.0", "X23461.0.0", "X23468.0.0",
                "X23408.0.0", "X23481.0.0", "X23482.0.0", "X23483.0.0", "X23484.0.0", 
                "X23486.0.0", "X23487.0.0", "X23488.0.0", "X23494.0.0", "X30290.0.0")

# Assign gene names to data frame
df$Gene <- gene_names 

# Convert data frame to long format for plotting
df_melted <- melt(df, id.vars = "Gene")

############################
# Subset analysis: Positive and Negative Fold Changes

# Extract genes with positive fold changes
df_positive <- df[df$FoldChange > 0, ]

# Extract genes with negative fold changes
df_negative <- df[df$FoldChange < 0, ]

# Convert positive subset to long format
df_melted_positive <- melt(df_positive, id.vars = "Gene")

# Generate heatmap for positive fold changes (white to red gradient)
ggplot(df_melted_positive, aes(x = Gene, y = variable, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "white", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Genes") +
  ylab("Fold Change") +
  labs(title = "Positive Fold Change Heatmap",
       fill = "logFC")

# Convert negative subset to long format
df_melted_negative <- melt(df_negative, id.vars = "Gene")

# Generate heatmap for negative fold changes (blue to white gradient)
p=ggplot(df_melted_negative, aes(x = Gene, y = variable, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "white", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Genes") +
  ylab("Fold Change") +
  labs(title = "Negative Fold Change Heatmap",
       fill = "logFC")

ggsave("/results/metabolomics Heatmap.png", plot = p, width = 6, height = 4, dpi = 300)

