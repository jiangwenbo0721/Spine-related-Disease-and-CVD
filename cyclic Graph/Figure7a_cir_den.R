# back --------------------------------------------------------------------
# Set working directory and load required libraries
install.packages("ggraph")
install.packages("igraph")
install.packages("ggsci")
library(dplyr)        # Data manipulation
library(tidyr)        # Data tidying
library(stringr)      # String operations
library(data.table)   # Fast data manipulation
library(ggplot2)      # Data visualization
library(ggraph)       # Graph visualization
library(igraph)       # Network analysis
library(tidyverse)    # Collection of tidyverse packages
library(ggsci)        # Color palettes for ggplot2
library(RColorBrewer) # Color palettes

# Custom function to capitalize first word of a string
capitalize_first_word <- function(text) { 
  str_to_title(substring(text, 1, 1)) %>% paste0(substring(text, 2))
}

# data --------------------------------------------------------------------
# Import interaction data file
res_file = rio::import("/code/Figure7a/Spine_CVD_circular_diagram.txt")

# Process interaction data
inter_file = res_file %>%
  dplyr::select(c("phen1","phen2","Region","rsID","uniqID","GenomicLocus"))
inter_file$phen2 = inter_file$phen2 %>% gsub("ARR","Arr",.)  # Standardize phenotype names
inter_file = inter_file %>% separate("uniqID", c("CHR","POS","other","effect"))  # Split unique ID into components
inter_file$map_rs = paste0(inter_file$phen1," ",inter_file$phen2,"_",inter_file$rsID)  # Create mapping ID

# Load filtered high data and merge with interaction data
load("/code/Figure7a/df.Rdata")
df_filtered_high_int = df_filtered_high[which(df_filtered_high$map_rs %in% inter_file$map_rs),]
df_filtered_high_int$Direction = ifelse(df_filtered_high_int$T.placo < 0, -1, 1)  # Add direction based on T.placo
df_filtered_high_int = df_filtered_high_int %>% dplyr::select(c("map_rs","Direction"))

# Merge interaction data with direction information
inter_file = inter_file %>% left_join(df_filtered_high_int)
inter_file$LeadSNPs = inter_file$rsID  # Create LeadSNPs column
data_all = inter_file  # Final data for analysis

# process -----------------------------------------------------------------
# Process phenotype names and SNP information
data_all$phen2 = data_all$phen2 %>% capitalize_first_word()  # Capitalize first letter
data_all$LeadSNPs_org = data_all$LeadSNPs  # Keep original SNP names
n_to_change = grep("rs", data_all$LeadSNPs)  # Identify rsIDs to modify
data_all$LeadSNPs[n_to_change] = data_all$LeadSNPs[n_to_change] %>% 
  paste0(":", data_all$other[n_to_change], ":", data_all$effect[n_to_change])  # Add effect info to rsIDs

# Import and process final dataframe with gene annotations
final_df = rio::import(file="/code/Figure7a/final_df.txt")
final_df = unique(final_df)  # Remove duplicates

# Handle missing values in gene annotation columns
final_df$eqtl_proofed_symbols[which(final_df$eqtl_proofed_symbols == "")] = NA
final_df$posMap_proofed_symbols[which(final_df$posMap_proofed_symbols == "")] = NA
final_df$ci_proofed_symbols[which(final_df$ci_proofed_symbols == "")] = NA

# Split multiple gene annotations into separate rows
final_df = final_df %>% 
  separate_longer_delim(c(eqtl_proofed_symbols, posMap_proofed_symbols), delim="; ") %>% 
  separate_longer_delim(c(ci_proofed_symbols), delim="; ")

# Add prefixes to gene annotation types
final_df$eqtl_proofed_symbols[!is.na(final_df$eqtl_proofed_symbols)] = paste0("eQTL: ", final_df$eqtl_proofed_symbols[!is.na(final_df$eqtl_proofed_symbols)])
final_df$ci_proofed_symbols[!is.na(final_df$ci_proofed_symbols)] = paste0("ci: ", final_df$ci_proofed_symbols[!is.na(final_df$ci_proofed_symbols)])
final_df$posMap_proofed_symbols[!is.na(final_df$posMap_proofed_symbols)] = paste0("posMap: ", final_df$posMap_proofed_symbols[!is.na(final_df$posMap_proofed_symbols)])

# Create long format of gene annotations
final_df_longer = final_df %>% 
  dplyr::select(c("pheid","GenomicLocus","posMap_proofed_symbols",
                  "eqtl_proofed_symbols","ci_proofed_symbols")) %>%
  pivot_longer(c("posMap_proofed_symbols","eqtl_proofed_symbols","ci_proofed_symbols"), values_to = "Nearest_genes1") %>% 
  dplyr::select(c(pheid, GenomicLocus, Nearest_genes1))
final_df_longer$final_id = paste0(final_df_longer$pheid, "_", final_df_longer$GenomicLocus)
final_df_longer = final_df_longer %>% dplyr::select(final_id, Nearest_genes1) %>% na.omit()
final_df_longer$final_id = final_df_longer$final_id %>% gsub("stroke","Stroke",.)  # Standardize phenotype name

# plot --------------------------------------------------------------------
# Prepare data for plotting
data_test = data_all
data_test$final_id = paste0(data_test$phen1, "_", data_test$phen2, "_", data_test$GenomicLocus)
data_test = data_test %>% left_join(final_df_longer)  # Merge with gene annotations
data_test$Nearest_genes1[which(is.na(data_test$Nearest_genes1))] = "None"  # Fill missing genes
data_test = data_test[-which(data_test$Nearest_genes1 == "None"),]  # Remove rows without gene annotations
col_to_plot = c("phen1","phen2","Region","LeadSNPs","Nearest_genes1")  # Columns to include in plot

# Custom function to rename and format data for hierarchical plotting
rename_multi_subpoint = function(data_test, col_to_plot = c()) {
  data_renamed = data_test %>% dplyr::select(all_of(col_to_plot))
  
  # Function to generate new columns with concatenated values
  generate_new_columns <- function(df, n) {
    for (i in 2:n) {
      new_col_name <- paste("new_col", i, sep = "_")
      df[[new_col_name]] <- apply(df[, 1:i], 1, paste, collapse = "_")
    }
    return(df)
  }
  
  data_renamed <- generate_new_columns(data_renamed, 4)
  data_renamed$new_col_1 = data_renamed$phen1
  
  # Function to add spaces to duplicate regions for better visualization
  add_spaces <- function(df) {
    all_category = unique(df$category)
    duped_region = df$region[duplicated(df$region)] %>% unique()
    dup_category_inter = list()
    
    for (i in 1:length(duped_region)) {
      dup_category_inter[[i]] = df$category[which(df$region == duped_region[i])] %>% unique()
      names(dup_category_inter)[i] = duped_region[i]
    }
    
    # Add spaces to duplicate regions to separate them visually
    for (i in 1:nrow(df)) {
      region <- df$region[i]
      category <- df$category[i]
      if(region %in% duped_region){
        dup_category_sub = dup_category_inter[[which(names(dup_category_inter) == region)]] %>% unique()
        n_cate = which(dup_category_sub == category) - 1
        spaces <- strrep(" ", n_cate)
        df$region[i] = paste0(df$region[i], spaces)
      }
    }
    return(df)
  }
  
  # Process each level of the hierarchy
  i = 1
  data_last = data.frame()
  for (i in 1:4) {
    colum_to_select <- c(paste("new_col", i, sep = "_"), col_to_plot[i+1])
    data_inter = data_renamed[, colum_to_select]
    colnames(data_inter) = c("category","region")
    data_inter <- data_inter %>% add_spaces()
    colnames(data_inter) = c("main", col_to_plot[i+1])
    
    if(i == 1){
      data_last = data_inter
    } else {
      data_last = data_last %>% cbind(data_inter[,2])
    }
    colnames(data_last)[i+1] = col_to_plot[i+1]
  }
  return(data_last)
}

# Apply renaming function to prepare final plotting data
data_test_final = data_test %>% rename_multi_subpoint(col_to_plot)
data_test_final$Direction = data_test$Direction
colnames(data_test_final)[1] = "phen1"

# plot preparation--------------------------------------------------------------------
# Prepare data for network visualization
data <- data_test_final

# Melt data to long format for edge creation
melted_data <- data %>%
  dplyr::select(c("phen1","phen2","Region","LeadSNPs","Nearest_genes1")) %>% 
  pivot_longer(c("phen2","Region","LeadSNPs","Nearest_genes1"), values_to = "to") %>% 
  dplyr::select(c(phen1, to))

# Add self-references for root nodes
melted_data = melted_data %>% 
  rbind(data.frame(phen1 = unique(data$phen1), to = unique(data$phen1)))

# Create edges for each level of the hierarchy
edges1 <- data.frame(from = "Spine Disorders", to = unique(data$phen1))  # Root to phenotypes
edges2 <- data.frame(from = data$phen1, to = data$phen2)                # Phenotypes to sub-phenotypes
edges3 <- data.frame(from = data$phen2, to = data$Region)               # Sub-phenotypes to regions
edges4 <- data.frame(from = data$Region, to = data$LeadSNPs)            # Regions to SNPs
edges5 <- data.frame(from = data$LeadSNPs, to = data$Nearest_genes1)    # SNPs to genes

# Combine all edges and remove duplicates
edges <- rbind(edges1, edges2, edges3, edges4, edges5)
edges <- distinct(edges)

# plot -------------------------------------------------------------------
# Create node data frame with visualization attributes
vertices_name <- unique(c(as.character(edges$from), as.character(edges$to)))
node_size_scale = 2

# Set visual attributes for different node types
shape_pre = rep(16, length(vertices_name))  # Default shape
shape_pre[which(vertices_name %in% data$LeadSNPs)] = 15  # Different shape for SNPs
size_pre = rep(1, length(vertices_name))    # Default size
size_pre[which(vertices_name %in% data$LeadSNPs)] = 1    # Size for SNPs
text_size_pre = rep(1.03, length(vertices_name))  # Default text size
text_size_pre[which(vertices_name %in% data$LeadSNPs)] = 1.05  # Larger text for SNPs

vertices <- data.frame(
  name = vertices_name,
  node_size = size_pre,
  text_size = text_size_pre,
  shape = shape_pre,
  name_org = vertices_name %>% gsub(" ","",.)  # Name without spaces
)

# Assign groups for coloring
rownames(vertices) <- vertices$name
vertices$group <- melted_data$phen1[match(vertices$name, melted_data$to)]
n_rs = which(vertices$name %in% data$LeadSNPs)
vertices$group[n_rs] <- data$Direction[match(vertices$name[n_rs], data$LeadSNPs)]  # Color SNPs by direction

# Create hierarchical graph object
hierarchy <- graph_from_data_frame(edges, vertices = vertices, directed = TRUE)

# Define color palette
library(colorspace)
expand_n = 1.2
colorspecs = c("blue", "red", "#80B1D3", "#FDB462", "#B3DE69", 
               "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F","darkgreen")

# Create adjusted version of the plot with scaled node positions
p <- ggraph(hierarchy, layout = 'dendrogram', circular = TRUE)

# Adjust positions of different node types for better spacing
scale_region = 0.72
scale_rsid = 0.67
p$data$x[which(p$data$name %in% data$phen1)] = 0.7 * p$data$x[which(p$data$name %in% data$phen1)]
p$data$y[which(p$data$name %in% data$phen1)] = 0.7 * p$data$y[which(p$data$name %in% data$phen1)]
p$data$x[which(p$data$name %in% data$phen2)] = scale_rsid * p$data$x[which(p$data$name %in% data$phen2)]
p$data$y[which(p$data$name %in% data$phen2)] = scale_rsid * p$data$y[which(p$data$name %in% data$phen2)]
p$data$x[which(p$data$name %in% data$Region)] = 0.68 * p$data$x[which(p$data$name %in% data$Region)]
p$data$y[which(p$data$name %in% data$Region)] = 0.68 * p$data$y[which(p$data$name %in% data$Region)]
p$data$x[which(p$data$name %in% data$LeadSNPs)] = scale_rsid * p$data$x[which(p$data$name %in% data$LeadSNPs)]
p$data$y[which(p$data$name %in% data$LeadSNPs)] = scale_rsid * p$data$y[which(p$data$name %in% data$LeadSNPs)]
p$data$x[which(p$data$name %in% data$Nearest_genes1)] = scale_region * p$data$x[which(p$data$name %in% data$Nearest_genes1)]
p$data$y[which(p$data$name %in% data$Nearest_genes1)] = scale_region * p$data$y[which(p$data$name %in% data$Nearest_genes1)]
p$data$name_org = p$data$name %>% gsub(" ","",.)  # Remove spaces from names

# Create final adjusted plot
expand_n = 1.1
p1 = p + geom_edge_diagonal(aes(x = x*1, y = y*1), color = 'grey', alpha = 1) +
  geom_node_text(aes(x = x*1.06, y = y*1.06, label = name_org,
                     color = "black", size = 1,
                     hjust = "outward",
                     angle = -((-node_angle(x,y)+90) %% 180)+ 90), 
                 color = "black", alpha = 1, fontface = "bold") +
  geom_node_point(aes(x = x*1.03, y = y*1.03, shape = as.factor(shape), 
                      color = group, size = 1), stroke = 0.5) +
  scale_color_manual(values = colorspecs) +
  scale_shape_manual(values = c(15, 16)) +
  theme_void() +
  theme(legend.position = "right") +
  expand_limits(x = c(-1*expand_n, expand_n), y = c(-1*expand_n, expand_n)) +
  guides(size = guide_legend("Node Size"), shape = guide_legend("Shape"))

# Remove legends for final output
p1 = p1 + guides(size = "none", shape = "none", colour = "none")

# Save adjusted plot
ggsave(p1, file = "/results/final.pdf", height = 50, width = 50, units = "cm")

