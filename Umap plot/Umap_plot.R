rm(list=ls())
gc()
# Load
# devtools::install_github("sajuukLyu/ggunchull", type = "source")
library(ggplot2)
library(readxl)
library(dplyr)
library(umap)
library(ggunchull)
library(cluster)
library(clusterCrit) 
library(factoextra)
library(mclust)
library(ggrepel)
library(tidyverse)
library(gplots)
library(proxy)
library(fpc)
library(clusterSim)
library(parallel)
library(doParallel)
library(missRanger)
library(ranger)
library(Rtsne) 
library(dbscan)
library(plotly)
library(ggpubr)
library(ggforce)
library(concaveman)
library(uwot)
# -------------------- UMAP visualization of original data clustering
getwd()
setwd("")

mean_data_5group <- read_excel("spinediseaseM45_M51_estimate_nocov.xlsx", sheet=1)

# Data Preprocessing
clustering_data <- mean_data_5group %>%
  select(-protein) 

# Elbow Rule
wss <- numeric(10)
for (k in 1:10) {
  km_temp <- kmeans(clustering_data, centers = k, nstart = 25, iter.max = 50)
  wss[k] <- km_temp$tot.withinss
}

# Draw the elbow rule diagram (Observe the optimal K value)
#ggplot(data.frame(K = 1:10, WSS = wss), aes(x = K, y = WSS)) +
  #geom_line() + geom_point() +
  #labs(title = "Elbow Method on Standardized Data (K=4)", 
       #x = "Number of Clusters", 
       #y = "Total Within-Cluster SS") +
  #theme_minimal()

# Perform K-means clustering
set.seed(123)
km <- kmeans(clustering_data, centers = 4, nstart = 25)

# UMAP Dimension Reduction (also using standardized data)
set.seed(123)  
umap_res <- umap(clustering_data, 
                 n_components = 2,
                 n_neighbors = 15,
                 min_dist = 0.1)

# Define the save path
save_path <- "D:/clustering_umap.pdf"

# Open the PDF device and set the canvas size
pdf(save_path, width = 10, height = 8)

# Draw the clustering results
ggplot(data = data.frame(UMAP1 = umap_res[, 1],  
                         UMAP2 = umap_res[, 2],  
                         Cluster = factor(km$cluster)),  
       aes(x = UMAP1, y = UMAP2, color = Cluster)) +
  geom_point(size = 5, alpha = 0.8) +  
  theme_minimal() +
  labs(title = "K-means Clustering of Proteins (UMAP)",
       x = "UMAP Dimension 1",
       y = "UMAP Dimension 2") +
  scale_color_manual(name = "Cluster", 
                     values = c("1" = "#CA2928",  
                                "2" = "#2171A9",  
                                "3" = "#96CB85",  
                                "4" = "#FFBB78"))

dev.off()


# Add the UMAP dimensionality reduction results and clustering labels back to mean_data_5group
mean_data_5group_with_cluster <- mean_data_5group %>%
  dplyr::mutate(UMAP1 = umap_res[, 1],
                UMAP2 = umap_res[, 2],
                Cluster = km$cluster)

#savemean_data_5group_with_cluster
write.csv(mean_data_5group_with_cluster, 
          file = "mean_data_5group_with_cluster.csv", 
          row.names = FALSE)


#########It is classified into five categories according to protein
#########It is classified into five categories according to protein
#########It is classified into five categories according to protein

# Read Data
M45 <- read.csv("unique_proteins_M45_status.csv")
M45_proteins <- unique(M45$Protein)

M46 <- read.csv("unique_proteins_M46_status.csv")
M46_proteins <- unique(M46$Protein)

M47 <- read.csv("unique_proteins_M47_status.csv")
M47_proteins <- unique(M47$Protein)

M48 <- read.csv("unique_proteins_M48_status.csv")
M48_proteins <- unique(M48$Protein)

M50 <- read.csv("unique_proteins_M50_status.csv")
M50_proteins <- unique(M50$Protein)

M51 <- read.csv("unique_proteins_M51_status.csv")
M51_proteins <- unique(M51$Protein)

common_proteins <- read.csv("common_proteins.csv")
common_protein <- unique(common_proteins$Protein)

# Read the master data
mean_data_5group <- read_excel("spinediseaseM45_M51_estimate_nocov.xlsx", sheet = 1)

# Assign a category label to each protein
mean_data_5group <- mean_data_5group %>%
  mutate(Category = case_when(
    protein %in% M45_proteins ~ "M45",
    protein %in% M46_proteins ~ "M46",
    protein %in% M47_proteins ~ "M47",
    protein %in% M48_proteins ~ "M48",
    protein %in% M50_proteins ~ "M50",
    protein %in% M51_proteins ~ "M51",
    protein %in% common_protein ~ "Common",
    TRUE ~ "Other"  # Proteins that do not fall into any of the above categories
  ))

# data preprocessing
clustering_data <- mean_data_5group %>%
  select(-protein, -Category) 

# Elbow Rule
wss <- numeric(10)
for (k in 1:10) {
  km_temp <- kmeans(clustering_data, centers = k, nstart = 25, iter.max = 50)
  wss[k] <- km_temp$tot.withinss
}

#Draw the elbow rule diagram (observe the optimal K value)
#ggplot(data.frame(K = 1:10, WSS = wss), aes(x = K, y = WSS)) +
  #geom_line() + geom_point() +
  #labs(title = "Elbow Method on Standardized Data (K=4)", 
      # x = "Number of Clusters", 
      # y = "Total Within-Cluster SS") +
  #theme_minimal()

# Perform K-means clustering
set.seed(123)
km <- kmeans(clustering_data, centers = 4, nstart = 25)

# UMAP dimension reduction (also using standardized data)
set.seed(123)  
umap_res <- umap(clustering_data, 
                 n_components = 2,
                 n_neighbors = 15,
                 min_dist = 0.1)

# Define color mapping
color_mapping <- c(
  "M45" = "#B9C19A",       
  "M46" = "#BAADC6",      
  "M47" = "#5B8D71",       
  "M48" = "#BA5050",      
  "M50" = "#6376A0",      
  "M51" = "#DDA7A7",      
  "Common" = "#57AF37",   
  "Other" = "#F0F0F0"      
)

color_mapping <- c(
  "M45" = "#7A8C6E",      
  "M46" = "#9E8FB2",      
  "M47" = "#3A7A5F",      
  "M48" = "#BA5050",      
  "M50" = "#6376A0",      
  "M51" = "#DDA7A7",      
  "Common" = "#FFBB78",   
  "Other" = "#F0F0F0"     
)



color_mapping <- c(
  "M45" = "#B9C19A",      
  "M46" = "#9E8FB2",       
  "M47" = "#2A6A4F",       
  "M48" = "#BA5050",       
  "M50" = "#6376A0",       
  "M51" = "#DDA7A7",       
  "Common" = "#FFBB78",    
  "Other" = "#F0F0F0"      
)
# Define the save path
save_path <- "D:/5_clustering_umap04220516.pdf"

# Open the PDF device and set the canvas size
pdf(save_path, width = 10, height = 8)  # Width: 10 inches, Height: 8 inches

# Draw the clustering results
# Define transparency mapping
alpha_mapping <- ifelse(mean_data_5group$Category == "Other", 0.2, 1)  # Other is 0.2, and that of others is 1

# Draw the clustering results
ggplot(data = data.frame(UMAP1 = umap_res[, 1],  # Extract the first dimension of UMAP
                         UMAP2 = umap_res[, 2],  # Extract the second dimension of UMAP
                         Category = mean_data_5group$Category),  # Protein Categories
       aes(x = UMAP1, y = UMAP2, color = Category)) +
  geom_point(size = 5, alpha = alpha_mapping) +  
  ggunchull::stat_unchull(
    aes(fill = Category),  
    alpha = 0,             
    linewidth = 1,         
    delta = .5,            
    lty = 2                
  ) +
  theme_minimal() +
  labs(title = "Protein Clustering by Disease Category (UMAP)",
       x = "UMAP Dimension 1",
       y = "UMAP Dimension 2") +
  scale_color_manual(values = color_mapping) +  
  scale_fill_manual(values = color_mapping)    

# Turn off the PDF device
dev.off()



#It is classified into four types according to protein
#It is classified into four types according to protein
#It is classified into four types according to protein

# Read Data
M45 <- read.csv("unique_proteins_M45_status.csv")
M45_proteins <- unique(M45$Protein)

M46 <- read.csv("unique_proteins_M46_status.csv")
M46_proteins <- unique(M46$Protein)

M47 <- read.csv("unique_proteins_M47_status.csv")
M47_proteins <- unique(M47$Protein)

M48 <- read.csv("unique_proteins_M48_status.csv")
M48_proteins <- unique(M48$Protein)

M50 <- read.csv("unique_proteins_M50_status.csv")
M50_proteins <- unique(M50$Protein)

M51 <- read.csv("unique_proteins_M51_status.csv")
M51_proteins <- unique(M51$Protein)

common_proteins <- read.csv("common_proteins.csv")
common_protein <- unique(common_proteins$Protein)

M47_M51 <- read.csv("M47_M51common_proteins.csv")
M47_M51_proteins <- unique(M47_M51$Protein)

# Read Data
mean_data_5group <- read_excel("spinediseaseM45_M51_estimate_nocov.xlsx", sheet = 1)

# Assign a category label to each protein
mean_data_5group <- mean_data_5group %>%
  mutate(Category = case_when(
    protein %in% M45_proteins ~ "M45",
    protein %in% M46_proteins ~ "M46",
    protein %in% M47_M51_proteins ~ "M47_M51",
    protein %in% common_protein ~ "Common",
    TRUE ~ "Other"  #Proteins that do not fall into any of the above categories
  ))

# Data Preprocessing
clustering_data <- mean_data_5group %>%
  select(-protein, -Category)  # Remove the protein name and category column

# Elbow Rule
wss <- numeric(10)
for (k in 1:10) {
  km_temp <- kmeans(clustering_data, centers = k, nstart = 25, iter.max = 50)
  wss[k] <- km_temp$tot.withinss
}

#Draw the elbow rule diagram (observe the optimal K value)
#ggplot(data.frame(K = 1:10, WSS = wss), aes(x = K, y = WSS)) +
  #geom_line() + geom_point() +
  #labs(title = "Elbow Method on Standardized Data (K=4)", 
       #x = "Number of Clusters", 
       #y = "Total Within-Cluster SS") +
 # theme_minimal()

# Perform K-means clustering
set.seed(123)
km <- kmeans(clustering_data, centers = 4, nstart = 25)

# UMAP dimension reduction (also using standardized data)
set.seed(123)  
umap_res <- umap(clustering_data, 
                 n_components = 2,
                 n_neighbors = 15,
                 min_dist = 0.1)

# Define color mapping
color_mapping <- c(
  "M45" = "#B9C19A",
  "M46" = "#BAADC6",
  "M47_M51" = "#DC143C",
  "Common" = "#57AF37",
  "Other" = "#F0F0F0"  
)

# Define color mapping
color_mapping <- c(
  "M45" = "#7A8C6E",
  "M46" = "#9E8FB2",
  "M47_M51" = "#DC143C",
  "Common" = "#FFBB78",
  "Other" = "#F0F0F0"  
)


# Define color mapping
color_mapping <- c(
  "M45" = "#B9C19A",
  "M46" = "#9E8FB2",
  "M47_M51" = "#DC143C",
  "Common" = "#FFBB78",
  "Other" = "#F0F0F0"  
)
# Define the save path
save_path <- "D:/4_clustering_umap04220516.pdf"

# Open the PDF device and set the canvas size
pdf(save_path, width = 10, height = 8)  # Width: 10 inches, Height: 8 inches

# Draw the clustering results
# Define transparency mapping
alpha_mapping <- ifelse(mean_data_5group$Category == "Other", 0.2, 1)  # Other is 0.2, and that of others is 1

# Draw the clustering results
ggplot(data = data.frame(UMAP1 = umap_res[, 1],  # Extract the first dimension of UMAP
                         UMAP2 = umap_res[, 2],  # Extract the second dimension of UMAP
                         Category = mean_data_5group$Category),  # Protein Categories
       aes(x = UMAP1, y = UMAP2, color = Category)) +
  geom_point(size = 5, alpha = alpha_mapping) +  
  ggunchull::stat_unchull(
    aes(fill = Category),  
    alpha = 0,             
    linewidth = 1,        
    delta = .5,            
    lty = 2                
  ) +
  theme_minimal() +
  labs(title = "Protein Clustering by Disease Category (UMAP)",
       x = "UMAP Dimension 1",
       y = "UMAP Dimension 2") +
  scale_color_manual(values = color_mapping) + 
  scale_fill_manual(values = color_mapping)    

# Turn off the PDF device
dev.off()



