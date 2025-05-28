install.packages("plotly")
install.packages("dplyr")
install.packages("rio")
library(plotly)       # For interactive 3D plotting
library(ggplot2)      # For general plotting
library(dplyr)        # For data manipulation

# Load data files
load(file="/code/manhattan/chr_lengths.Rdata")          # Chromosome length data
load(file="/code/manhattan/df.Rdata")                   # GWAS results to plot
data_all_lead_snp=rio::import("/code/manhattan/data_all_lead_snp.txt")  # Lead SNP information

# Clean and format trait names in lead SNP data
data_all_lead_snp$trait=data_all_lead_snp$trait%>% 
  gsub("-1e-5","",.) %>%    # Remove scientific notation suffixes
  gsub("_1e-5","",.) %>% 
  gsub("ARR","Arr",.) %>%   # Standardize abbreviation
  gsub("-"," ",.) %>%       # Replace hyphens with spaces
  gsub("_"," ",.)           # Replace underscores with spaces

# Create mapping ID combining trait and rsID
data_all_lead_snp$map_rs=paste0(data_all_lead_snp$trait,"_",data_all_lead_snp$rsID)

# Filter GWAS results to only include lead SNPs and select top 10 per trait by p-value
add_info=df_filtered_high %>% 
  filter(map_rs %in% data_all_lead_snp$map_rs)%>%
  group_by(trait) %>% 
  slice_max(logP,n=10)

# Create lead SNP dataframe with sequence numbers within each trait-chromosome group
lead_df=df_filtered_high[which(df_filtered_high$map_rs %in% add_info$map_rs),] %>% 
  group_by(trait, CHR) %>% 
  mutate(sequence = row_number())

# Set shape and size parameters for plotting
int_shape="diamond"  # Shape for lead SNPs
shape_map <- c("circle", "square", "diamond", "cross", "x")  # Available shapes
df_filtered_high$shape="circle"  # Default shape
df_filtered_high$shape[which(df_filtered_high$map_rs %in% add_info$map_rs)]=int_shape  # Set lead SNPs to diamond
df_filtered_high$size=1          # Default size
df_filtered_high$size[which(df_filtered_high$map_rs %in% add_info$map_rs)]=10  # Larger size for lead SNPs
# set order of 3d manhattan
trait_arrange=c("M50_ukb_CAD_meta.38_placo", 
                "M47_ukb_ARR_FEN.38_placo", "M47_ukb_CAD_meta.38_placo",
                "M47_ukb_CVD_FEN.38_placo", "M47_ukb_HF_meta.38_placo", 
                
                "M48_ukb_ARR_FEN.38_placo", "M48_ukb_CAD_meta.38_placo",
                "M48_ukb_CVD_FEN.38_placo", "M48_ukb_HF_meta.38_placo", 
                "M48_ukb_IHD_FEN.38_placo", "M48_ukb_MI_FEN.38_placo",
                
                "M51_ukb_ARR_FEN.38_placo","M51_ukb_CAD_meta.38_placo",
                "M51_ukb_CVD_FEN.38_placo","M51_ukb_HF_meta.38_placo",
                "M51_ukb_IHD_FEN.38_placo") %>% 
  gsub(".38_placo","",.) %>% 
  gsub("_ukb","",.) %>% 
  gsub("_FEN","",.) %>% 
  gsub("ARR","Arr",.) %>% 
  gsub("_meta","",.) %>% 
  gsub("_"," ",.)
# Create 3D Manhattan plot
p2=plot_ly(
  data = df_filtered_high,
  x = ~cum_BP,      # Cumulative base pair position
  y = ~trait,       # Phenotype/trait
  z = ~logP,        # -log10(p-value)
  color = ~trait,   # Color by trait
  type = "scatter3d",
  mode = "markers",
  marker = list(symbol = ~shape,  # Set marker shape
                size = ~size      # Set marker size
  )
) %>% layout(
  title = NULL,
  scene = list(
    xaxis = list(
      title = "CHR",
      tickvals = chr_starts,  # Chromosome start positions
      ticktext = c(1:22," "), # Chromosome labels
      tickangle=40,           # Label angle
      range = c(0, max(df_filtered_high$cum_BP)+100000000),
      gridwidth=1             # Grid line width
    ),
    yaxis = list(title = "trait",
                 tickvals = trait_arrange,  # Trait positions
                 ticktext=trait_arrange     # Trait labels
    ),
    zaxis = list(title = "-log10(P)"),     # Z-axis label
    aspectratio= list(x= 2.8, y= 3.8, z= 1),  # Aspect ratio
    camera = list(eye = list(x = 3, y =4, z = 3.5))  # Camera/view angle
  ),
  showlegend = FALSE  # Hide legend
)
# Prepare data for labeling top SNPs
high_logP_points <- df_filtered_high[which(df_filtered_high$map_rs %in% add_info$map_rs),] %>% 
  group_by(trait, CHR) %>% 
  mutate(plot_seq = row_number())  # Create sequence numbers within trait-chromosome groups

# Create backup of original data
high_logP_points_org=high_logP_points
high_logP_points=high_logP_points_org

# Calculate base pair differences between consecutive SNPs
high_logP_points$BP_dif <- c(0, diff(high_logP_points$cum_BP))

# Adjust plot sequence numbers to avoid overlapping labels
for(i in 2:nrow(high_logP_points)) {
  if(high_logP_points$trait[i] == high_logP_points$trait[i-1] && 
     high_logP_points$CHR[i] != high_logP_points$CHR[i-1] && 
     high_logP_points$BP_dif[i] < 200000000) {
    n_same=which(high_logP_points$trait==high_logP_points$trait[i] & 
                   high_logP_points$CHR==high_logP_points$CHR[i])
    high_logP_points$plot_seq[n_same] <- high_logP_points$plot_seq[i-1]+(1:length(n_same))
  }
}
high_logP_points=high_logP_points %>% arrange(trait, cum_BP, logP)

# Function to generate annotation positions for SNP labels
generate_seq_all_in_one=function(x_start_seq=-10,
                                 x_sep_seq=10,
                                 y_start_seq=-30,
                                 y_sep_seq=10){
  n <- nrow(high_logP_points)
  full_sequence <- seq_len(n)
  
  # Create list of annotation positions for each SNP
  annotations_list_last=lapply(full_sequence, function(i) {
    ax_now = x_start_seq+(x_sep_seq*high_logP_points$plot_seq[i])  # X offset
    ay_now = y_start_seq+(y_sep_seq*high_logP_points$plot_seq[i])  # Y offset
    
    # Determine anchor positions based on offset direction
    xanchor_now = if_else(ax_now>0,"left","right")
    yanchor_now = if_else(ay_now>0,"top","bottom")
    if(ay_now==0) yanchor_now = "middle"
    if(ax_now==0) xanchor_now = "middle"
    
    # Return annotation parameters for this SNP
    res=list(
      x = high_logP_points$cum_BP[i],      # SNP position on x-axis
      y = high_logP_points$trait[i],       # SNP position on y-axis
      z = high_logP_points$logP[i],        # SNP position on z-axis
      text = high_logP_points$SNP[i],      # SNP ID to display
      textangle = 0,                       # Text angle
      ax = ax_now,                         # Arrow x offset
      ay = ay_now,                         # Arrow y offset
      font = list(color = "black", size = 12),  # Text formatting
      arrowcolor = "grey",                 # Arrow color
      arrowsize = 1,                       # Arrow size
      arrowwidth = 1,                      # Arrow width
      arrowhead = 1,                       # Arrow head size
      xanchor = xanchor_now,                # Text x anchor
      yanchor = yanchor_now                 # Text y anchor
    )
    return(res)
  })
  return(annotations_list_last)
}

# Generate annotation positions with specific spacing parameters
fine_tune1=generate_seq_all_in_one(x_start_seq=12,
                                   x_sep_seq=-10,
                                   y_start_seq=12,
                                   y_sep_seq=-10)

# Add annotations to the plot and adjust camera/view angle
p7=p2 %>% layout(
  scene=list(annotations = fine_tune1,  # Add SNP labels
             xaxis = list(tickangle=0), # Adjust x-axis label angle
             yaxis = list(title = "trait", tickangle=0),  # Adjust y-axis
             camera = list(eye = list(x = 3, y =4, z = 3.5))
  )
)
p7
htmlwidgets::saveWidget(widget = p7, file = "/results/manhattanplot.html", selfcontained = TRUE)
