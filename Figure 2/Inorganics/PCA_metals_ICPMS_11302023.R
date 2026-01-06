#calling the necessary packages:
library(ggplot2)
library(stringr)
library(data.table)
library(readr)
library(dplyr)
library(tidyr)
library(reshape)     
library(gplots)
library(pheatmap)
library(heatmap3)
library(magrittr) 
library(tibble) 
library(RColorBrewer)
library(paletteer)
library(clr)
library(tidyverse) #used for data science. The eight core packages inside this library are: ggplot2 (data visualisation), dplyr (data manipulation), tidyr, readr, purrr, tibble, stringr, and forcats
library(KODAMA) # to use the normalisation function
library(ggrepel) #mainly used to repel overlapping text labels in ggplots
library(vegan) #popular library for analysing ecological diversity and for multivariate analysis of community data. Here, we use it for PCoA
library(svglite) #to save the plots in support vector graphics (svg) format
library(factoextra) #for extracting and visualizing outputs of multivariate analyses such as PCA, k-means
library(ggsci) #provides color palettes for ggplot2 that can be used for scientific journals
library(matrixStats) #contains highly optimized functions to perform statistics on matrix data
library(cowplot) #efficient functions to arrange several plots

################## PCA - all timepoints - O16 and 018 ############################
setwd("C:/Snowbowl_Data/center scaled metals and PCA code")



metadata<- read.csv("metadata_snowOnly_11292023_metals.csv")
clr_metals <- read.csv("snow_only_metals_CLR_scaled_snowbowl_11302023.csv")

data_merge <- metadata %>% dplyr::select("filename", "ATTRIBUTE_Snow_Type") %>% 
  left_join(clr_metals) %>% 
  column_to_rownames("filename")

# filter data in sample type column
data_merge <- data_merge %>% filter(ATTRIBUTE_Snow_Type != "SMRW_Low")
  

#specify the number of colors needed for snow type variable
num_colors <- length(unique(data_merge$ATTRIBUTE_SnowType))

#choose color blind pallete
color_palette <- brewer.pal(num_colors, "Set1")

#make scree plot 
res.pca <- data_merge  %>% dplyr::select(-ATTRIBUTE_Snow_Type) %>% 
  prcomp(scale = FALSE) #TRUE if scaling was not performed before
fviz_eig(res.pca)

data_merge$ATTRIBUTE_Snow_Type <- factor(data_merge$ATTRIBUTE_Snow_Type, levels = c("High Snowmaking", "Low Snowmaking","Natural Snowfall"))


#plot pca biplot
biplot <- fviz_pca_biplot(res.pca, repel = TRUE, 
                col.ind = data_merge$ATTRIBUTE_Snow_Type, #individual color
                label = "var",
                addEllipses = TRUE,
                select.var = list(contrib = 20), 
                col.var = "gray42", # Variables color
                palette = color_palette,
                legend.title = "Snow Type"
)
biplot

new_biplot <- biplot + 
  labs(title = "Metal and metalloid in Snow Types PCA Biplot Analysis", # Add a descriptive title
       x = "Principal Component 1 (30.7%)", # Replace with the actual explained variance
       y = "Principal Component 2 (26.8%)") +
  theme_minimal(base_size=12) + # Adjust the theme if needed
  theme(legend.position = "right",  # Adjust legend position if needed
  panel.grid = element_blank())
new_biplot

ggsave(file="pca_biplot_05232024_METALS_from_R.svg", plot=new_biplot, width=10, height=8)

# PERMANOVA
metadata <- data_merge[, 1]
metadata_df <- data.frame(metadata = metadata)
metabolites <- data_merge[, 2:ncol(data_merge)]
dist_metabolites <- vegdist(metabolites, method = "euclidean", na.rm = TRUE)
permanova_all <- adonis2(dist_metabolites ~ metadata_df$metadata, metadata_df, na.action = na.omit)

write.csv(permanova_all, "Permanova results 05232024.csv", row.names=TRUE)






#plot pca plot
pca <- fviz_pca_ind(res.pca, col.ind = data_merge$ATTRIBUTE_Snow_Type, 
                    addEllipses = TRUE, label = "none", 
                    legend.title= "Snow Type")
pca

data_merge$Snow_Type_with_ellipse <- ifelse(data_merge$ATTRIBUTE_Snow_Type == "Low Snowmaking", NA, data_merge$ATTRIBUTE_Snow_Type)

# Subset only rows where ellipses are desired
data_for_ellipses <- data_merge[!is.na(data_merge$Snow_Type_with_ellipse), ]

# Plot PCA with ellipses excluded for "Low Snowmaking"
pca_1 <- fviz_pca_ind(
  res.pca, 
  col.ind = data_merge$ATTRIBUTE_Snow_Type, # Individual color
  addEllipses = TRUE, 
  ellipse.group = data_for_ellipses$Snow_Type_with_ellipse, # Use filtered groups for ellipses
  label = "none", 
  legend.title = "Snow Type"
)

plot(pca_1)

new_pca <- pca_1 + 
  labs(title = "", # Add a descriptive title
       x = "Principal Component 1 (30.7%)", # Replace with the actual explained variance
       y = "Principal Component 2 (26.8%)") +
  theme_minimal(base_size=12) + # Adjust the theme if needed
  theme(legend.position = "right",  # Adjust legend position if needed
        panel.grid = element_blank())
new_pca

ggsave(file="pca_01222025.svg", plot=new_pca, width=10, height=8)

#output loadings plot
pc_loadings <- res.pca$rotation
pc_loadings <- as_tibble(pc_loadings, rownames = "feature")
write.csv(pc_loadings, "metals_CLR_loadings.csv",row.names = TRUE)

#plot driving features
fviz_pca_var(res.pca, select.var = list(contrib = 25), col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


