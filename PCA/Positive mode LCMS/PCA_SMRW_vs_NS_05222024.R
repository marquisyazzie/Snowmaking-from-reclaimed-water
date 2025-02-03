setwd("C:/Aron_Lab_Data/SMRW_vs_NS_project/PCA/Averaged_technical_replicate_pos_SMRW_v_NS")

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


metadata <- read.csv("average_technical_replicates_metadata_05222024_pos_mode.csv")
clr <- read.csv("average_technical_replicates_mzmine_3_10292023_pos_CLR_Scaled_positive.csv")

data_merge <- metadata %>% dplyr::select("filename", "ATTRIBUTE_Sample_Type") %>% 
  left_join(clr) %>% 
  column_to_rownames("filename")

# filter data in sample type column
data_merge <- data_merge %>% filter(ATTRIBUTE_Snow_Type != "SMRW_Low")
  

#specify the number of colors needed for snow type variable
num_colors <- length(unique(data_merge$ATTRIBUTE_SnowType))

#choose color blind pallete
color_palette <- brewer.pal(num_colors, "Set1")

#make scree plot 
res.pca <- data_merge  %>% dplyr::select(-ATTRIBUTE_Sample_Type) %>% 
  prcomp(scale = FALSE) #TRUE if scaling was not performed before
fviz_eig(res.pca)

#plot pca biplot
biplot <- fviz_pca_biplot(res.pca, repel = TRUE, 
                col.ind = data_merge$ATTRIBUTE_Sample_Type, #individual color
                label = "var",
                addEllipses = TRUE,
                select.var = list(contrib = 20), 
                col.var = "gray42", # Variables color
                palette = color_palette,
                legend.title = "Sample Type"
)
biplot

new_biplot <- biplot + 
  labs(title = "", # Add a descriptive title
       x = "Principal Component 1 (40.8%)", # Replace with the actual explained variance
       y = "Principal Component 2 (20%)") +
  theme_minimal(base_size=12) + # Adjust the theme if needed
  theme(legend.position = "right",  # Adjust legend position if needed
  panel.grid = element_blank())
new_biplot

ggsave(file="pca_biplot_05222024_organic_SMRW_vs_NS.svg", plot=new_biplot, width=10, height=8)

# PERMANOVA
metadata <- data_merge[, 1]
metadata_df <- data.frame(metadata = metadata)
metabolites <- data_merge[, 2:ncol(data_merge)]
dist_metabolites <- vegdist(metabolites, method = "euclidean", na.rm = TRUE)
permanova_all <- adonis2(dist_metabolites ~ metadata_df$metadata, metadata_df, na.action = na.omit)

write.csv(permanova_all, "Permanova results SMRW vs NS 05222024.csv", row.names=TRUE)






#plot pca plot
pca <- fviz_pca_ind(res.pca, col.ind = data_merge$ATTRIBUTE_Sample_Type, 
                    addEllipses = TRUE, label = "none", 
                    legend.title= "")
pca

new_pca <- pca + 
  labs(title = "", # Add a descriptive title
       x = "Principal Component 1 (40.8%)", # Replace with the actual explained variance
       y = "Principal Component 2 (20%)") +
  theme_minimal(base_size=12) + # Adjust the theme if needed
  theme(legend.position = "right",  # Adjust legend position if needed
        panel.grid = element_blank())
new_pca

ggsave(file="pca_SMRW_vs_NS_05222024.svg", plot=new_pca, width=10, height=8)

#output loadings plot
pc_loadings <- res.pca$rotation
pc_loadings <- as_tibble(pc_loadings, rownames = "feature")
write.csv(pc_loadings, "metals_CLR_loadings.csv",row.names = TRUE)

#plot driving features
fviz_pca_var(res.pca, select.var = list(contrib = 25), col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


