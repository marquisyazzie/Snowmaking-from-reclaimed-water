setwd("C:/Users/Marquis.Yazzie/OneDrive - University of Denver/Desktop/Aron_Lab_Data/SMRW_vs_NS_project/Volcano_Plot_Averaged_technical_replicates/averaged_features_allSMRW_NS")

#read in libraries
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(broom)
library(AICcmodavg)
library(dplyr)

# Read the data into R
dset<-read.csv("average_technical_replicates_10292023_pos_Normalised_Quant.csv", header = TRUE)
data<- as.matrix(dset[,4:18])
featureIDs <- dset[,1]
# Check the loaded dataset
# Dimension of the dataset
dim(data) 
head(data)

# Separate the two conditions into two smaller data frames
SMRW = data[,1:9]
NS = data[,10:15]

# Compute the means of the samples of each condition
SMRW.mean = apply(SMRW, 1, mean)
NS.mean = apply(NS, 1, mean)

# Just get the maximum of all the means for plotting
limit = max(SMRW.mean, NS.mean)
means <- cbind.data.frame(SMRW.mean, NS.mean)

###############################################################
######### SMRW VS NS ANALYSIS FROM THIS POINT ONWARD ###########
DEres.df <- as.data.frame(data)
DEres.df <- cbind.data.frame(featureIDs, DEres.df)
#compute fold
fold = SMRW.mean / NS.mean
log2fold <- log2(fold)
#add to results dataframe
DEres.df$logFC <- log2fold

pvalue = NULL # Empty list for the p-values

for (i in 1:nrow(data)) { # For each gene : 
  x = SMRW[i,] # SMRW of gene number i
  y = NS[i,]   # NS of gene number i
  
  # Compute Mann-Whitney U test (Wilcoxon rank sum test)
  w = wilcox.test(x, y, na.rm = TRUE, exact = FALSE, paired = FALSE)
  
  # Store results
  pvalue[i] = w$p.value
}

# Add results to results dataframe
DEres.df$pvalue <- pvalue

# FDR correction
DEres.df$FDR <- p.adjust(DEres.df$pvalue, method = "fdr")


# Generate -log10(pvalues)
pvallogged <- -log10(pvalue)
# add results to results dataframe
DEres.df$logpval <- pvallogged
# Generate histogram of -log10(pvalues)
pval_hist <- ggplot(DEres.df, aes(x=logpval)) +
  geom_histogram( binwidth=0.05, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  ggtitle("-log10(pvalue) Distribution") +
  labs(x= "Significance levels occuring in D.E. results", y = "Count") +
  xlim(0,4) +
  theme_minimal() 
pval_hist



# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)<br /><br /><br />
DEres.df$Sample_Type <- "Not Significant"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "Significant in SMRW"
DEres.df$Sample_Type[DEres.df$logFC > 0.6 & DEres.df$pvalue < 0.05] <- "Significant in SMRW"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "Significant in NS"
DEres.df$Sample_Type[DEres.df$logFC < -0.6 & DEres.df$pvalue < 0.05] <- "Significant in NS"

write.csv(DEres.df, "All_SMRW_vs_NS_volcano_plot_MAnn_whitney_09222025.csv", row.names = FALSE)

library(svglite)
library(plotly)

# Biostatsquid theme
theme_set(theme_classic(base_size = 20) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))

myvolcano <- ggplot(data = DEres.df, 
                    aes(x = logFC, y = -log10(pvalue), 
                        col = Sample_Type, label = featureIDs)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 1, alpha = 0.2) +
  scale_color_manual(
    values = c("Not Significant" = "grey", "Significant in NS" = "#00AFBB", "Significant in SMRW" = "#bb0c00"),
    labels = c("Not Significant" = "Not Significant", 
               "Significant in NS" = "Significant in Natural Snowfall", 
               "Significant in SMRW" = "Significant in Snowmaking")
  ) +
  coord_cartesian(ylim = c(0, 5), xlim = c(-15, 15)) +
  labs(color = 'Differential Abundance') +
  scale_x_continuous(breaks = seq(-15, 15, 5)) +
  ggtitle('Difference in Pharmaceuticals in Snow Type') +
  theme(plot.title = element_text(size=20), 
        axis.title = element_text(size=15))

ggplotly(myvolcano)


ggsave(file="All_SMRW_NS_s_with_alpha0.2_Mann Whitney.png", plot=myvolcano, width=8, height=6)
