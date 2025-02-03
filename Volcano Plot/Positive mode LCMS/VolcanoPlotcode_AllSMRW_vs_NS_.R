setwd("C:/Users/marqu/OneDrive - University of Denver/Desktop/Aron_Lab_Data/SMRW_vs_NS_project/Volcano_Plot_Averaged_technical_replicates/averaged_features_allSMRW_NS")
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

#compute stat sig - t-test
pvalue = NULL # Empty list for the p-values
tstat = NULL # Empty list of the t test statistics
for(i in 1 : nrow(data)) { # For each gene : 
  x = SMRW[i,] # SMRW of gene number i
  y = NS[i,] # NS of gene number i
  
  # Compute t-test between the two conditions
  t = t.test(x, y, na.rm=TRUE, paired=FALSE, exact=FALSE)
  
  # Put the current p-value in the pvalues list
  pvalue[i] = t$p.value
  # Put the current t-statistic in the tstats list
  tstat[i] = t$statistic
}
# add results to results dataframe
DEres.df$pvalue <- pvalue

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

write.csv(DEres.df, "All_SMRW_vs_NS_Averaged_technical_replicates_volcano_plot_with_significance_01252024.csv", row.names = FALSE)



setwd("C:/Aron_Lab_Data/SMRW_vs_NS_project/GNPS/pos_updated_FBMNannotations_SMRW_vs_10022023")
mydf <- read.csv("All_SMRW_vs_NS_Averaged_technical_replicates_volcano_plot_with_significance_01252024.csv", 
                 sep = ",", quote="\"")
mylookup <- read.csv("pos_updated_gnps_10022023_FBMN_with_annotations_SMRW_vs_NS.csv", 
                     quote="\"", sep = "," )
head(mydf)
head(mylookup)
joined_df <- merge(mydf, mylookup, by.x = "featureIDs", 
                   by.y = "X.Scan.", all.x = TRUE, all.y = FALSE)
write.csv(joined_df, "joined_df.csv", row.names = FALSE)


######## Plotting the annotated compound names with the t test pvalue and log FC change of High SMRW_vs_NS######3

data <-read.csv("AllSMRW_average_technical_replicates_replicates_plot_with_annotations_Of_interest_01252024.csv", header = TRUE)


# Loading relevant libraries 
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(broom)
library(AICcmodavg)
library(dplyr)

library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(ggrepel) # for nice annotations


# Biostatsquid theme
theme_set(theme_classic(base_size = 20) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))

myvolcano <- ggplot(data = data, aes(x = logFC, y = -log10(pvalue), col = Sample_Type, label = Compound_Name)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 1, alpha = 0.2) +
  geom_text(size=3,  col = "black", nudge_y= 0.01) +
  scale_color_manual(values = c("grey","#00AFBB", "#bb0c00"), # to set the colours of our variable
                     labels = c("Not Significant", "Significant in NS", "Significant in SMRW")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 5), xlim = c(-15, 15)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Differential Abundance') +
  scale_x_continuous(breaks = seq(-15, 15, 5)) + # to customise the breaks in the x axis
  ggtitle('Difference in Pharmaceuticals in Snow Type')+
  theme(plot.title =element_text(size=20), axis.title=element_text(size=15))



library(svglite)
require("ggplot2")
library(plotly)
ggplotly(myvolcano)


ggsave(file="All_SMRW_NS_averaged_technical_replicate_annotates_with_alpha0.2_05232024.svg", plot=myvolcano, width=10, height=8)

######## Plot volcano plot with no annotations ########


data <-read.csv("All_SMRW_vs_NS_Averaged_technical_replicates_volcano_plot_with_significance_01252024.csv", header = TRUE)
quant <- read.csv("average_technical_replicates_10292023_pos_Normalised_Quant.csv", header = TRUE)

data <- data %>%
  left_join(select(quant, FeatureID, mz), by = c("featureIDs" = "FeatureID"))

data <- data %>%
  mutate(mz = ifelse(logFC < 12 & logpval < 3, NA, mz))

data <- data %>%
  mutate(featureIDs = ifelse(logFC < 12 & logpval < 3, NA, featureIDs))

# Biostatsquid theme
theme_set(theme_classic(base_size = 20) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))

myvolcano <- ggplot(data = data, aes(x = logFC, y = logpval, col = Sample_Type)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 1, alpha = 0.2) +
  scale_color_manual(values = c("grey","#00AFBB", "#bb0c00"), # to set the colours of our variable
                     labels = c("Not Significant", "Significant in NS", "Significant in SMRW")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 5), xlim = c(-15, 15)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Differential Abundance') +
  scale_x_continuous(breaks = seq(-15, 15, 5)) + # to customise the breaks in the x axis
  ggtitle('Difference in Pharmaceuticals in Snow Type')+
  theme(plot.title =element_text(size=20), axis.title=element_text(size=15)) +
  geom_text(
    aes(label = ifelse(Sample_Type == "Significant in SMRW", featureIDs, NA)),
    hjust = 0, vjust = 0, size = 3
  )


library(svglite)
require("ggplot2")
library(plotly)
ggplotly(myvolcano)

ggsave(file="more abundant features in SMRW 11262024 volcano plot featureIDs.svg", plot=myvolcano, width=10, height=8)
