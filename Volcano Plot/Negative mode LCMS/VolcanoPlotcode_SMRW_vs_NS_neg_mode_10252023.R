setwd("C:/Aron_Lab_Data/SMRW_vs_NS_project/Volcano_plot_SMRW_vs_NS_neg_mode")
setwd("C:/Users/Owner/OneDrive - University of Denver/Desktop/Aron_Lab_Data/SMRW_vs_NS_project/Volcano_plot_SMRW_vs_NS_neg_mode")

#read in libraries
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(broom)
library(AICcmodavg)
library(dplyr)

# Read the data into R
dset<-read.csv("mzmine_2_1016_2023_Normalised_Quant_table_neg_seperated.csv", header = TRUE)
data<- as.matrix(dset[,4:19])
featureIDs <- dset[,1]
# Check the loaded dataset
# Dimension of the dataset
dim(data) 
head(data)

# Separate the two conditions into two smaller data frames
SMRW = data[,1:9]
NS = data[,10:15]
QC = data[,16]

# Compute the means of the samples of each condition
SMRW.mean = apply(SMRW, 1, mean)
NS.mean = apply(NS, 1, mean)


# Just get the maximum of all the means for plotting
limit = max(SMRW.mean, NS.mean)
means <- cbind.data.frame(SMRW.mean, NS.mean)

###############################################################
######### SMRW_High VS NS ANALYSIS FROM THIS POINT ONWARD ###########
DEres.df <- as.data.frame(data)
DEres.df <- cbind.data.frame(featureIDs, DEres.df)
#compute fold
fold = SMRW.mean / NS.mean
log2fold <- log2(fold)
#add to results dataframe
DEres.df$logFC <- log2fold

# Generate histogram of the fold differences
foldhist <- ggplot(DEres.df, aes(x=logFC)) +
  geom_histogram( binwidth=0.01, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  ggtitle("Fold Change Distribution") +
  labs(x= "Log fold change in abundance between groups", y = "Count") +
  theme_minimal() 
foldhist

#compute stat sig - t-test
pvalue = NULL # Empty list for the p-values
tstat = NULL # Empty list of the t test statistics
for(i in 1 : nrow(data)) { # For each gene : 
  x = SMRW[i,] # SMRW_High of gene number i
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

write.csv(DEres.df, "All_SMRW_vs_NS_Averaged_technical_replicates_volcano_plot_with_significance_12172024.csv", row.names = FALSE)


#### plotting volcano with no annotations

data <-read.csv("All_SMRW_vs_NS_Averaged_technical_replicates_volcano_plot_with_significance_12172024.csv", header = TRUE)
quant <- read.csv("mzmine_2_1016_2023_Normalised_Quant_table_neg_with_seperated_feature_mass_retentiontime_columns.csv", header = TRUE)

data <- data %>%
  left_join(select(quant, FeatureID, Mass), by = c("featureIDs" = "FeatureID"))

data <- data %>%
  mutate(Mass = ifelse(logFC < 5 & logpval < 2, NA, Mass))

data <- data %>%
  mutate(featureIDs = ifelse(logFC < 5 & logpval < 2, NA, featureIDs))


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
    aes(label = ifelse(Sample_Type == "Significant in SMRW", Mass, NA)),
    hjust = 0, vjust = 0, size = 3
  )


library(svglite)
require("ggplot2")
library(plotly)
ggplotly(myvolcano)

ggsave(file="more abundant features in SMRW 12172024 volcano plot featureIDs.svg", plot=myvolcano, width=10, height=8)


## volcano plot with no labels mass 
# Biostatsquid theme
theme_set(theme_classic(base_size = 20) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))

myvolcano2 <- ggplot(data = data, aes(x = logFC, y = logpval, col = Sample_Type)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 1, alpha = 0.4) +
  scale_color_manual(values = c("grey","#00AFBB", "#bb0c00"), # to set the colours of our variable
                     labels = c("Not Significant", "Significant in NS", "Significant in SMRW")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 5), xlim = c(-15, 15)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Differential Abundance') +
  scale_x_continuous(breaks = seq(-15, 15, 5)) + # to customise the breaks in the x axis
  ggtitle('Difference in Pharmaceuticals in Snow Type')+
  theme(plot.title =element_text(size=20), axis.title=element_text(size=15))


ggplotly(myvolcano2)

ggsave(file="more abundant features in SMRW 12172024 volcano plot.svg", plot=myvolcano2, width=10, height=8)

# Volcano Plots:
# put the biological significance (fold changes)
# and statistical significance (p-value) in one plot
# Generate the volcano plot
volcanoPlot <- function(results.df, fold_cutoff=0.5,
                        pvalue_cutoff=0.05, title = 
                          "Volcano Plot of Features"){
  # Inputs: ##############################################
  
  # results.df: dataframe containing columns: 
  #             1) "featureIDs": feature IDs
  #             2) "logFC": foldchange values from DE analysis
  #             3) "logpval": log-transformed p-values from DE analysis
  # fold cutoff: significance threshold to render visually on plot;
  #               denotes fold difference between mutant and wildtype
  #               also referred to as "biologial signal"
  # p_value_cutoff: significance threshold to render visually on plot;
  #                 denotes statistical significance
  ########################################################
  # create factor denoting differential expression status
  stats.df <- results.df %>% mutate(diffexpressed = 
                                      ifelse(logFC < -(fold_cutoff) & 
                                               logpval >= -log10(pvalue_cutoff),
                                             "Significant, NS", 
                                             ifelse(logFC > fold_cutoff & 
                                                      logpval >= -log10(pvalue_cutoff), 
                                                    "Significant, SMRW", "No difference in snow type")))
  
  stats.df$diffexpressed <- as.factor(stats.df$diffexpressed)
  # Base plot
  plot <- ggplot(stats.df, aes(x = logFC, y = logpval, col = diffexpressed, text = featureIDs)) +
    geom_point() +
    ggtitle(title) + 
    ylab("-log10(pvalue)") +
    geom_vline(xintercept=c(-fold_cutoff, fold_cutoff), col="red") +
    geom_hline(yintercept=-log10(pvalue_cutoff), col="darkturquoise") +
    scale_color_manual(values=c("black", "blue", "red"), 
                       name = "Differential feature abundance") +
    theme_minimal()
  return(plot)
  
}

VolcanoPlot
library("plotly")
# Volcano: put the biological significance (fold-change)
# and statistical significance (p-value) in one plot
plot <- volcanoPlot(stats.df, fold_cutoff=0.5, pvalue_cutoff=0.05)
ggplotly(plot)

install.packages('svglite')
library(svglite)
ggsave(file="SMRW_vs_NS_neg_volcano_12172024.svg", plot=plot, width=10, height=8)
write.csv(DEres.df, "SMRW_vs_NS_volcano_t-test_neg_12172024.csv", row.names=FALSE)


############################################################################
#####box plot with wilcoxon and krustal-wallis stats for feature m/z 797####
library("ggpubr")
library("datatable")
library("dplyr")
library("tidyverse")
library("readr")
my_data2 <- read.csv("~/Dropbox/2023-Lanthanide_Nature/R-Code_Positive/Volcano_Inputs/Normalized_Imp_Quant_799.csv")
##note: this was created from copying feature ID 10001_797.404_8.868 into the new table above##
p <- ggboxplot(my_data2,x="group", y="norm_peak_area", color="group", add = "jitter", shape = "group", yscale="log10")
p

my_comparisons <- list( c("A7", "KO"), c("KO", "WT"), c("A7", "WT"), c("WT", "WT_limited") )
fig <- p + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
  stat_compare_means(label.y = .0012)
fig
BH_corr <- pairwise.wilcox.test(my_data2$norm_peak_area, my_data2$group, p.adjust.method = "BH", paired = FALSE)
m1m2sup <- BH_corr[["p.value"]]
p_m1m2sup <- c(m1m2sup[1,1], m1m2sup[2,1], m1m2sup[2,2])

ggsave(file="~/Dropbox/2023-Lanthanide_Nature/R-Code_Positive/Volcano_Inputs/wilcox_boxplot.svg", plot=fig, width=10, height=8)

############################################################################
#####box plot with anova and stats for feature m/z 797####

library("ggpubr")
library("data.table")
my_data2 <- read.csv("~/Dropbox/2023-Lanthanide_Nature/R-Code_Positive/Volcano_Inputs/Normalized_Imp_Quant_799.csv")
##note: this was created from copying feature ID 10001_797.404_8.868 into the new table above##
p <- ggboxplot(my_data2,x="group", y="norm_peak_area", color="group", add = "jitter", shape = "group", yscale="log10")
p

my_comparisons <- list( c("A7", "KO"), c("KO", "WT"), c("A7", "WT"), c("WT", "WT_limited") )
fig <- p + stat_compare_means(comparisons = my_comparisons, method = "t.test") + stat_compare_means(method = "anova", label.y = .0012)
fig

ggsave(file="~/Dropbox/2023-Lanthanide_Nature/R-Code_Positive/Volcano_Inputs/anova_boxplot.svg", plot=fig, width=10, height=8)
