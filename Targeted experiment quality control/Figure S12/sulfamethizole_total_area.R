setwd("C:/Users/marqu/OneDrive - University of Denver/Desktop/DU Research/SMRW project/Targeted analysis of SMRW and NS/SMRW_PRM_Pharma_retry_06102024")
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(svglite)
library(extrafont)

#load csv
quant_file <- read.csv("SMRW_PRM_06112024_transitions_for_quant.csv")


# subset compounds by compound name 

sulfamethizole_d4_quant <- quant_file[quant_file$Molecule_Name =="Sulfamethizole-d4",]

# filter out QC 

sulfamethizole_d4_quant<- sulfamethizole_d4_quant %>% filter(Sample_Name != "QC_Pooled")

# filter NA values 

sulfamethizole_quant <- subset(sulfamethizole_d4_quant, Sample_Name != "Other")




# IS stability 
sulfamethizole_quant$Sample_Name <- factor(sulfamethizole_quant$Sample_Name, levels = c("Calibration_Blank", "0.1_ppb_calibrant",
                                                                                         "0.5_ppb_calibrant","1ppb_calibrant", 
                                                                                         "2.5ppb_calibrant","5ppb_calibrant",
                                                                                         "10ppb_calibrant","25ppb_calibrant", 
                                                                                         "50ppb_calibrant", "100ppb_calibrant", 
                                                                                         "Extraction_Blank", "Natural_Snowfall", "Low_Snowmaking",
                                                                                         "High_Snowmaking_Mar_2022", "High_Snowmaking_Dec_2022", 
                                                                                         "High_Snowmaking_Dil_Dec_2022"))


# Plot carbamazepine d10 total area 
sulfamethizole_quant_2 <- sulfamethizole_quant %>%
  tail(40) %>%
  ggplot( aes(x= Sample_Name, y= Total_Area)) +
  geom_line(linetype = "dashed") +
  geom_point(size= 3) +
  theme_ipsum() + 
  theme_classic() + 
  labs(y= "Peak Area (AU)", x= expression(paste("")))+
  theme(axis.title.x = element_text(family = "sans", size = 16, hjust = 0.5),  # Centering x-axis title
        axis.title.y = element_text(family = "sans", size = 16, hjust = 0.5),  # Centering y-axis title
        axis.text.x = element_text(family = "sans", size = 16, color = "black", angle = 45, hjust =1), 
        axis.text.y = element_text(family = "sans", size = 16, color= "black"),
        axis.line=element_line(size =1),
        axis.ticks=element_line(size=1), 
        axis.ticks.length=unit(.25, "cm"))

sulfamethizole_quant_2

ggsave(file="sulfamethizole_d6_area.svg", plot=sulfamethizole_quant_2, width =8, height = 6 )
