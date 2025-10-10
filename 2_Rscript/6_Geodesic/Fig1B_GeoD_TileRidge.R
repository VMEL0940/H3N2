##### Library  #####
library(tidyverse)
library(corrplot)
library(ggridges)
library(RColorBrewer)

##### Import data ####
## Set working directory ##
setwd("~/H3N2/")

### import "FBQ" data 
GeoD_Summary <- read.csv("1_Data/6_Geodesic/1_50csv_new/All_Pair_Summary.csv", header = TRUE, na.string = "") # 

## Round mean value 
GeoD_Summary$DisMean <- round(GeoD_Summary$DisMean, 3)

## Subset of value
table <- GeoD_Summary[,c(2,3)]

## split by gene
table2 <- separate(data = table, col = GenePair, into = c("Gene1", "Gene2"), sep = "-")

table2$Gene1 <- factor(table2$Gene1, levels = c("PB2","PB1","PA","HA",
                                                "NP", "NA", "M", "NS"))

table2$Gene2 <- factor(table2$Gene2, levels = c("PB2","PB1","PA","HA",
                                                "NP", "NA", "M", "NS"))

table2 <- table2 %>% 
  arrange(Gene1, Gene2)

table2$Gene1[16] <- "PB2"
table2$Gene2[16] <- "HA"

table2$Gene1[17] <- "PB1"
table2$Gene2[17] <- "HA"

table2$Gene1[18] <- "PA"
table2$Gene2[18] <- "HA"

table2$Gene1[18] <- "PA"
table2$Gene2[18] <- "HA"

table2$Gene1[28] <- "M"
table2$Gene2[28] <- "NS"

table2 <- table2 %>% 
  arrange(Gene1, Gene2)

table2 <- table2 %>% 
  mutate(DisMean_s = scale(DisMean))

p <- table2 %>% 
  ggplot(aes(Gene1, Gene2, fill = DisMean_s,
             label = sprintf("%.3f", DisMean))) +
  geom_tile() +
  geom_text(col = "black") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_fill_gradient2(low = "#008080", mid = "#b0e0e6", high = "white")


ggsave("HA_GeoD.svg", plot = p, device = "svg", 
       path = "3_Figures/Fig_1/Fig1B_GeoD",
       width = 300, height = 300, units = "mm", dpi = 720)


################ For reference ################
### import "FBQ" data 
GeoD_All <- read.csv("1_Data/6_Geodesic/1_50csv_new/All_Pair_Results.csv", header = TRUE, na.string = "") # 

HA_GeoD <- GeoD_All %>% 
  filter(GenePair == "HA-M" | GenePair == "HA-NA" | GenePair == "HA-PB1" | GenePair == "HA-PB2" | 
           GenePair == "HA-PA" | GenePair == "HA-NS" | GenePair == "HA-NP")

HA_GeoD$GenePair <- factor(HA_GeoD$GenePair, levels = c("HA-PB2","HA-PB1","HA-PA",
                                                        "HA-NP", "HA-NA", "HA-M", "HA-NS"))


GeneCol <- c("#C780E8", "#9D94FF","#59ADF6","#42D6A4","#F8F38D","#FFB480","#FF6961")


p <- ggplot(HA_GeoD, aes(x = value, y = GenePair, fill = GenePair)) +
    geom_density_ridges() +
    scale_fill_manual(values = GeneCol, guide="none") +
    xlab("Tree-Tree Geodesic distance") +
    ylab("Gene Pair with the HA gene") +
    theme_bw() + 
    theme(legend.position="none")

