### Influ H3N2 GenDist Proj ###
## Package ##
library(tidyverse)
library(RColorBrewer)


## Set working directory ##
setwd("/01_GenDisFlu/")

## Data import ##
tMRCA <- read.csv("Data/tMRCA/tMRCA.csv", header = T, na.strings = "")

tMRCA$Gene <- factor(tMRCA$Gene, levels = c("PB2","PB1","PA","HA","NP", "NA", "M", "NS"))

GeneCol <- c("#C780E8", "#9D94FF","#59ADF6","#08CAD1","#42D6A4","#F8F38D","#FFB480","#FF6961")

tMRCA$VGroup <- factor(tMRCA$Group, levels = c(1,2,3,4,5,6,7,8,9),
                               labels = c("1.Mos99","2.Fuj02","3.Cal04",
                                          "4.Wis05","5.Bris07","6.Prth09",
                                          "7.Vic11", "8.Swtz13", "9.HK15"))

VacCoCol <- c("#716172", "#EFA95F","#3EA3C5","#6D9431","#800000","#9E6D13","#092E50","#769F74","#954495")


## Visualization


jpeg(filename = "/01_GenDisFlu/Fig/Desc/tMRCA/HgR_tMRCA.jpeg", width = 280, height = 180, units = "mm", res = 320)

# Solo plot
tMRCA %>% 
  ggplot(aes(x=tMRCA_Hg, y = reorder(Gene, desc(Gene)), color = VGroup)) + 
  geom_errorbar(aes(x= tMRCA_Hg, y= reorder(Gene, desc(Gene)), xmin=tMRCA_Hg_L95, xmax=tMRCA_Hg_H95), color = "black", width=0.3, size = 0.7) +
  geom_point(aes(x=tMRCA_Hg, y = reorder(Gene, desc(Gene)), fill = VGroup),  size = 4, color = "black", shape =21) +
  ylab("Gene Segment") +
  xlab("Time of most recent common ancestor (years)") +
  scale_fill_manual(values = VacCoCol, name = "Vaccine Strains") +
  scale_x_continuous(breaks = seq(1990, 2020, by = 2)) +
  theme_bw()  +
    theme(axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 9))

dev.off()

jpeg(filename = "/01_GenDisFlu/Fig/Desc/tMRCA/Gene_Vac_tMRCA.jpeg", width = 280, height = 220, units = "mm", res = 320)

# Facet_grid Figure
tMRCA %>% 
  ggplot(aes(x=tMRCA_Hg, y = reorder(Gene, desc(Gene)), color = Gene)) + 
  geom_errorbar(aes(x= tMRCA_Hg, y= reorder(Gene, desc(Gene)), xmin=tMRCA_Hg_L95, xmax=tMRCA_Hg_H95), color = "black", width=0.3, size = 0.7) +
  geom_point(aes(x=tMRCA_Hg, y = reorder(Gene, desc(Gene)), fill = Gene), size = 3, color = "black", shape =21) +
  ylab("") +
  scale_fill_manual(values = GeneCol, name = "Gene segments") +
  scale_x_continuous(breaks = seq(1990, 2020, by = 2)) +
  xlab("Time of most recent common ancestor (years)") +
  theme_bw()  +
  theme(axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 9)) +
  facet_grid(rows = vars(VGroup))

dev.off()


# Solo plot
tMRCA %>% 
  ggplot(aes(x=tMRCA_Hg, y = reorder(Gene, desc(Gene)), color = VGroup)) + 
  geom_errorbar(aes(x= tMRCA_Hg, y= reorder(Gene, desc(Gene)), xmin=tMRCA_Hg_L95, xmax=tMRCA_Hg_H95), color = "black", width=0.3, size = 0.7) +
  geom_jitter(aes(x=tMRCA_Hg, y = reorder(Gene, desc(Gene)), fill = VGroup),  width = 0, height = 0.13, size = 3.2, color = "black", shape =21) +
  ylab("Gene Segment") +
  scale_fill_manual(values = VacCoCol, name = "Vaccine Strains") +
  scale_x_continuous(breaks = seq(1990, 2020, by = 2)) +
  xlab("Time of most recent common ancestor (years)") +
  theme_bw()  +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9))


