## Package ##
library(tidyverse)
library(treeio)
library(ggtree)

## Set working directory ##
setwd("~/H3N2/")

## Data import ##
tree <- read.nexus("1_Data/1_Phylogeny/2_MCC/1_Before2020_601/Work3_HA_MCC_edited_OctPatch.trees")
ggtree(tree)

## Data import ##
All_HA_725 <- read.csv("1_Data/10_Final_MetaData_Submission/H3N2_Testdata_724strains_1_Metadata.csv", header = T, na.strings = "")

colnames(All_HA_725)[2] <- "taxa"

## Factor color set ##
All_HA_725$vaccine_code_F <- factor(All_HA_725$vaccine_code, 
                            levels = c("Dar21",
                                       "Cam20",
                                       "HK19",
                                       "Kan17",
                                       "HK15",
                                       "Swtz13",
                                       "Vic11",
                                       "Prth09",
                                       "Bris07",
                                       "Wis05",
                                       "Cal04",
                                       "Fuj02",
                                       "Mos99"),
                            labels = c("13.Dar21",
                                       "12.Cam20",
                                       "11.HK19",
                                       "10.Kan17",
                                       "9.HK15",
                                       "8.Swtz13",
                                       "7.Vic11",
                                       "6.Prth09",
                                       "5.Bris07",
                                       "4.Wis05",
                                       "3.Cal04",
                                       "2.Fuj02",
                                       "1.Mos99"))

trend_colors_2024 <- c("#ffbe98", "#e2552c", "#515b87",  "#f8d974", "#42bdcb", 
                       "#00a170", "#a23c26", "#9eb4d3", "#748c69", "#ff9687" )

VacCoCol3 <- c("#00a170","#f8d974","#A4D9E5","#e2552c",
               "#954495","#769F74","#092E50","#9E6D13","#800000","#6D9431","#3EA3C5","#EFA95F","#716172")


color.df <- All_HA_725


## Data import ##
ggtree(tree, mrsd="2024-04-01") %<+% color.df +
  theme_tree2() +
  geom_tippoint(aes(color = vaccine_code_F), size = 1.0) +
  scale_color_manual(values = VacCoCol3, name = "Vaccine Period", guide = "none") 


## Ad dominance 
All_HA_725$Dom_F <- factor(All_HA_725$Dom, 
                         levels = c(2,0,1),
                         labels = c("No Prediction","Extinct","Dominant"))

color.df <- All_HA_725[,c(2,9,10)]

Dom <- c("#CCCCCC","#333333","#96fd29")

## Data import ##
p <- ggtree(tree, mrsd="2024-04-01") %<+% color.df +
    theme_tree2() +
    geom_tippoint(aes(fill = vaccine_code_F, color = Dom_F), size = 1.8, shape = 21, stroke = 0.65) +
    scale_fill_manual(values = VacCoCol3, name = "Vaccine Period") +
    scale_color_manual(values = Dom, name = "Dominance") 


ggsave("Dominance_Tree_601.svg", plot = p, device = "svg",
       path = "3_Figures/Fig_3/Fig3A_DomTree601",
       width = 350, height = 700, units = "mm", dpi = 720)

