## Package ##
library(tidyverse)
library(ggtree)
library(treeio)

## Set working directory ##
setwd("~/H3N2/")

## Data import ##
PB2_tree <- read.nexus("1_Data/1_Phylogeny/2_MCC/1_Before2020_601/Work3_PB2_MCC_edited_OctPatch.trees")
PB1_tree <- read.nexus("1_Data/1_Phylogeny/2_MCC/1_Before2020_601/Work3_PB1_MCC_edited_OctPatch.trees")
PA_tree <- read.nexus("1_Data/1_Phylogeny/2_MCC/1_Before2020_601/Work3_PA_MCC_edited_OctPatch.trees")
HA_tree <- read.nexus("1_Data/1_Phylogeny/2_MCC/1_Before2020_601/Work3_HA_MCC_edited_OctPatch.trees")
NP_tree <- read.nexus("1_Data/1_Phylogeny/2_MCC/1_Before2020_601/Work3_NP_MCC_edited_OctPatch.trees")
NA_tree <- read.nexus("1_Data/1_Phylogeny/2_MCC/1_Before2020_601/Work3_NA_MCC_edited_OctPatch.trees")
M_tree <- read.nexus("1_Data/1_Phylogeny/2_MCC/1_Before2020_601/Work3_M_MCC_edited_OctPatch.trees")
NS_tree <- read.nexus("1_Data/1_Phylogeny/2_MCC/1_Before2020_601/Work3_NS_MCC_edited_OctPatch.trees")

## Data import ##
All_HA_601 <- read.csv("1_Data/10_Final_MetaData_Submission/H3N2_Traindata_600strains_1_Metadata.csv", header = T, na.strings = "")

colnames(All_HA_601)[2] <- "taxa"

## Factor color set ##
All_HA_601$vaccine_code_F <- factor(All_HA_601$vaccine_code, 
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


VacCoCol3 <- c("#954495","#769F74","#092E50","#9E6D13","#800000","#6D9431","#3EA3C5","#EFA95F","#716172")


color.df <- All_HA_601[,c(2,8)]


## Ad dominance 
All_HA_601$Dom_F <- factor(All_HA_601$Dom, 
                           levels = c(2,0,1),
                           labels = c("No Prediction","Extinct","Dominant"))

color.df <- All_HA_601[,c(2,8,9)]

Dom <- c("#CCCCCC","#333333","#96fd29")


## 1. PB2 ##

## Data import ##
p <- ggtree(PB2_tree, mrsd = "2024-04-01") %<+% color.df +
  theme_tree2() +
  geom_tippoint(aes(color = vaccine_code_F), size = 1.8) +
  scale_color_manual(values = VacCoCol3, name = "Vaccine Period", guide = "none") +
  vexpand(.015)


ggsave("01_Veri601_PB2_noDom.svg", plot = p, device = "svg", 
       path = "3_Figures/Fig_1/Fig1A_trees",
       width = 210, height = 300, units = "mm", dpi = 720)


## 2. PB1 ##

## Data import ##
p <- ggtree(PB1_tree, mrsd="2024-04-01") %<+% color.df +
    theme_tree2() +
    geom_tippoint(aes(color = vaccine_code_F), size = 1.8) +
    scale_color_manual(values = VacCoCol3, name = "Vaccine Period", guide = "none") +
    vexpand(.015)

ggsave("02_Veri601_PB1_noDom.svg", plot = p, device = "svg", 
       path = "3_Figures/Fig_1/Fig1A_trees",
       width = 210, height = 300, units = "mm", dpi = 720)


## 3. PA ##

## Data import ##
p <- ggtree(PA_tree, mrsd="2024-04-01") %<+% color.df +
    theme_tree2() +
    geom_tippoint(aes(color = vaccine_code_F), size = 1.8) +
    scale_color_manual(values = VacCoCol3, name = "Vaccine Period", guide = "none") +
    vexpand(.015)

ggsave("03_Veri601_PA_noDom.svg", plot = p, device = "svg", 
       path = "3_Figures/Fig_1/Fig1A_trees",
       width = 210, height = 300, units = "mm", dpi = 720)


## 4. HA ##

## Data import ##
p <- ggtree(HA_tree, mrsd="2024-04-01") %<+% color.df +
    theme_tree2() +
    geom_tippoint(aes(color = vaccine_code_F), size = 1.8) +
    scale_color_manual(values = VacCoCol3, name = "Vaccine Period", guide = "none") +
    vexpand(.015)

ggsave("04_Veri601_HA_noDom.svg", plot = p, device = "svg", 
       path = "3_Figures/Fig_1/Fig1A_trees",
       width = 210, height = 300, units = "mm", dpi = 720)


## 5. NP ##

## Data import ##
p <- ggtree(NP_tree, mrsd="2024-04-01") %<+% color.df +
    theme_tree2() +
    geom_tippoint(aes(color = vaccine_code_F), size = 1.8) +
    scale_color_manual(values = VacCoCol3, name = "Vaccine Period", guide = "none") +
    vexpand(.015)

ggsave("05_Veri601_NP_noDom.svg", plot = p, device = "svg", 
       path = "3_Figures/Fig_1/Fig1A_trees",
       width = 210, height = 300, units = "mm", dpi = 720)



## 6. NA ##

## Data import ##
p <- ggtree(NA_tree, mrsd="2024-04-01") %<+% color.df +
    theme_tree2() +
    geom_tippoint(aes(color = vaccine_code_F), size = 1.8) +
    scale_color_manual(values = VacCoCol3, name = "Vaccine Period", guide = "none") +
    vexpand(.015)

ggsave("06_Veri601_NA_noDom.svg", plot = p, device = "svg", 
       path = "3_Figures/Fig_1/Fig1A_trees",
       width = 210, height = 300, units = "mm", dpi = 720)


## 7. M ##

## Data import ##
p <- ggtree(M_tree, mrsd="2024-04-01") %<+% color.df +
    theme_tree2() +
    geom_tippoint(aes(color = vaccine_code_F), size = 1.8) +
    scale_color_manual(values = VacCoCol3, name = "Vaccine Period", guide = "none") +
    vexpand(.015)

ggsave("07_Veri601_M_noDom.svg", plot = p, device = "svg", 
       path = "3_Figures/Fig_1/Fig1A_trees",
       width = 210, height = 300, units = "mm", dpi = 720)


## 8. NS ##

## Data import ##
p <- ggtree(NS_tree, mrsd="2024-04-01") %<+% color.df +
    theme_tree2() +
    geom_tippoint(aes(color = vaccine_code_F), size = 1.8) +
    scale_color_manual(values = VacCoCol3, name = "Vaccine Period", guide = "none") +
    vexpand(.015)

ggsave("08_Veri601_NS_noDom.svg", plot = p, device = "svg", 
       path = "3_Figures/Fig_1/Fig1A_trees",
       width = 210, height = 300, units = "mm", dpi = 720)
