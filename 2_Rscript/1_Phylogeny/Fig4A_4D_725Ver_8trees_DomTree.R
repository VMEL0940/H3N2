## Package ##
library(tidyverse)
library(ggtree)
library(treeio)

## Set working directory ##
setwd("~/H3N2/")

## Data import ##
PB2_tree <- read.nexus("1_Data/1_Phylogeny/2_MCC/3_2nd_CrossValidation_725/01_PB2_MCC_Ordering_V1.tree")
PB1_tree <- read.nexus("1_Data/1_Phylogeny/2_MCC/3_2nd_CrossValidation_725/02_PB1_MCC_Ordering_V1.tree")
PA_tree <- read.nexus("1_Data/1_Phylogeny/2_MCC/3_2nd_CrossValidation_725/03_PA_MCC_Ordering_V1.tree")
HA_tree <- read.nexus("1_Data/1_Phylogeny/2_MCC/3_2nd_CrossValidation_725/04_HA_MCC_Ordering_V1.tree")
NP_tree <- read.nexus("1_Data/1_Phylogeny/2_MCC/3_2nd_CrossValidation_725/05_NP_MCC_Ordering_V1.tree")
NA_tree <- read.nexus("1_Data/1_Phylogeny/2_MCC/3_2nd_CrossValidation_725/06_NA_MCC_Ordering_V1.tree")
M_tree <- read.nexus("1_Data/1_Phylogeny/2_MCC/3_2nd_CrossValidation_725/07_M_MCC_Ordering_V1.tree")
NS_tree <- read.nexus("1_Data/1_Phylogeny/2_MCC/3_2nd_CrossValidation_725/08_NS_MCC_Ordering_V1.tree")

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


## 1. PB2 ##

## Data import ##
p <- ggtree(PB2_tree, mrsd = "2024-04-01") %<+% color.df +
  theme_tree2() +
  geom_tippoint(aes(color = vaccine_code_F), size = 1.8) +
  scale_color_manual(values = VacCoCol3, name = "Vaccine Period", guide = "none") +
  vexpand(.015)


ggsave("01_Veri725_PB2_noDom.svg", plot = p, device = "svg", 
       path = "3_Figures/Fig_4/Fig4A_Trees725",
       width = 300, height = 700, units = "mm", dpi = 720)


## 2. PB1 ##

## Data import ##
p <- ggtree(PB1_tree, mrsd="2024-04-01") %<+% color.df +
    theme_tree2() +
    geom_tippoint(aes(color = vaccine_code_F), size = 1.8) +
    scale_color_manual(values = VacCoCol3, name = "Vaccine Period", guide = "none") +
    vexpand(.015)

ggsave("02_Veri725_PB1_noDom.svg", plot = p, device = "svg", 
       path = "3_Figures/Fig_4/Fig4A_Trees725",
       width = 210, height = 300, units = "mm", dpi = 720)


## 3. PA ##

## Data import ##
p <- ggtree(PA_tree, mrsd="2024-04-01") %<+% color.df +
    theme_tree2() +
    geom_tippoint(aes(color = vaccine_code_F), size = 1.8) +
    scale_color_manual(values = VacCoCol3, name = "Vaccine Period", guide = "none") +
    vexpand(.015)

ggsave("03_Veri725_PA_noDom.svg", plot = p, device = "svg", 
       path = "3_Figures/Fig_4/Fig4A_Trees725",
       width = 210, height = 300, units = "mm", dpi = 720)


## 4. HA ##

## Data import ##
p <- ggtree(HA_tree, mrsd="2024-04-01") %<+% color.df +
    theme_tree2() +
    geom_tippoint(aes(color = vaccine_code_F), size = 1.8) +
    scale_color_manual(values = VacCoCol3, name = "Vaccine Period", guide = "none") +
    vexpand(.015)

ggsave("04_Veri725_HA_noDom.svg", plot = p, device = "svg", 
       path = "3_Figures/Fig_4/Fig4A_Trees725",
       width = 210, height = 300, units = "mm", dpi = 720)


## 5. NP ##

## Data import ##
p <- ggtree(NP_tree, mrsd="2024-04-01") %<+% color.df +
    theme_tree2() +
    geom_tippoint(aes(color = vaccine_code_F), size = 1.8) +
    scale_color_manual(values = VacCoCol3, name = "Vaccine Period", guide = "none") +
    vexpand(.015)

ggsave("05_Veri725_NP_noDom.svg", plot = p, device = "svg", 
       path = "3_Figures/Fig_4/Fig4A_Trees725",
       width = 210, height = 300, units = "mm", dpi = 720)



## 6. NA ##

## Data import ##
p <- ggtree(NA_tree, mrsd="2024-04-01") %<+% color.df +
    theme_tree2() +
    geom_tippoint(aes(color = vaccine_code_F), size = 1.8) +
    scale_color_manual(values = VacCoCol3, name = "Vaccine Period", guide = "none") +
    vexpand(.015)

ggsave("06_Veri725_NA_noDom.svg", plot = p, device = "svg", 
       path = "3_Figures/Fig_4/Fig4A_Trees725",
       width = 210, height = 300, units = "mm", dpi = 720)


## 7. M ##

## Data import ##
p <- ggtree(M_tree, mrsd="2024-04-01") %<+% color.df +
    theme_tree2() +
    geom_tippoint(aes(color = vaccine_code_F), size = 1.8) +
    scale_color_manual(values = VacCoCol3, name = "Vaccine Period", guide = "none") +
    vexpand(.015)

ggsave("07_Veri725_M_noDom.svg", plot = p, device = "svg", 
       path = "3_Figures/Fig_4/Fig4A_Trees725",
       width = 210, height = 300, units = "mm", dpi = 720)


## 8. NS ##

## Data import ##
p <- ggtree(NS_tree, mrsd="2024-04-01") %<+% color.df +
    theme_tree2() +
    geom_tippoint(aes(color = vaccine_code_F), size = 1.8) +
    scale_color_manual(values = VacCoCol3, name = "Vaccine Period", guide = "none") +
    vexpand(.015)

ggsave("08_Veri725_NS_noDom.svg", plot = p, device = "svg", 
       path = "3_Figures/Fig_4/Fig4A_Trees725",
       width = 210, height = 300, units = "mm", dpi = 720)




## 4.1 HA_Dom ##



## Data import ##
p <- ggtree(HA_tree, mrsd="2024-04-01") %<+% color.df +
    theme_tree2() +
    geom_tippoint(aes(fill = vaccine_code_F, color = Dom_F), size = 1.5, shape = 21, stroke = 0.65) +
    scale_fill_manual(values = VacCoCol3, guide = "none") +
    scale_color_manual(values = Dom, name = "Dominance") +
    vexpand(.015)

ggsave("Dominance_Tree_725.svg", plot = p, device = "svg",
       path = "3_Figures/Fig_4/Fig4D_DomTree725",
       width = 350, height = 700, units = "mm", dpi = 720)
