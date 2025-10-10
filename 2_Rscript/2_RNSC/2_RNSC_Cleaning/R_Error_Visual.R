### Influ H3N2 GenDist Proj ###
## Package ##
library(tidyverse)
library(factoextra)
library(FactoMineR)
library(ggpubr)

## Set working directory ##
setwd("/1_GenDisFlu/")

#### HA #####
## Data import ##
HA_G1 <- read.table("Data/Genetic Distances/HA distances/HA_G1_distances_2.txt", header = T)

## Pretest - Check Distribution of N value by Tree
HA_G1 %>% 
  filter(vaccineStrain == "EPI103320_A_Moscow_10_1999_NA_NA_NA_NA" & 
           compareStrain == "EPI103320_A_Moscow_10_1999_NA_NA_NA_NA") %>% 
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count", fill = "#99CCCC", color = "black") +
  xlab("Number of Non-synonymous between the same vaccine strain in the tree") +
  theme_bw()

## Comparing the Impact of Error
Fig_HAG1_N2 <- HA_G1 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count", fill = "#3333FF", color = "black") +
  xlab("Number of Non-synonymous between vaccine & wild strains") +
  theme_bw()

Fig_HAG1_N <- HA_G1 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count", fill = "#FF6666", color = "black") +
  xlab("Number of Non-synonymous between vaccine & wild strains") +
  theme_bw()

ggarrange(Fig_HAG1_N, Fig_HAG1_N2, 
          labels = c("Before Cleaning - HA Group 1","After Cleaning - HA Group 1"),
          nrow = 2)






## Pretest - Check Distribution of N value by Tree
NS_G1 %>% 
  filter(vaccineStrain == "EPI103320_A_Moscow_10_1999_NA_NA_NA_NA" & 
           compareStrain == "EPI103320_A_Moscow_10_1999_NA_NA_NA_NA") %>% 
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count", fill = "#99CCCC", color = "black") +
  xlab("Number of Non-synonymous between the same vaccine strain in the tree") +
  theme_bw()

## Comparing the Impact of Error
Fig_NSG1_N2 <- NS_G1 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count", fill = "#3333FF", color = "black") +
  xlab("Number of Non-synonymous between vaccine & wild strains") +
  theme_bw()

Fig_NSG1_N <- NS_G1 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count", fill = "#FF6666", color = "black") +
  xlab("Number of Non-synonymous between vaccine & wild strains") +
  theme_bw()

ggarrange(Fig_NSG1_N, Fig_NSG1_N2, 
          labels = c("Before Cleaning - NS Group 1","After Cleaning - NS Group 1"),
          nrow = 2)
