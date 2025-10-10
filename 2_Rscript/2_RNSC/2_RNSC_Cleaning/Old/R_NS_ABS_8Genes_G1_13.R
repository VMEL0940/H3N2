### Influ H3N2 GenDist Proj ###
## Package ##
library(tidyverse)

## Set working directory ##
setwd("01_GenDisFlu/")
#setwd("01_GenDisFlu/")


########### Import Data #################

## Data import -- RNSC ##
G1 <- read.delim("Data/Genetic Distances/Verification_Final_725/03_median/G1_median.tsv", header = T)
G2 <- read.delim("Data/Genetic Distances/Verification_Final_725/03_median/G2_median.tsv", header = T)
G3 <- read.delim("Data/Genetic Distances/Verification_Final_725/03_median/G3_median.tsv", header = T)
G4 <- read.delim("Data/Genetic Distances/Verification_Final_725/03_median/G4_median.tsv", header = T)
G5 <- read.delim("Data/Genetic Distances/Verification_Final_725/03_median/G5_median.tsv", header = T)
G6 <- read.delim("Data/Genetic Distances/Verification_Final_725/03_median/G6_median.tsv", header = T)
G7 <- read.delim("Data/Genetic Distances/Verification_Final_725/03_median/G7_median.tsv", header = T)
G8 <- read.delim("Data/Genetic Distances/Verification_Final_725/03_median/G8_median.tsv", header = T)
G9 <- read.delim("Data/Genetic Distances/Verification_Final_725/03_median/G9_median.tsv", header = T)
G10 <- read.delim("Data/Genetic Distances/Verification_Final_725/03_median/G10_median.tsv", header = T)
G11 <- read.delim("Data/Genetic Distances/Verification_Final_725/03_median/G11_median.tsv", header = T)
G12 <- read.delim("Data/Genetic Distances/Verification_Final_725/03_median/G12_median.tsv", header = T)
G13 <- read.delim("Data/Genetic Distances/Verification_Final_725/03_median/G13_median.tsv", header = T)

## Data import -- RNSC ##
Index_725 <- read.csv("Data/Veri_725_2024/total_index_725.csv", header = T, na.strings = "")

### Drop the vaccine strains 
G1 <- G1 %>% 
  filter(compareStrain != "EPI103320_A_Moscow_10_1999_NA_NA_NA_NA")
G2 <- G2 %>% 
  filter(compareStrain != "EPI358781_A_Fujian_411_2002_NA_NA_NA_NA")
G3 <- G3 %>% 
  filter(compareStrain != "EPI367109_A_California_7_2004_NA_NA_NA_NA")
G4 <- G4 %>% 
  filter(compareStrain != "EPI502253_A_Wisconsin_67_2005_NA_NA_NA_NA")
G5 <- G5 %>% 
  filter(compareStrain != "EPI577980_A_Brisbane_10_2007_NA_NA_NA_NA")
G6 <- G6 %>% 
  filter(compareStrain != "EPI577969_A_Perth_16_2009_NA_NA_NA_NA")
G7 <- G7 %>% 
  filter(compareStrain != "EPI417234_A_Victoria_361_2011_NA_NA_NA_NA")
G8 <- G8 %>% 
  filter(compareStrain != "EPI614441_A_Switzerland_9715293_2013_NA_NA_NA_NA")
G9 <- G9 %>% 
  filter(compareStrain != "EPI686117_A_Hong_Kong_15611_2015_NA_NA_NA")
G10 <- G10 %>% 
  filter(compareStrain != "MG974447_A_Kansas_14_2017_NA_NA_NA_NA")
G11 <- G11 %>% 
  filter(compareStrain != "1592032_A_Hong_Kong_2671_2019_NA_NA_NA")
G12 <- G12 %>% 
  filter(compareStrain != "1841681_A_Cambodia_e0826360_2020_NA_NA_NA_NA")
G13 <- G13 %>% 
  filter(compareStrain != "EPI2415906_A_Darwin_9_2021_NA_NA_NA_NA")

AllG <- rbind(G1, G2, G3, G4, G5, G6, G7, G8, G9, G10, G11, G12, G13)

## Clear the name of compare strain 
AllG$compareStrain_V1 <- gsub("_NA_NA_NA_NA", "", AllG$compareStrain)
AllG$compareStrain_V1 <- gsub("_NA_NA_NA", "", AllG$compareStrain_V1)
AllG$compareStrain_V1 <- gsub("_NA_NA", "", AllG$compareStrain_V1)
AllG$compareStrain_V1 <- gsub("2023_NA", "2023", AllG$compareStrain_V1)
AllG$compareStrain <- gsub("2009_NA", "2009", AllG$compareStrain_V1)

## linking two data
full_data <- left_join(Index_725, AllG, by = "compareStrain")

write.csv(full_data[-1,-1], "Data/FullReg_724.csv")
