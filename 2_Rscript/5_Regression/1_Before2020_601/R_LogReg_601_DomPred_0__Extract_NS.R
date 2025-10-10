## Package ##
library(tidyverse)
library(lme4)
library(pROC)
library(lmtest)
library(corrplot)

## Set working directory ##
setwd("01_GenDisFlu/")
#setwd("01_GenDisFlu/") ## for mac

## Data import ##
AllG_NS <- read.csv("Data/Genetic Distances/GeneCompete/AllG_NS.csv", header = T, na.strings = "")

AllG_NS$compareStrain_V1 <- gsub("_NA_NA_NA_NA", "", AllG_NS$compareStrain)
AllG_NS$compareStrain_V1 <- gsub("_NA_NA_NA", "", AllG_NS$compareStrain_V1)
AllG_NS$compareStrain_V1 <- gsub("_NA_NA", "", AllG_NS$compareStrain_V1)
AllG_NS$compareStrain <- gsub("2009_NA", "2009", AllG_NS$compareStrain_V1)

Ind <- AllG_NS[,c(3,2,4,6)] %>% 
  spread(key = Gene, value = Nmedian)

colnames(Ind) <- c("ID","VacStr","HA_Nonsyn","M_Nonsyn","NA_Nonsyn","NP_Nonsyn","NS_Nonsyn","PA_Nonsyn","PB1_Nonsyn","PB2_Nonsyn")

Ind2 <- AllG_NS[,c(3,4,9)] %>% 
  spread(key = Gene, value = Smedian)

colnames(Ind2) <- c("ID","HA_Syn","M_Syn","NA_Syn","NP_Syn","NS_Syn","PA_Syn","PB1_Syn","PB2_Syn")

All_Ind <- right_join(Ind, Ind2, by = "ID")

write.csv(All_Ind, "Data/H3N2_600_GenD.csv")
