## Package ##
library(tidyverse)

## Set working directory ##
setwd("01_GenDisFlu/")
#setwd("01_GenDisFlu/") ## for mac

## Evaluate the 601 Data -- Data 1 ###

## Data import ##
Meta <- read.csv("Data/AllData/H3N2_Traindata_600strains_1_Metadata.csv", header = T, na.strings = "")
NonSyn <- read.csv("Data/AllData/H3N2_Traindata_600strains_2_NonsynDist.csv", header = T, na.strings = "")
Syn <- read.csv("Data/AllData/H3N2_Traindata_600strains_3_SynDist.csv", header = T, na.strings = "")
Others <- read.csv("Data/AllData/H3N2_Traindata_600strains_4_OtherPredictors.csv", header = T, na.strings = "")

## Combine Data ##
All600 <- left_join(Meta, NonSyn[,c(-2,-3)], by = "ID")
All600 <- left_join(All600, Syn[,c(-2,-3)], by = "ID")
All600 <- left_join(All600, Others[,c(-2,-3)], by = "ID")

## Drop HK15 ##
All600_Train <- All600 %>% 
  filter(vaccine_code != "HK15")

## Association between HA Reassortment and   ##
### (1) PB2 -- No Association ###
chisq.test(table(All600_Train$Dom, All600_Train$HA_PB2_Reassort))
#fisher.test(table(All600_Train$Dom, All600_Train$HA_PB2_Reassort))
### (2) PB1 -- No Association ###
chisq.test(table(All600$Dom, All600$HA_PB1_Reassort))
#fisher.test(table(All600$Dom, All600$HA_PB1_Reassort))
### (3) PA -- No Association ###
chisq.test(table(All600$Dom, All600$HA_PA_Reassort))
#fisher.test(table(All600$Dom, All600$HA_PA_Reassort))
### (4) NP -- No Association ###
chisq.test(table(All600$Dom, All600$HA_NP_Reassort))
#fisher.test(table(All600$Dom, All600$HA_NP_Reassort))
### (5) NA -- No Association ###
chisq.test(table(All600$Dom, All600$HA_NA_Reassort))
#fisher.test(table(All600$Dom, All600$HA_NA_Reassort))
### (6) M -- No Association ###
chisq.test(table(All600$Dom, All600$HA_M_Reassort))
#fisher.test(table(All600$Dom, All600$HA_M_Reassort))
### (7) NS -- No Association ###
chisq.test(table(All600$Dom, All600$HA_NS_Reassort))
#fisher.test(table(All600$Dom, All600$HA_NS_Reassort))
### (8) All -- No Association ###
chisq.test(table(All600$Dom, All600$All_Reassort))
#fisher.test(table(All600$Dom, All600$All_Reassort))

