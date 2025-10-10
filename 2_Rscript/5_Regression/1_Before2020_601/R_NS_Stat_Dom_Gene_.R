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

All600_Train$Dominance <- factor(All600_Train$Dom, 
                            levels = c(0,1),
                            labels = c("Extinct", "Dominant"))

wilcox.test(PB2_Nonsyn ~ Dom, data = All600_Train)
wilcox.test(PB1_Nonsyn ~ Dom, data = All600_Train)
wilcox.test(PA_Nonsyn ~ Dom, data = All600_Train)
wilcox.test(HA_Nonsyn ~ Dom, data = All600_Train)
wilcox.test(NP_Nonsyn ~ Dom, data = All600_Train)
wilcox.test(NA_Nonsyn ~ Dom, data = All600_Train)
wilcox.test(M_Nonsyn ~ Dom, data = All600_Train)
wilcox.test(NS_Nonsyn ~ Dom, data = All600_Train)

wilcox.test(PB2_Syn ~ Dom, data = All600_Train)
wilcox.test(PB1_Syn ~ Dom, data = All600_Train)
wilcox.test(PA_Syn ~ Dom, data = All600_Train)
wilcox.test(HA_Syn ~ Dom, data = All600_Train)
wilcox.test(NP_Syn ~ Dom, data = All600_Train)
wilcox.test(NA_Syn ~ Dom, data = All600_Train)
wilcox.test(M_Syn ~ Dom, data = All600_Train)
wilcox.test(NS_Syn ~ Dom, data = All600_Train)

wilcox.test(HA_RBD ~ Dom, data = All600_Train)
wilcox.test(HA_15A ~ Dom, data = All600_Train)
wilcox.test(HA_Koel ~ Dom, data = All600_Train)

