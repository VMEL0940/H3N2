## Package ##
library(tidyverse)
library(lme4)
library(pROC)
library(lmtest)
library(corrplot)

## Set working directory ##
setwd("/01_GenDisFlu/")
#setwd("01_GenDisFlu/") ## for mac

## Data import ##
AllG_NS <- read.csv("Data/Genetic Distances/GeneCompete/AllG_NS.csv", header = T, na.strings = "")

AllG_NS$compareStrain_V1 <- gsub("_NA_NA_NA_NA", "", AllG_NS$compareStrain)
AllG_NS$compareStrain_V1 <- gsub("_NA_NA_NA", "", AllG_NS$compareStrain_V1)
AllG_NS$compareStrain_V1 <- gsub("_NA_NA", "", AllG_NS$compareStrain_V1)
AllG_NS$compareStrain <- gsub("2009_NA", "2009", AllG_NS$compareStrain_V1)

Ind <- AllG_NS[,c(3,2,4,6)] %>% 
  spread(key = Gene, value = Nmedian)

colnames(Ind)[1] <- "ID"

## leftjoin * subset ##
Dom <- read.csv("Data/Dom.csv", header = T, na.strings = "")

regdata <- right_join(Ind, Dom, by = "ID")

colnames(regdata)[5] <- "N_A"

### Subset data - Exclude HK15 ##
train_regdata <- regdata %>% 
  filter(vaccineStrain != "EPI686117_A_Hong_Kong_15611_2015_NA_NA_NA")

train_regdata$Vaccine_code <- factor(train_regdata$vaccineStrain, 
                               levels = c("EPI103320_A_Moscow_10_1999_NA_NA_NA_NA",
                                          "EPI358781_A_Fujian_411_2002_NA_NA_NA_NA", 
                                          "EPI367109_A_California_7_2004_NA_NA_NA_NA",
                                          "EPI502253_A_Wisconsin_67_2005_NA_NA_NA_NA",
                                          "EPI577980_A_Brisbane_10_2007_NA_NA_NA_NA",
                                          "EPI577969_A_Perth_16_2009_NA_NA_NA_NA", 
                                          "EPI417234_A_Victoria_361_2011_NA_NA_NA_NA",
                                          "EPI614441_A_Switzerland_9715293_2013_NA_NA_NA_NA"),
                               labels = c("1.Mos99",
                                          "2.Fuj02", 
                                          "3.Cal04",
                                          "4.Wis05",
                                          "5.Bris07",
                                          "6.Prth09", 
                                          "7.Vic11",
                                          "8.Swtz13"))

train_regdata$Dominance <- factor(train_regdata$subtMRCA_1stDesc, 
                         levels = c(1,0),
                         labels = c("Dominant","Extinct"))

subtMRCA_Neb_Prop <- train_regdata %>% 
  group_by(Vaccine_code, subtMRCA_1stDesc) %>% 
  summarise(n = n())

subtMRCA_Prop <- train_regdata %>% 
  group_by(subtMRCA_1stDesc) %>% 
  summarise(n = n())

#### 1. Chisq.test - Reassortment #####

chisq.test(train_regdata$subtMRCA_1stDesc, train_regdata$Reasrt)

train_regdata %>% 
  group_by(Reasrt, subtMRCA_1stDesc) %>% 
  summarise(n = n())

#chisq.test(train_regdata$subtMRCA_Neibhor_Dom, train_regdata$CoEvo_PB2)
#chisq.test(train_regdata$subtMRCA_Neibhor_Dom, train_regdata$CoEvo_PB1) ## NO 1 Value 
#chisq.test(train_regdata$subtMRCA_Neibhor_Dom, train_regdata$CoEvo_PA) ## One 1 value
#chisq.test(train_regdata$subtMRCA_Neibhor_Dom, train_regdata$CoEvo_NP) ## One 1 value
#chisq.test(train_regdata$subtMRCA_Neibhor_Dom, train_regdata$CoEvo_NA) ## One 1 value
#chisq.test(train_regdata$subtMRCA_Neibhor_Dom, train_regdata$CoEvo_M) ## One 1 value
#chisq.test(train_regdata$subtMRCA_Neibhor_Dom, train_regdata$CoEvo_NS) ## Twp 1 value


