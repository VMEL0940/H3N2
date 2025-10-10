## Package ##
library(tidyverse)
library(lme4)
library(pROC)
library(lmtest)
library(corrplot)

## Set working directory ##
setwd("/01_GenDisFlu/")

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

VacCoCol <- c("#716172","#EFA95F","#3EA3C5","#6D9431","#800000","#9E6D13","#092E50","#769F74")

train_regdata$subtMRCA_1stDesc[is.na(train_regdata$subtMRCA_1stDesc)] <- 2
train_regdata$Dominance <- factor(train_regdata$subtMRCA_1stDesc, 
                         levels = c(1,0),
                         labels = c("Dominant","Extinct"))

Dom <- c("#96fd29","#000000")


subtMRCA_Neb_Prop <- train_regdata %>% 
  group_by(Vaccine_code, subtMRCA_1stDesc) %>% 
  summarise(n = n())


## Koel ##
train_regdata %>% 
  ggplot() +
  geom_jitter(aes(x= Vaccine_code, y= HA_Koel, fill = Auspice), size = 1.5, shape = 21, color = "black", width = 0.25) +
  geom_boxplot(aes(x= Vaccine_code, y= HA_Koel, color = Auspice), size = 1, alpha = 0) +
  scale_fill_manual(values = Dom) +
  scale_color_manual(values = Dom) +
  ylab("The number of Amino acid change in the Koel Position")+
  xlab("Vaccine Period")+
  theme_bw()

## RBS ##
train_regdata %>% 
  ggplot() +
  geom_jitter(aes(x= Vaccine_code, y= HA_RBD, fill = Auspice), size = 1.5, shape = 21, color = "black", width = 0.25) +
  geom_boxplot(aes(x= Vaccine_code, y= HA_RBD, color = Auspice), size = 1, alpha = 0) +
  scale_fill_manual(values = Dom) +
  scale_color_manual(values = Dom) +
  ylab("The number of Amino acid change in the RBS")+
  xlab("Vaccine Period")+
  theme_bw()


## 15A ##
train_regdata %>% 
  ggplot() +
  geom_jitter(aes(x= Vaccine_code, y= HA_15A, fill = Auspice), size = 1.5, shape = 21, color = "black", width = 0.25) +
  geom_boxplot(aes(x= Vaccine_code, y= HA_15A, color = Auspice), size = 1, alpha = 0) +
  scale_color_manual(values = Dom) +
  scale_fill_manual(values = Dom) +
  ylab("The number of Amino acid change around 15A near the RBS")+
  xlab("Vaccine Period")+
  theme_bw()



jpeg(filename = "/01_GenDisFlu/Fig/Desc/Dom/Dom_HA_N_Box.jpg", width = 150, height = 150, units = "mm", res = 520)

## HA ##
train_regdata %>% 
  ggplot() +
  geom_jitter(aes(x= Vaccine_code, y= HA, fill = Dominance), size = 1.5, shape = 21, color = "black", width = 0.25) +
  geom_boxplot(aes(x= Vaccine_code, y= HA, color = Dominance), size = 1, alpha = 0) +
  scale_color_manual(values = Dom) +
  scale_fill_manual(values = Dom) +
  ylab("Nonsynonymous genetic distance of the HA gene")+
  xlab("Period of strain spread")+
  theme_bw()

dev.off()


jpeg(filename = "/01_GenDisFlu/Fig/Desc/Dom/Dom_NA_N_Box.jpg", width = 150, height = 150, units = "mm", res = 520)

## NA ##
train_regdata %>% 
  ggplot() +
  geom_jitter(aes(x= Vaccine_code, y= N_A, fill = Dominance), size = 1.5, shape = 21, color = "black", width = 0.25) +
  geom_boxplot(aes(x= Vaccine_code, y= N_A, color = Dominance), size = 1, alpha = 0) +
  scale_color_manual(values = Dom) +
  scale_fill_manual(values = Dom) +
  ylab("Nonsynonymous genetic distance of the NA gene")+
  xlab("Period of strain spread")+
  theme_bw()

dev.off()

jpeg(filename = "/01_GenDisFlu/Fig/Desc/Dom/Dom_PA_N_Box.jpg", width = 150, height = 150, units = "mm", res = 520)

## PA ##
train_regdata %>% 
  ggplot() +
  geom_jitter(aes(x= Vaccine_code, y= PA, fill = Dominance), size = 1.5, shape = 21, color = "black", width = 0.25) +
  geom_boxplot(aes(x= Vaccine_code, y= PA, color = Dominance), size = 1, alpha = 0) +
  scale_color_manual(values = Dom) +
  scale_fill_manual(values = Dom) +
  ylab("Nonsynonymous genetic distance  of the PA gene")+
  xlab("Period of strain spread")+
  theme_bw()

dev.off()

jpeg(filename = "/01_GenDisFlu/Fig/Desc/Dom/Dom_M_N_Box.jpg", width = 150, height = 150, units = "mm", res = 520)


## M ##
train_regdata %>% 
  ggplot() +
  geom_jitter(aes(x= Vaccine_code, y= M, fill = Dominance), size = 1.5, shape = 21, color = "black", width = 0.25) +
  geom_boxplot(aes(x= Vaccine_code, y= M, color = Dominance), size = 1, alpha = 0) +
  scale_color_manual(values = Dom) +
  scale_fill_manual(values = Dom) +
  ylab("Nonsynonymous genetic distance of the M gene")+
  xlab("Period of strain spread")+
  theme_bw()


dev.off()

jpeg(filename = "/01_GenDisFlu/Fig/Desc/Dom/Dom_NP_N_Box.jpg", width = 150, height = 150, units = "mm", res = 520)
## NP ##
train_regdata %>% 
  ggplot() +
  geom_jitter(aes(x= Vaccine_code, y= NP, fill = Dominance), size = 1.5, shape = 21, color = "black", width = 0.25) +
  geom_boxplot(aes(x= Vaccine_code, y= NP, color = Dominance), size = 1, alpha = 0) +
  scale_color_manual(values = Dom) +
  scale_fill_manual(values = Dom) +
  ylab("Nonsynonymous genetic distance of the NP gene")+
  xlab("Period of strain spread")+
  theme_bw()
dev.off()

jpeg(filename = "/01_GenDisFlu/Fig/Desc/Dom/Dom_PB2_N_Box.jpg", width = 150, height = 150, units = "mm", res = 520)
## PB2 ##
train_regdata %>% 
  ggplot() +
  geom_jitter(aes(x= Vaccine_code, y= PB2, fill = Dominance), size = 1.5, shape = 21, color = "black", width = 0.25) +
  geom_boxplot(aes(x= Vaccine_code, y= PB2, color = Dominance), size = 1, alpha = 0) +
  scale_color_manual(values = Dom) +
  scale_fill_manual(values = Dom) +
  ylab("Nonsynonymous genetic distance of the PB2 gene")+
  xlab("Period of strain spread")+
  theme_bw()
dev.off()

jpeg(filename = "/01_GenDisFlu/Fig/Desc/Dom/Dom_PB1_N_Box.jpg", width = 150, height = 150, units = "mm", res = 520)
## PB1 ##
train_regdata %>% 
  ggplot() +
  geom_jitter(aes(x= Vaccine_code, y= PB1, fill = Dominance), size = 1.5, shape = 21, color = "black", width = 0.25) +
  geom_boxplot(aes(x= Vaccine_code, y= PB1, color = Dominance), size = 1, alpha = 0) +
  scale_color_manual(values = Dom) +
  scale_fill_manual(values = Dom) +
  ylab("Nonsynonymous genetic distance of the PB1 gene")+
  xlab("Period of strain spread")+
  theme_bw()
dev.off()

jpeg(filename = "/01_GenDisFlu/Fig/Desc/Dom/Dom_NS_N_Box.jpg", width = 150, height = 150, units = "mm", res = 520)
## NS ##
train_regdata %>% 
  ggplot() +
  geom_jitter(aes(x= Vaccine_code, y= NS, fill = Dominance), size = 1.5, shape = 21, color = "black", width = 0.25) +
  geom_boxplot(aes(x= Vaccine_code, y= NS, color = Dominance), size = 1, alpha = 0) +
  scale_color_manual(values = Dom) +
  scale_fill_manual(values = Dom) +
  ylab("Nonsynonymous genetic distance of the NS gene")+
  xlab("Period of strain spread")+
   theme_bw()
dev.off()



jpeg(filename = "/01_GenDisFlu/Fig/Desc/Dom/Cor_N_Box.jpg", width = 150, height = 150, units = "mm", res = 520)

# CorMatrix for the selection
cormat <- train_regdata[,c(3:10,24)]

colnames(cormat)[9] <- "Dominance"
colnames(cormat)[3] <- "NA"

cormat_2 <- cormat[c(8,7,6,1,4,3,2,5,9)]

M <- cor(cormat_2)
testRes = cor.mtest(cormat_2, conf.level = 0.95)

## Correlation plot
corrplot(M, p.mat = testRes$p, method = 'circle', insig='blank', tl.srt = 0, 
         sig.level = 0.05, addCoef.col ='black', number.cex = 0.8)

dev.off()





