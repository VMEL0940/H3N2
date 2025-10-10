### Seasonal H3N2 Prediction - N and S Genetic distance descriptive ###
## Package ##
library(tidyverse)

## Set working directory ##
setwd("01_GenDisFlu/Data/")


AllG_NS <- read.csv("Genetic Distances/Verification_702/Summary/AllG_NS.csv", header = T, na.strings = "")

### All Group ###

AllG_NS$Gene <- factor(AllG_NS$Gene, levels = c("PB2", "PB1", "PA", "HA", 
                              "NP", "NA", "M", "NS"))

GeneCol <- c("#C780E8", "#9D94FF","#59ADF6","#08CAD1","#42D6A4","#F8F38D","#FFB480","#FF6961")


AllG_NS$Vaccine_code <- factor(AllG_NS$vaccineStrain, 
                                levels = c("EPI103320_A_Moscow_10_1999_NA_NA_NA_NA",
                                           "EPI358781_A_Fujian_411_2002_NA_NA_NA_NA", 
                                           "EPI367109_A_California_7_2004_NA_NA_NA_NA",
                                           "EPI502253_A_Wisconsin_67_2005_NA_NA_NA_NA",
                                           "EPI577980_A_Brisbane_10_2007_NA_NA_NA_NA",
                                           "EPI577969_A_Perth_16_2009_NA_NA_NA_NA", 
                                           "EPI417234_A_Victoria_361_2011_NA_NA_NA_NA",
                                           "EPI614441_A_Switzerland_9715293_2013_NA_NA_NA_NA",
                                           "EPI686117_A_Hong_Kong_15611_2015_NA_NA_NA",
                                           "MG974447_A_Kansas_14_2017_NA_NA_NA_NA",
                                           "1592032_A_Hong_Kong_2671_2019_NA_NA_NA", 
                                           "1846767_A_Cambodia_E0826360_2020_NA_NA_NA_NA",
                                           "2415906_A_Darwin_9_2021_NA_NA_NA_NA"),
                                labels = c("1.Mos99",
                                           "2.Fuj02", 
                                           "3.Cal04",
                                           "4.Wis05",
                                           "5.Bris07",
                                           "6.Prth09", 
                                           "7.Vic11",
                                           "8.Swtz13",
                                           "9.HK15",
                                           "10.Kan17",
                                           "11.HK19",
                                           "12.Cam20",
                                           "13.Dar21"))

## Descriptive Summary of N & S value by Gene and Vaccine period ##
Desc_all <- AllG_NS %>% 
  group_by(Gene, Vaccine_code) %>% 
  summarise(Nmean = mean(Nmedian), Nse = 1.96*sd(Nmedian)/sqrt(n()),
            Smean = mean(Smedian), Sse = 1.96*sd(Smedian)/sqrt(n()),
            n = n())

Desc_all_N_All <- AllG_NS %>% 
  group_by(Gene) %>% 
  summarise(DisMean = mean(Nmedian), DisMedian = median(Nmedian),
            lowerCrIn = quantile(Nmedian, na.rm = T, probs = 0.05),
            UpperCrIn = quantile(Nmedian, na.rm = T, probs = 0.95),
            lower = mean(Nmedian) - 1.96*sd(Nmedian)/sqrt(n()),
            upper = mean(Nmedian) + 1.96*sd(Nmedian)/sqrt(n()),
            n = n())


Desc_all_S_All <- AllG_NS %>% 
  group_by(Gene) %>% 
  summarise(DisMean = mean(Smedian), DisMedian = median(Smedian),
            lowerCrIn = quantile(Smedian, na.rm = T, probs = 0.05),
            UpperCrIn = quantile(Smedian, na.rm = T, probs = 0.95),
            lower = mean(Smedian) - 1.96*sd(Smedian)/sqrt(n()),
            upper = mean(Smedian) + 1.96*sd(Smedian)/sqrt(n()),
            n = n())



Desc_all_N <- AllG_NS %>% 
  group_by(Gene, Vaccine_code) %>% 
  summarise(DisMean = mean(Nmedian), DisMedian = median(Nmedian),
            lowerCrIn = quantile(Nmedian, na.rm = T, probs = 0.05),
            UpperCrIn = quantile(Nmedian, na.rm = T, probs = 0.95),
            lower = mean(Nmedian) - 1.96*sd(Nmedian)/sqrt(n()),
            upper = mean(Nmedian) + 1.96*sd(Nmedian)/sqrt(n()),
            n = n())


Desc_all_S <- AllG_NS %>% 
  group_by(Gene, Vaccine_code) %>% 
  summarise(DisMean = mean(Smedian), DisMedian = median(Smedian),
            lowerCrIn = quantile(Smedian, na.rm = T, probs = 0.05),
            UpperCrIn = quantile(Smedian, na.rm = T, probs = 0.95),
            lower = mean(Smedian) - 1.96*sd(Smedian)/sqrt(n()),
            upper = mean(Smedian) + 1.96*sd(Smedian)/sqrt(n()),
            n = n())

#write.csv(Desc_all_N, file="Data/Genetic Distances/GeneCompete/Desc_all_N.csv")
#write.csv(Desc_all_S, file="Data/Genetic Distances/GeneCompete/Desc_all_S.csv")
#write.csv(Desc_all_N_All, file="Data/Genetic Distances/GeneCompete/Desc_all_N_All.csv")
#write.csv(Desc_all_S_All, file="Data/Genetic Distances/GeneCompete/Desc_all_S_All.csv")
#write.csv(AllG_NS,file="Data/Genetic Distances/GeneCompete/AllG_NS.csv")


### Only After 2020 ##
jpeg(filename = "01_GenDisFlu/Fig/Desc/Onlyafter2020/All_N_Boxjitter.jpeg", width = 150, height = 130, units = "mm", res = 320)

## Box plot for N value by Gene ##
AllG_NS %>% 
  filter(Group == "G10" | Group == "G11" | Group == "G12" | Group == "G13") %>% 
  ggplot() +
  geom_jitter(aes(x= Gene, y= Nmedian), width = 0.25, size = 0.50, color = "gray40") +
  geom_boxplot(aes(x= Gene, y = Nmedian), color = GeneCol, size = 1.3, outlier.shape = NA, alpha = 0) +
  scale_color_manual(values = GeneCol, guide="none") +
  scale_fill_manual(values = GeneCol, guide="none") +
  ylab("Nonsynonymous genetic distance")+
  xlab("Gene segment")+
  ylim(-1, 60) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12)) 

dev.off()



jpeg(filename = "01_GenDisFlu/Fig/Desc/Onlyafter2020/All_S_Boxjitter.jpeg", width = 150, height = 130, units = "mm", res = 320)

## Box plot for S value by Gene ##
AllG_NS %>% 
  filter(Group == "G10" | Group == "G11" | Group == "G12" | Group == "G13") %>% 
  ggplot() +
  geom_jitter(aes(x= Gene, y= Smedian), width = 0.25, size = 0.6, color = "gray40") +
  geom_boxplot(aes(x= Gene, y = Smedian), color = GeneCol, size = 1.3, outlier.shape = NA, alpha = 0) +
  scale_color_manual(values = GeneCol, guide="none") +
  scale_fill_manual(values = GeneCol, guide="none") +
  ylab("Synonymous genetic distance")+
  xlab("Gene segment")+
  ylim(-1, 65) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12)) 

dev.off()


