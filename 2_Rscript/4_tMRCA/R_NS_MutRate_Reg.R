### Influ H3N2 GenDist Proj ###
## Package ##
library(tidyverse)
library(lubridate)
library(lme4)
library(lmerTest)
## Set working directory ##
setwd("/01_GenDisFlu/")

## Data import ##
AllG_NS <- read.csv("Data/Genetic Distances/GeneCompete/AllG_NS.csv", header = T, na.strings = "")
Subsetlist <- read.csv("Data/H3N2_Index_20230331.csv")

## Set up the Factors
AllG_NS$Gene <- factor(AllG_NS$Gene, levels = c("PB2", "PB1", "PA", "HA", 
                                                "NP", "NA", "M", "NS"))
AllG_NS$Vaccine_code <- factor(AllG_NS$Vaccine_code,levels = c("1.Mos99","2.Fuj02","3.Cal04",
                                                               "4.Wis05","5.Bris07","6.Prth09",
                                                               "7.Vic11", "8.Swtz13", "9.HK15"))

VacCoCol <- c("#716172", "#EFA95F","#3EA3C5","#6D9431","#800000","#9E6D13","#092E50","#769F74","#954495")


## Matching compare strains
AllG_NS$compareStrain_V1 <- gsub("_NA_NA_NA_NA", "", AllG_NS$compareStrain)
AllG_NS$compareStrain_V1 <- gsub("_NA_NA_NA", "", AllG_NS$compareStrain_V1)
AllG_NS$compareStrain_V1 <- gsub("_NA_NA", "", AllG_NS$compareStrain_V1)
AllG_NS$compareStrain <- gsub("2009_NA", "2009", AllG_NS$compareStrain_V1)
AllG_NS <- AllG_NS[,-12]

## adding date + Modify date
AllG_NS_Date <- full_join(AllG_NS, Subsetlist[,c(1,4:6)], by = "compareStrain")
AllG_NS_Date$ColDate <- paste0(AllG_NS_Date$Year, "/", AllG_NS_Date$Month, "/", AllG_NS_Date$Day)
AllG_NS_Date$ColDate <- as.Date(AllG_NS_Date$ColDate,"%Y/%m/%d")

## convert decimate date and align the date 
AllG_NS_Date$Ydate <- decimal_date(AllG_NS_Date$ColDate)
AllG_NS_Date$Regdate <- 0 

AllG_NS_Date[AllG_NS_Date$Vaccine_code == "1.Mos99",]$Regdate <- AllG_NS_Date[AllG_NS_Date$Vaccine_code == "1.Mos99",]$Ydate - 1999.748
AllG_NS_Date[AllG_NS_Date$Vaccine_code == "2.Fuj02",]$Regdate <- AllG_NS_Date[AllG_NS_Date$Vaccine_code == "2.Fuj02",]$Ydate - 2004.748
AllG_NS_Date[AllG_NS_Date$Vaccine_code == "3.Cal04",]$Regdate <- AllG_NS_Date[AllG_NS_Date$Vaccine_code == "3.Cal04",]$Ydate - 2005.748
AllG_NS_Date[AllG_NS_Date$Vaccine_code == "4.Wis05",]$Regdate <- AllG_NS_Date[AllG_NS_Date$Vaccine_code == "4.Wis05",]$Ydate - 2006.748
AllG_NS_Date[AllG_NS_Date$Vaccine_code == "5.Bris07",]$Regdate <- AllG_NS_Date[AllG_NS_Date$Vaccine_code == "5.Bris07",]$Ydate - 2008.748
AllG_NS_Date[AllG_NS_Date$Vaccine_code == "6.Prth09",]$Regdate <- AllG_NS_Date[AllG_NS_Date$Vaccine_code == "6.Prth09",]$Ydate - 2010.748
AllG_NS_Date[AllG_NS_Date$Vaccine_code == "7.Vic11",]$Regdate <- AllG_NS_Date[AllG_NS_Date$Vaccine_code == "7.Vic11",]$Ydate - 2012.748
AllG_NS_Date[AllG_NS_Date$Vaccine_code == "8.Swtz13",]$Regdate <- AllG_NS_Date[AllG_NS_Date$Vaccine_code == "8.Swtz13",]$Ydate - 2015.748
AllG_NS_Date[AllG_NS_Date$Vaccine_code == "9.HK15",]$Regdate <- AllG_NS_Date[AllG_NS_Date$Vaccine_code == "9.HK15",]$Ydate - 2016.748

## Separate value of two lineages of Swtz13 and HK15
AllG_NS_Date$Vaccine_code_Reg <- as.character(AllG_NS_Date$Vaccine_code)

vec_swtz13 <- AllG_NS_Date[AllG_NS_Date$Vaccine_code == "8.Swtz13" & AllG_NS_Date$Gene == "HA" & AllG_NS_Date$Nmedian >= 12,]$compareStrain
vec_HK15 <- AllG_NS_Date[AllG_NS_Date$Vaccine_code == "9.HK15" & AllG_NS_Date$Gene == "HA" & AllG_NS_Date$Nmedian >= 12,]$compareStrain

AllG_NS_Date$Vaccine_code_Reg <- ifelse(AllG_NS_Date$compareStrain %in% vec_swtz13, "8.Swtz13_2", AllG_NS_Date$Vaccine_code_Reg)
AllG_NS_Date$Vaccine_code_Reg <- ifelse(AllG_NS_Date$compareStrain %in% vec_HK15, "9.HK15_2", AllG_NS_Date$Vaccine_code_Reg)

unique(AllG_NS_Date$Vaccine_code_Reg)
  
jpeg(filename = "/01_GenDisFlu/Fig/Desc/NS_Time/HA_N_Time_VaccG.jpeg", width = 180, height = 180, units = "mm", res = 400)

## Multi-intercept concept of 1. HA - N
AllG_NS_Date %>%
  filter(Gene == "HA") %>% 
  ggplot(aes(x = Regdate, y = Nmedian, color = Vaccine_code)) +
  geom_point(aes(fill = Vaccine_code), size = 2, ) +
  scale_x_continuous(name = "Year") +
  ylab("Number of Non-synonymous mutation of the HA gene") +
  scale_fill_manual(values = VacCoCol, name = "Vaccine Strains", guide="none") +
  scale_color_manual(values = VacCoCol, guide="none") +
  scale_y_continuous(breaks = seq(0, 35, by = 5)) +
  geom_hline(yintercept = 0) +
  geom_smooth(se = FALSE, method = lm) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9))

dev.off()

HA_N_model <- lmer(Nmedian ~ Regdate + (1 | Vaccine_code), data=AllG_NS_Date[AllG_NS_Date$Gene == "HA",])
summary(HA_N_model)


## Multi-intercept concept of 1. HA - S
AllG_NS_Date %>%
  filter(Gene == "HA") %>% 
  ggplot(aes(x = Regdate, y = Smedian, color = Vaccine_code)) +
  geom_point( size = 1.5) +
  scale_x_continuous(name = "Year") +
  ylab("Number of Synonymous mutation") +
  scale_color_discrete(name = "Vaccine Strains") +
  scale_y_continuous(breaks = seq(0, 35, by = 5)) +
  geom_hline(yintercept = 0) +
  geom_smooth(se = FALSE, method = lm) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9))

HA_S_model <- lmer(Smedian ~ Regdate + (1 | Vaccine_code), data=AllG_NS_Date[AllG_NS_Date$Gene == "HA",])
summary(HA_S_model)



unique(AllG_NS_Date[AllG_NS_Date$Gene == "PB2",]$Vaccine_code_Reg)


jpeg(filename = "/01_GenDisFlu/Fig/Desc/NS_Time/PB2_N_Time_VaccG.jpeg", width = 180, height = 180, units = "mm", res = 400)

## Multi-intercept concept of 1. PB2 - N
AllG_NS_Date %>%
  filter(Gene == "PB2") %>% 
  ggplot(aes(x = Regdate, y = Nmedian, color = Vaccine_code)) +
  geom_point(aes(fill = Vaccine_code), size = 2) +
  scale_x_continuous(name = "Year") +
  ylab("Number of Non-synonymous mutation of the PB2 gene") +
  scale_fill_manual(values = VacCoCol, name = "Vaccine Strains", guide="none") +
   scale_color_manual(values = VacCoCol, guide="none") +
  scale_y_continuous(breaks = seq(0, 35, by = 5)) +
  geom_hline(yintercept = 0) +
  geom_smooth(se = FALSE, method = lm) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9))

dev.off()

PB2_N_model <- lmer(Nmedian ~ Regdate + (1 | Vaccine_code_Reg), data=AllG_NS_Date[AllG_NS_Date$Gene == "PB2",])
summary(PB2_N_model)

PB2_N_model <- lmer(Nmedian ~ Regdate + (1 | Vaccine_code), data=AllG_NS_Date[AllG_NS_Date$Gene == "PB2",])
summary(PB2_N_model)

jpeg(filename = "/01_GenDisFlu/Fig/Desc/NS_Time/PB1_N_Time_VaccG.jpeg", width = 180, height = 180, units = "mm", res = 520)


## Multi-intercept concept of 2. PB1 - N
AllG_NS_Date %>%
  filter(Gene == "PB1") %>% 
  ggplot(aes(x = Regdate, y = Nmedian, color = Vaccine_code)) +
  geom_point(aes(fill = Vaccine_code), size = 2) +
  scale_x_continuous(name = "Year") +
  ylab("Number of Non-synonymous mutation of the PB1 gene") +
  scale_fill_manual(values = VacCoCol, name = "Vaccine Strains", guide="none") +
   scale_color_manual(values = VacCoCol, guide="none") +
  scale_y_continuous(breaks = seq(0, 35, by = 5)) +
  geom_hline(yintercept = 0) +
  geom_smooth(se = FALSE, method = lm) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9))

dev.off()

PB1_N_model <- lmer(Nmedian ~ Regdate + (1 | Vaccine_code_Reg), data=AllG_NS_Date[AllG_NS_Date$Gene == "PB1",])
summary(PB1_N_model)


PB1_N_model <- lmer(Nmedian ~ Regdate + (1 | Vaccine_code), data=AllG_NS_Date[AllG_NS_Date$Gene == "PB1",])
summary(PB1_N_model)

jpeg(filename = "/01_GenDisFlu/Fig/Desc/NS_Time/PA_N_Time_VaccG.jpeg", width = 180, height = 180, units = "mm", res = 520)

## Multi-intercept concept of 3. PA - N
AllG_NS_Date %>%
  filter(Gene == "PA") %>% 
  ggplot(aes(x = Regdate, y = Nmedian, color = Vaccine_code)) +
  geom_point(aes(fill = Vaccine_code), size = 2) +
  scale_x_continuous(name = "Year") +
  ylab("Number of Non-synonymous mutation of the PA gene") +
  scale_fill_manual(values = VacCoCol, name = "Vaccine Strains", guide="none") +
   scale_color_manual(values = VacCoCol, guide="none") +
  scale_y_continuous(breaks = seq(0, 35, by = 5)) +
  geom_hline(yintercept = 0) +
  geom_smooth(se = FALSE, method = lm) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9))

dev.off()


PA_N_model <- lmer(Nmedian ~ Regdate + (1 | Vaccine_code_Reg), data=AllG_NS_Date[AllG_NS_Date$Gene == "PA",])
summary(PA_N_model)

PA_N_model <- lmer(Nmedian ~ Regdate + (1 | Vaccine_code), data=AllG_NS_Date[AllG_NS_Date$Gene == "PA",])
summary(PA_N_model)



## Multi-intercept concept of 4. HA - N
AllG_NS_Date %>%
  filter(Gene == "HA") %>% 
  ggplot(aes(x = Regdate, y = Nmedian, color = Vaccine_code)) +
  geom_point(aes(fill = Vaccine_code), size = 2) +
  scale_x_continuous(name = "Year") +
  ylab("Number of Non-synonymous mutation of the HA gene") +
  scale_fill_manual(values = VacCoCol, name = "Vaccine Strains", guide="none") +
   scale_color_manual(values = VacCoCol, guide="none") +
  scale_y_continuous(breaks = seq(0, 35, by = 5)) +
  geom_hline(yintercept = 0) +
  geom_smooth(se = FALSE, method = lm) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9))

HA_N_model <- lmer(Nmedian ~ Regdate + (1 | Vaccine_code_Reg), data=AllG_NS_Date[AllG_NS_Date$Gene == "HA",])
summary(HA_N_model)

HA_N_model <- lmer(Nmedian ~ Regdate + (1 | Vaccine_code), data=AllG_NS_Date[AllG_NS_Date$Gene == "HA",])
summary(HA_N_model)



jpeg(filename = "/01_GenDisFlu/Fig/Desc/NS_Time/NP_N_Time_VaccG.jpeg", width = 180, height = 180, units = "mm", res = 400)

## Multi-intercept concept of 5. NP - N
AllG_NS_Date %>%
  filter(Gene == "NP") %>% 
  ggplot(aes(x = Regdate, y = Nmedian, color = Vaccine_code)) +
  geom_point(aes(fill = Vaccine_code), size = 2) +
  scale_x_continuous(name = "Year") +
  ylab("Number of Non-synonymous mutation of the NP gene") +
  scale_fill_manual(values = VacCoCol, name = "Vaccine Strains", guide="none") +
   scale_color_manual(values = VacCoCol, guide="none") +
  scale_y_continuous(breaks = seq(0, 35, by = 5)) +
  geom_hline(yintercept = 0) +
  geom_smooth(se = FALSE, method = lm) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9))

dev.off()

NP_N_model <- lmer(Nmedian ~ Regdate + (1 | Vaccine_code_Reg), data=AllG_NS_Date[AllG_NS_Date$Gene == "NP",])
summary(NP_N_model)

NP_N_model <- lmer(Nmedian ~ Regdate + (1 | Vaccine_code), data=AllG_NS_Date[AllG_NS_Date$Gene == "NP",])
summary(NP_N_model)

jpeg(filename = "/01_GenDisFlu/Fig/Desc/NS_Time/NA_N_Time_VaccG.jpeg", width = 180, height = 180, units = "mm", res = 400)

## Multi-intercept concept of 6. NA - N
AllG_NS_Date %>%
  filter(Gene == "NA") %>% 
  ggplot(aes(x = Regdate, y = Nmedian, color = Vaccine_code)) +
  geom_point(aes(fill = Vaccine_code), size = 2) +
  scale_x_continuous(name = "Year") +
  ylab("Number of Non-synonymous mutation of the NA gene") +
  scale_fill_manual(values = VacCoCol, name = "Vaccine Strains", guide="none") +
  scale_color_manual(values = VacCoCol, guide="none") +
  scale_y_continuous(breaks = seq(0, 35, by = 5)) +
  geom_hline(yintercept = 0) +
  geom_smooth(se = FALSE, method = lm) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9))
dev.off()

NA_N_model <- lmer(Nmedian ~ Regdate + (1 | Vaccine_code_Reg), data=AllG_NS_Date[AllG_NS_Date$Gene == "NA",])
summary(NA_N_model)

NA_N_model <- lmer(Nmedian ~ Regdate + (1 | Vaccine_code), data=AllG_NS_Date[AllG_NS_Date$Gene == "NA",])
summary(NA_N_model)

jpeg(filename = "/01_GenDisFlu/Fig/Desc/NS_Time/M_N_Time_VaccG.jpeg", width = 180, height = 180, units = "mm", res = 400)

## Multi-intercept concept of 7. M - N
AllG_NS_Date %>%
  filter(Gene == "M") %>% 
  ggplot(aes(x = Regdate, y = Nmedian, color = Vaccine_code)) +
  geom_point(aes(fill = Vaccine_code), size = 2) +
  scale_x_continuous(name = "Year") +
  ylab("Number of Non-synonymous mutation of the M gene") +
  scale_fill_manual(values = VacCoCol, name = "Vaccine Strains", guide="none") +
   scale_color_manual(values = VacCoCol, guide="none") +
  scale_y_continuous(breaks = seq(0, 35, by = 5)) +
  geom_hline(yintercept = 0) +
  geom_smooth(se = FALSE, method = lm) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9))
dev.off()

M_N_model <- lmer(Nmedian ~ Regdate + (1 | Vaccine_code_Reg), data=AllG_NS_Date[AllG_NS_Date$Gene == "M",])
summary(M_N_model)

M_N_model <- lmer(Nmedian ~ Regdate + (1 | Vaccine_code), data=AllG_NS_Date[AllG_NS_Date$Gene == "M",])
summary(M_N_model)


jpeg(filename = "/01_GenDisFlu/Fig/Desc/NS_Time/NS_N_Time_VaccG.jpeg", width = 180, height = 180, units = "mm", res = 400)

## Multi-intercept concept of 8. NS - N
AllG_NS_Date %>%
  filter(Gene == "NS") %>% 
  ggplot(aes(x = Regdate, y = Nmedian, color = Vaccine_code)) +
  geom_point(aes(fill = Vaccine_code), size = 2) +
  scale_x_continuous(name = "Year") +
  ylab("Number of Non-synonymous mutation of the NS gene") +
  scale_fill_manual(values = VacCoCol, name = "Vaccine Strains", guide="none") +
   scale_color_manual(values = VacCoCol, guide="none") +
  scale_y_continuous(breaks = seq(0, 35, by = 5)) +
  geom_hline(yintercept = 0) +
  geom_smooth(se = FALSE, method = lm) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9))
dev.off()


NS_N_model <- lmer(Nmedian ~ Regdate + (1 | Vaccine_code_Reg), data=AllG_NS_Date[AllG_NS_Date$Gene == "NS",])
summary(NS_N_model)


NS_N_model <- lmer(Nmedian ~ Regdate + (1 | Vaccine_code), data=AllG_NS_Date[AllG_NS_Date$Gene == "NS",])
summary(NS_N_model)
