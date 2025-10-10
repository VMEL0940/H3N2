### Influ H3N2 GenDist Proj ###
## Package ##
library(tidyverse)

## Set working directory ##
setwd("01_GenDisFlu/")

## Data import ##
AllG_NS <- read.csv("Data/Genetic Distances/GeneCompete/AllG_NS.csv", header = T, na.strings = "")
Subsetlist <- read.csv("Data/H3N2_Index_20230331.csv")

AllG_NS$Gene <- factor(AllG_NS$Gene, levels = c("PB2", "PB1", "PA", "HA", 
                                            "NP", "NA", "M", "NS"))
AllG_NS$Vaccine_code <- factor(AllG_NS$Vaccine_code,levels = c("1.Mos99","2.Fuj02","3.Cal04",
                                                               "4.Wis05","5.Bris07","6.Prth09",
                                                               "7.Vic11", "8.Swtz13", "9.HK15"),
                                                        labels = c("1.Mos99","2.Fuj02","3.Cal04",
                                                                   "4.Wis05","5.Bris07","6.Prth09",
                                                                   "7.Vic11", "8.Swtz13", "9.HK15"))

VacCoCol <- c("#716172", "#EFA95F","#3EA3C5","#6D9431","#800000","#9E6D13","#092E50","#769F74","#954495")


AllG_NS$compareStrain_V1 <- gsub("_NA_NA_NA_NA", "", AllG_NS$compareStrain)
AllG_NS$compareStrain_V1 <- gsub("_NA_NA_NA", "", AllG_NS$compareStrain_V1)
AllG_NS$compareStrain_V1 <- gsub("_NA_NA", "", AllG_NS$compareStrain_V1)
AllG_NS$compareStrain <- gsub("2009_NA", "2009", AllG_NS$compareStrain_V1)
AllG_NS <- AllG_NS[,-12]

AllG_NS_Date <- full_join(AllG_NS, Subsetlist[,c(1,4:6)], by = "compareStrain")
AllG_NS_Date$ColDate <- paste0(AllG_NS_Date$Year, "/", AllG_NS_Date$Month, "/", AllG_NS_Date$Day)
AllG_NS_Date$ColDate <- as.Date(AllG_NS_Date$ColDate,"%Y/%m/%d")

##### 1. HA gene ##### 
 
jpeg(filename = "01_GenDisFlu/Fig/Desc/NS_Time/HA_N_Time_Lged.jpeg", width = 250, height = 130, units = "mm", res = 800)

AllG_NS_Date %>%
  filter(Gene == "HA") %>% 
  ggplot() +
  geom_point(aes(x = ColDate, y = Nmedian, fill = Vaccine_code), size = 2.5, shape = 21, color = "gray20") +
  geom_vline(xintercept = as.Date(c("2004-10-1", "2005-10-1",
                                    "2006-10-1", "2008-10-1",
                                    "2010-10-1", "2012-10-1", 
                                    "2015-10-1", "2016-10-1")), 
             color = "red", lty = 2, 
             lwd = 0.4) +
  scale_x_date(name = "Year", minor_breaks = "1 year") +
  ylab("Nonsynonymous genetic distance of the HA gene") +
  scale_fill_manual(values = VacCoCol, name = "Vaccine Period") +
  scale_y_continuous(breaks = seq(0, 45, by = 5), limits = c(0, 35)) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9))

dev.off()
  
jpeg(filename = "01_GenDisFlu/Fig/Desc/NS_Time/HA_S_Time_Lged.jpeg", width = 250, height = 130, units = "mm", res = 800)

AllG_NS_Date %>%
  filter(Gene == "HA") %>% 
  ggplot() +
  geom_point(aes(x = ColDate, y = Smedian, fill = Vaccine_code), size = 2.5, shape = 21, color = "gray20") +
  geom_vline(xintercept = as.Date(c("2004-10-1", "2005-10-1",
                                    "2006-10-1", "2008-10-1",
                                    "2010-10-1", "2012-10-1", 
                                    "2015-10-1", "2016-10-1")), 
             color = "red", lty = 2, 
             lwd = 0.4) +
  scale_x_date(name = "Year", minor_breaks = "1 year") +
  ylab("Synonymous genetic distance of the HA gene") +
  scale_fill_manual(values = VacCoCol, name = "Vaccine Period") +
  scale_y_continuous(breaks = seq(0, 45, by = 5), limits = c(0, 35)) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9))

dev.off()


##### 2. PB2 gene ##### 

jpeg(filename = "/01_GenDisFlu/Fig/Desc/NS_Time/PB2_N_Time.jpeg", width = 300, height = 150, units = "mm", res = 520)

AllG_NS_Date %>%
  filter(Gene == "PB2") %>% 
  ggplot() +
  geom_point(aes(x = ColDate, y = Nmedian, fill = Vaccine_code), size = 2.5, shape = 21, color = "gray20") +
  geom_vline(xintercept = as.Date(c("2004-10-1", "2005-10-1",
                                    "2006-10-1", "2008-10-1",
                                    "2010-10-1", "2012-10-1", 
                                    "2015-10-1", "2016-10-1")), 
             color = "red", lty = 2, 
             lwd = 0.4) +
  scale_x_date(name = "Year", minor_breaks = "1 year") +
  ylab("Renaissance count of Non-synonymous mutation in the PB2 gene") +
  scale_fill_manual(values = VacCoCol, name = "Vaccine Strains") +
  scale_y_continuous(breaks = seq(0, 35, by = 5)) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9))

dev.off()


jpeg(filename = "/01_GenDisFlu/Fig/Desc/NS_Time/PB2_S_Time.jpeg", width = 300, height = 150, units = "mm", res = 520)

AllG_NS_Date %>%
  filter(Gene == "PB2") %>% 
  ggplot() +
  geom_point(aes(x = ColDate, y = Smedian, fill = Vaccine_code), size = 2.5, shape = 21, color = "gray20") +
  geom_vline(xintercept = as.Date(c("2004-10-1", "2005-10-1",
                                    "2006-10-1", "2008-10-1",
                                    "2010-10-1", "2012-10-1", 
                                    "2015-10-1", "2016-10-1")), 
             color = "red", lty = 2, 
             lwd = 0.4) +
  scale_x_date(name = "Year", minor_breaks = "1 year") +
  ylab("Renaissance count of Synonymous mutation of the PB2 gene") +
  scale_fill_manual(values = VacCoCol, name = "Vaccine Strains") +
  scale_y_continuous(breaks = seq(0, 45, by = 5)) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9))

dev.off()



##### 3. PB1 gene ##### 


jpeg(filename = "/01_GenDisFlu/Fig/Desc/NS_Time/PB1_N_Time.jpeg", width = 300, height = 150, units = "mm", res = 520)

AllG_NS_Date %>%
  filter(Gene == "PB1") %>% 
  ggplot() +
  geom_point(aes(x = ColDate, y = Nmedian, fill = Vaccine_code), size = 2.5, shape = 21, color = "gray20") +
  geom_vline(xintercept = as.Date(c("2004-10-1", "2005-10-1",
                                    "2006-10-1", "2008-10-1",
                                    "2010-10-1", "2012-10-1", 
                                    "2015-10-1", "2016-10-1")), 
             color = "red", lty = 2, 
             lwd = 0.4) +
  scale_x_date(name = "Year", minor_breaks = "1 year") +
  ylab("Renaissance count of Non-synonymous mutation in the PB1 gene") +
  scale_fill_manual(values = VacCoCol, name = "Vaccine Strains") +
  scale_y_continuous(breaks = seq(0, 35, by = 5)) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9))

dev.off()


jpeg(filename = "/01_GenDisFlu/Fig/Desc/NS_Time/PB1_S_Time.jpeg", width = 300, height = 150, units = "mm", res = 520)

AllG_NS_Date %>%
  filter(Gene == "PB1") %>% 
  ggplot() +
  geom_point(aes(x = ColDate, y = Smedian, fill = Vaccine_code), size = 2.5, shape = 21, color = "gray20") +
  geom_vline(xintercept = as.Date(c("2004-10-1", "2005-10-1",
                                    "2006-10-1", "2008-10-1",
                                    "2010-10-1", "2012-10-1", 
                                    "2015-10-1", "2016-10-1")), 
             color = "red", lty = 2, 
             lwd = 0.4) +
  scale_x_date(name = "Year", minor_breaks = "1 year") +
  ylab("Renaissance count of Synonymous mutation of the PB1 gene") +
  scale_fill_manual(values = VacCoCol, name = "Vaccine Strains") +
  scale_y_continuous(breaks = seq(0, 45, by = 5)) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9))

dev.off()



##### 4. PA gene ##### 


jpeg(filename = "/01_GenDisFlu/Fig/Desc/NS_Time/PA_N_Time.jpeg", width = 300, height = 150, units = "mm", res = 520)

AllG_NS_Date %>%
  filter(Gene == "PA") %>% 
  ggplot() +
  geom_point(aes(x = ColDate, y = Nmedian, fill = Vaccine_code), size = 2.5, shape = 21, color = "gray20") +
  geom_vline(xintercept = as.Date(c("2004-10-1", "2005-10-1",
                                    "2006-10-1", "2008-10-1",
                                    "2010-10-1", "2012-10-1", 
                                    "2015-10-1", "2016-10-1")), 
             color = "red", lty = 2, 
             lwd = 0.4) +
  scale_x_date(name = "Year", minor_breaks = "1 year") +
  ylab("Renaissance count of Non-synonymous mutation in the PA gene") +
  scale_fill_manual(values = VacCoCol, name = "Vaccine Strains") +
  scale_y_continuous(breaks = seq(0, 35, by = 5)) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9))

dev.off()


jpeg(filename = "/01_GenDisFlu/Fig/Desc/NS_Time/PA_S_Time.jpeg", width = 300, height = 150, units = "mm", res = 520)

AllG_NS_Date %>%
  filter(Gene == "PA") %>% 
  ggplot() +
  geom_point(aes(x = ColDate, y = Smedian, fill = Vaccine_code), size = 2.5, shape = 21, color = "gray20") +
  geom_vline(xintercept = as.Date(c("2004-10-1", "2005-10-1",
                                    "2006-10-1", "2008-10-1",
                                    "2010-10-1", "2012-10-1", 
                                    "2015-10-1", "2016-10-1")), 
             color = "red", lty = 2, 
             lwd = 0.4) +
  scale_x_date(name = "Year", minor_breaks = "1 year") +
  ylab("Renaissance count of Synonymous mutation of the PA gene") +
  scale_fill_manual(values = VacCoCol, name = "Vaccine Strains") +
  scale_y_continuous(breaks = seq(0, 45, by = 5)) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9))

dev.off()



##### 5. NP gene ##### 


jpeg(filename = "/01_GenDisFlu/Fig/Desc/NS_Time/NP_N_Time.jpeg", width = 300, height = 150, units = "mm", res = 520)

AllG_NS_Date %>%
  filter(Gene == "NP") %>% 
  ggplot() +
  geom_point(aes(x = ColDate, y = Nmedian, fill = Vaccine_code), size = 2.5, shape = 21, color = "gray20") +
  geom_vline(xintercept = as.Date(c("2004-10-1", "2005-10-1",
                                    "2006-10-1", "2008-10-1",
                                    "2010-10-1", "2012-10-1", 
                                    "2015-10-1", "2016-10-1")), 
             color = "red", lty = 2, 
             lwd = 0.4) +
  scale_x_date(name = "Year", minor_breaks = "1 year") +
  ylab("Renaissance count of Non-synonymous mutation in the NP gene") +
  scale_fill_manual(values = VacCoCol, name = "Vaccine Strains") +
  scale_y_continuous(breaks = seq(0, 35, by = 5)) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9))

dev.off()


jpeg(filename = "/01_GenDisFlu/Fig/Desc/NS_Time/NP_S_Time.jpeg", width = 300, height = 150, units = "mm", res = 520)

AllG_NS_Date %>%
  filter(Gene == "NP") %>% 
  ggplot() +
  geom_point(aes(x = ColDate, y = Smedian, fill = Vaccine_code), size = 2.5, shape = 21, color = "gray20") +
  geom_vline(xintercept = as.Date(c("2004-10-1", "2005-10-1",
                                    "2006-10-1", "2008-10-1",
                                    "2010-10-1", "2012-10-1", 
                                    "2015-10-1", "2016-10-1")), 
             color = "red", lty = 2, 
             lwd = 0.4) +
  scale_x_date(name = "Year", minor_breaks = "1 year") +
  ylab("Renaissance count of Synonymous mutation of the NP gene") +
  scale_fill_manual(values = VacCoCol, name = "Vaccine Strains") +
  scale_y_continuous(breaks = seq(0, 45, by = 5)) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9))

dev.off()



##### 6. NA gene ##### 


jpeg(filename = "/01_GenDisFlu/Fig/Desc/NS_Time/NA_N_Time.jpeg", width = 300, height = 150, units = "mm", res = 520)

AllG_NS_Date %>%
  filter(Gene == "NA") %>% 
  ggplot() +
  geom_point(aes(x = ColDate, y = Nmedian, fill = Vaccine_code), size = 2.5, shape = 21, color = "gray20") +
  geom_vline(xintercept = as.Date(c("2004-10-1", "2005-10-1",
                                    "2006-10-1", "2008-10-1",
                                    "2010-10-1", "2012-10-1", 
                                    "2015-10-1", "2016-10-1")), 
             color = "red", lty = 2, 
             lwd = 0.4) +
  scale_x_date(name = "Year", minor_breaks = "1 year") +
  ylab("Renaissance count of Non-synonymous mutation in the NA gene") +
  scale_fill_manual(values = VacCoCol, name = "Vaccine Strains") +
  scale_y_continuous(breaks = seq(0, 35, by = 5)) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9))

dev.off()


jpeg(filename = "/01_GenDisFlu/Fig/Desc/NS_Time/NA_S_Time.jpeg", width = 300, height = 150, units = "mm", res = 520)

AllG_NS_Date %>%
  filter(Gene == "NA") %>% 
  ggplot() +
  geom_point(aes(x = ColDate, y = Smedian, fill = Vaccine_code), size = 2.5, shape = 21, color = "gray20") +
  geom_vline(xintercept = as.Date(c("2004-10-1", "2005-10-1",
                                    "2006-10-1", "2008-10-1",
                                    "2010-10-1", "2012-10-1", 
                                    "2015-10-1", "2016-10-1")), 
             color = "red", lty = 2, 
             lwd = 0.4) +
  scale_x_date(name = "Year", minor_breaks = "1 year") +
  ylab("Renaissance count of Synonymous mutation of the NA gene") +
  scale_fill_manual(values = VacCoCol, name = "Vaccine Strains") +
  scale_y_continuous(breaks = seq(0, 45, by = 5)) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9))

dev.off()



##### 7. M gene ##### 


jpeg(filename = "/01_GenDisFlu/Fig/Desc/NS_Time/M_N_Time.jpeg", width = 300, height = 150, units = "mm", res = 520)

AllG_NS_Date %>%
  filter(Gene == "M") %>% 
  ggplot() +
  geom_point(aes(x = ColDate, y = Nmedian, fill = Vaccine_code), size = 2.5, shape = 21, color = "gray20") +
  geom_vline(xintercept = as.Date(c("2004-10-1", "2005-10-1",
                                    "2006-10-1", "2008-10-1",
                                    "2010-10-1", "2012-10-1", 
                                    "2015-10-1", "2016-10-1")), 
             color = "red", lty = 2, 
             lwd = 0.4) +
  scale_x_date(name = "Year", minor_breaks = "1 year") +
  ylab("Renaissance count of Non-synonymous mutation in the M gene") +
  scale_fill_manual(values = VacCoCol, name = "Vaccine Strains") +
  scale_y_continuous(breaks = seq(0, 35, by = 5)) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9))

dev.off()


jpeg(filename = "/01_GenDisFlu/Fig/Desc/NS_Time/M_S_Time.jpeg", width = 300, height = 150, units = "mm", res = 520)

AllG_NS_Date %>%
  filter(Gene == "M") %>% 
  ggplot() +
  geom_point(aes(x = ColDate, y = Smedian, fill = Vaccine_code), size = 2.5, shape = 21, color = "gray20") +
  geom_vline(xintercept = as.Date(c("2004-10-1", "2005-10-1",
                                    "2006-10-1", "2008-10-1",
                                    "2010-10-1", "2012-10-1", 
                                    "2015-10-1", "2016-10-1")), 
             color = "red", lty = 2, 
             lwd = 0.4) +
  scale_x_date(name = "Year", minor_breaks = "1 year") +
  ylab("Renaissance count of Synonymous mutation of the M gene") +
  scale_fill_manual(values = VacCoCol, name = "Vaccine Strains") +
  scale_y_continuous(breaks = seq(0, 45, by = 5)) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9))

dev.off()



##### 8. NS gene ##### 


jpeg(filename = "/01_GenDisFlu/Fig/Desc/NS_Time/NS_N_Time.jpeg", width = 300, height = 150, units = "mm", res = 520)

AllG_NS_Date %>%
  filter(Gene == "NS") %>% 
  ggplot() +
  geom_point(aes(x = ColDate, y = Nmedian, fill = Vaccine_code), size = 2.5, shape = 21, color = "gray20") +
  geom_vline(xintercept = as.Date(c("2004-10-1", "2005-10-1",
                                    "2006-10-1", "2008-10-1",
                                    "2010-10-1", "2012-10-1", 
                                    "2015-10-1", "2016-10-1")), 
             color = "red", lty = 2, 
             lwd = 0.4) +
  scale_x_date(name = "Year", minor_breaks = "1 year") +
  ylab("Renaissance count of Non-synonymous mutation in the NS gene") +
  scale_fill_manual(values = VacCoCol, name = "Vaccine Strains") +
  scale_y_continuous(breaks = seq(0, 35, by = 5)) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9))

dev.off()


jpeg(filename = "/01_GenDisFlu/Fig/Desc/NS_Time/NS_S_Time.jpeg", width = 300, height = 150, units = "mm", res = 520)

AllG_NS_Date %>%
  filter(Gene == "NS") %>% 
  ggplot() +
  geom_point(aes(x = ColDate, y = Smedian, fill = Vaccine_code), size = 2.5, shape = 21, color = "gray20") +
  geom_vline(xintercept = as.Date(c("2004-10-1", "2005-10-1",
                                    "2006-10-1", "2008-10-1",
                                    "2010-10-1", "2012-10-1", 
                                    "2015-10-1", "2016-10-1")), 
             color = "red", lty = 2, 
             lwd = 0.4) +
  scale_x_date(name = "Year", minor_breaks = "1 year") +
  ylab("Renaissance count of Synonymous mutation of the NS gene") +
  scale_fill_manual(values = VacCoCol, name = "Vaccine Strains") +
  scale_y_continuous(breaks = seq(0, 45, by = 5)) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9))

dev.off()



### Gene size vs Syn distance
SynMxbyGene <- AllG_NS %>% 
  group_by(Gene) %>% 
  summarize(SynMax = max(Smedian), NSynMax = max(Nmedian))

SynMxbyGene$ntlength <- c(2283, 2274, 2151, 1701, 1497, 1410, 982, 838)
SynMxbyGene$SynLenRatio <- SynMxbyGene$ntlength/SynMxbyGene$SynMax
SynMxbyGene$NSynLenRatio <- SynMxbyGene$ntlength/SynMxbyGene$NSynMax
