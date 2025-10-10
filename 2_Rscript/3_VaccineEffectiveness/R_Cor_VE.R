### Influ H3N2 VaccineEffec vs HA N Ave ###
## Package ##
library(tidyverse)

## Set working directory ##
setwd("/01_GenDisFlu/")

## Data import ##
VE <- read.csv("Data/Vaccine Effectiveness/VE.csv", header = T, na.strings = "")
AllG_NS <- read.csv("Data/Genetic Distances/GeneCompete/AllG_NS.csv", header = T, na.strings = "")
Subsetlist <- read.csv("Data/H3N2_Index_Final_20220920.csv")

AllG_NS$compareStrain_V1 <- gsub("_NA_NA_NA_NA", "", AllG_NS$compareStrain)
AllG_NS$compareStrain_V1 <- gsub("_NA_NA_NA", "", AllG_NS$compareStrain_V1)
AllG_NS$compareStrain_V1 <- gsub("_NA_NA", "", AllG_NS$compareStrain_V1)
AllG_NS$compareStrain <- gsub("2009_NA", "2009", AllG_NS$compareStrain_V1)

## Extract HA gene 
NS_HA <- AllG_NS %>% 
  filter(Gene == "HA")

NS_HA_Date <- full_join(NS_HA, Subsetlist[,c(1,4:6)], by = "compareStrain")
NS_HA_Date <- NS_HA_Date[,c(-11)]

NS_HA_Date$Half <- ifelse(NS_HA_Date$Month >= 7, "ScdHalf", "FstHalf") 

VE$vaccineStrain

VE$HA_N <- 0
VE[VE$vaccineStrain == "EPI358781_A_Fujian_411_2002", ]$HA_N <- mean(NS_HA_Date[NS_HA_Date$vaccineStrain == "EPI358781_A_Fujian_411_2002", ]$Nmedian)
VE[VE$vaccineStrain == "EPI367109_A_California_7_2004", ]$HA_N <- mean(NS_HA_Date[NS_HA_Date$vaccineStrain == "EPI367109_A_California_7_2004", ]$Nmedian)  

VE[VE$vaccineStrain == "EPI502253_A_Wisconsin_67_2005" & VE$X2ndhalf == 2006, ]$HA_N <- 
  mean(c(NS_HA_Date[NS_HA_Date$vaccineStrain == "EPI502253_A_Wisconsin_67_2005" & NS_HA_Date$Year == 2006, ]$Nmedian, NS_HA_Date[NS_HA_Date$vaccineStrain == "EPI502253_A_Wisconsin_67_2005" & NS_HA_Date$Year == 2007 & NS_HA_Date$Half == "FstHalf", ]$Nmedian))
VE[VE$vaccineStrain == "EPI502253_A_Wisconsin_67_2005" & VE$X2ndhalf == 2007, ]$HA_N <-   
  mean(c(NS_HA_Date[NS_HA_Date$vaccineStrain == "EPI502253_A_Wisconsin_67_2005" & NS_HA_Date$Year == 2007 & NS_HA_Date$Half == "ScdHalf", ]$Nmedian, NS_HA_Date[NS_HA_Date$vaccineStrain == "EPI577980_A_Brisbane_10_2007" & NS_HA_Date$Year == 2008 & NS_HA_Date$Half == "FstHalf", ]$Nmedian))

VE[VE$vaccineStrain == "EPI577980_A_Brisbane_10_2007" & VE$X2ndhalf == 2008, ]$HA_N <- 
  mean(c(NS_HA_Date[NS_HA_Date$vaccineStrain == "EPI577980_A_Brisbane_10_2007" & NS_HA_Date$Year == 2008, ]$Nmedian, NS_HA_Date[NS_HA_Date$vaccineStrain == "EPI577980_A_Brisbane_10_2007" & NS_HA_Date$Year == 2009 & NS_HA_Date$Half == "FstHalf", ]$Nmedian))
VE[VE$vaccineStrain == "EPI577980_A_Brisbane_10_2007" & VE$X2ndhalf == 2009, ]$HA_N <-   
  mean(c(NS_HA_Date[NS_HA_Date$vaccineStrain == "EPI577980_A_Brisbane_10_2007" & NS_HA_Date$Year == 2009 & NS_HA_Date$Half == "ScdHalf", ]$Nmedian, NS_HA_Date[NS_HA_Date$vaccineStrain == "EPI577980_A_Brisbane_10_2007" & NS_HA_Date$Year == 2010 & NS_HA_Date$Half == "FstHalf", ]$Nmedian))

VE[VE$vaccineStrain == "EPI577969_A_Perth_16_2009" & VE$X2ndhalf == 2010, ]$HA_N <- 
  mean(c(NS_HA_Date[NS_HA_Date$vaccineStrain == "EPI577969_A_Perth_16_2009" & NS_HA_Date$Year == 2010, ]$Nmedian, NS_HA_Date[NS_HA_Date$vaccineStrain == "EPI577969_A_Perth_16_2009" & NS_HA_Date$Year == 2011 & NS_HA_Date$Half == "FstHalf", ]$Nmedian))
VE[VE$vaccineStrain == "EPI577969_A_Perth_16_2009" & VE$X2ndhalf == 2011, ]$HA_N <-   
  mean(c(NS_HA_Date[NS_HA_Date$vaccineStrain == "EPI577969_A_Perth_16_2009" & NS_HA_Date$Year == 2011 & NS_HA_Date$Half == "ScdHalf", ]$Nmedian, NS_HA_Date[NS_HA_Date$vaccineStrain == "EPI577969_A_Perth_16_2009" & NS_HA_Date$Year == 2012 & NS_HA_Date$Half == "FstHalf", ]$Nmedian))

VE[VE$vaccineStrain == "EPI417234_A_Victoria_361_2011" & VE$X2ndhalf == 2012, ]$HA_N <- 
  mean(c(NS_HA_Date[NS_HA_Date$vaccineStrain == "EPI417234_A_Victoria_361_2011" & NS_HA_Date$Year == 2012, ]$Nmedian, NS_HA_Date[NS_HA_Date$vaccineStrain == "EPI417234_A_Victoria_361_2011" & NS_HA_Date$Year == 2013 & NS_HA_Date$Half == "FstHalf", ]$Nmedian))
VE[VE$vaccineStrain == "EPI417234_A_Victoria_361_2011" & VE$X2ndhalf == 2013, ]$HA_N <-   
  mean(c(NS_HA_Date[NS_HA_Date$vaccineStrain == "EPI417234_A_Victoria_361_2011" & NS_HA_Date$Year == 2013 & NS_HA_Date$Half == "ScdHalf", ]$Nmedian, NS_HA_Date[NS_HA_Date$vaccineStrain == "EPI417234_A_Victoria_361_2011" & NS_HA_Date$Year == 2014 & NS_HA_Date$Half == "FstHalf", ]$Nmedian))
VE[VE$vaccineStrain == "EPI417234_A_Victoria_361_2011" & VE$X2ndhalf == 2014, ]$HA_N <-   
  mean(c(NS_HA_Date[NS_HA_Date$vaccineStrain == "EPI417234_A_Victoria_361_2011" & NS_HA_Date$Year == 2014 & NS_HA_Date$Half == "ScdHalf", ]$Nmedian, NS_HA_Date[NS_HA_Date$vaccineStrain == "EPI417234_A_Victoria_361_2011" & NS_HA_Date$Year == 2015 & NS_HA_Date$Half == "FstHalf", ]$Nmedian))

VE[VE$vaccineStrain == "EPI614441_A_Switzerland_9715293_2013", ]$HA_N <- mean(NS_HA_Date[NS_HA_Date$vaccineStrain == "EPI614441_A_Switzerland_9715293_2013", ]$Nmedian)  

VE[VE$vaccineStrain == "EPI686117_A_Hong_Kong_15611_2015" & VE$X2ndhalf == 2016, ]$HA_N <- 
  mean(c(NS_HA_Date[NS_HA_Date$vaccineStrain == "EPI686117_A_Hong_Kong_15611_2015" & NS_HA_Date$Year == 2016, ]$Nmedian, NS_HA_Date[NS_HA_Date$vaccineStrain == "EPI686117_A_Hong_Kong_15611_2015" & NS_HA_Date$Year == 2017 & NS_HA_Date$Half == "FstHalf", ]$Nmedian))
VE[VE$vaccineStrain == "EPI686117_A_Hong_Kong_15611_2015" & VE$X2ndhalf == 2017, ]$HA_N <-   
  mean(c(NS_HA_Date[NS_HA_Date$vaccineStrain == "EPI686117_A_Hong_Kong_15611_2015" & NS_HA_Date$Year == 2017 & NS_HA_Date$Half == "ScdHalf", ]$Nmedian, NS_HA_Date[NS_HA_Date$vaccineStrain == "EPI686117_A_Hong_Kong_15611_2015" & NS_HA_Date$Year == 2018 & NS_HA_Date$Half == "FstHalf", ]$Nmedian))
VE[VE$vaccineStrain == "EPI686117_A_Hong_Kong_15611_2015" & VE$X2ndhalf == 2018, ]$HA_N <-   
  mean(c(NS_HA_Date[NS_HA_Date$vaccineStrain == "EPI686117_A_Hong_Kong_15611_2015" & NS_HA_Date$Year == 2018 & NS_HA_Date$Half == "ScdHalf", ]$Nmedian, NS_HA_Date[NS_HA_Date$vaccineStrain == "EPI686117_A_Hong_Kong_15611_2015" & NS_HA_Date$Year == 2019 & NS_HA_Date$Half == "FstHalf", ]$Nmedian))

cor.test(x= VE$VE, y= VE$HA_N, method = c("pearson")) 

## run linear regression with fixed intercept = 0
Nave_VE <- lm(VE ~ HA_N, data = VE)
summary(Nave_VE) ## Summary

Nave_VE$coefficients[2] 

jpeg(filename = "/01_GenDisFlu/Presentation/Fig/VE/VE_HA.jpeg", width = 150, height = 200, units = "mm", res = 420)

## Jitter plot with Linear line
VE %>%
  ggplot(aes(x=HA_N, y=VE)) +
  geom_point(size=1.2) +
  geom_smooth(method = "lm", color='red', se=FALSE ) +
  xlab("Number of Non-synonymous Mutation in the HA gene") +
  ylab("Vaccine Effectiveness") +
#  geom_hline(yintercept = 0, size = 0.6) +
#  geom_vline(xintercept = 0, size = 0.7) +
  scale_x_continuous(breaks=seq(0, 30, 5)) +
  theme_bw()

dev.off()
