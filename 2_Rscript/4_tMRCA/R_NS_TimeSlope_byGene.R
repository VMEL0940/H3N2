### Influ H3N2 GenDist Proj ###
## Package ##
library(tidyverse)
library(lme4)
## Set working directory ##
setwd("01_GenDisFlu/")

## Data import ##
## Data import ##
Meta <- read.csv("Data/AllData/H3N2_Traindata_600strains_1_Metadata.csv", header = T, na.strings = "")
NonSyn <- read.csv("Data/AllData/H3N2_Traindata_600strains_2_NonsynDist.csv", header = T, na.strings = "")
Syn <- read.csv("Data/AllData/H3N2_Traindata_600strains_3_SynDist.csv", header = T, na.strings = "")
Others <- read.csv("Data/AllData/H3N2_Traindata_600strains_4_OtherPredictors.csv", header = T, na.strings = "")

## Combine Data ##
All600 <- left_join(Meta, NonSyn[,c(-2,-3)], by = "ID")
All600 <- left_join(All600, Syn[,c(-2,-3)], by = "ID")
All600 <- left_join(All600, Others[,c(-2,-3)], by = "ID")

## Make decimal year column
# Create date object
All600$Date <- make_date(All600$Year, All600$Month, All600$Day)

# Compute decimal year
All600$DecimalYear <- All600$Year + (yday(All600$Date) - 1) / ifelse(leap_year(All600$Year), 366, 365)
All600$DecimalYear_0 <- All600$DecimalYear - 2000

## Calculation per amino acid
All600$PB2_Nonsyn <- All600$PB2_Nonsyn/761
All600$PB2_Syn <- All600$PB2_Syn/761
All600$PB1_Nonsyn <- All600$PB1_Nonsyn/758
All600$PB1_Syn <- All600$PB1_Syn/758
All600$PA_Nonsyn <- All600$PA_Nonsyn/717
All600$PA_Syn <- All600$PA_Syn/717
All600$HA_Nonsyn <- All600$HA_Nonsyn/567
All600$HA_Syn <- All600$HA_Syn/567
All600$NP_Nonsyn <- All600$NP_Nonsyn/499
All600$NP_Syn <- All600$NP_Syn/499
All600$NA_Nonsyn <- All600$NA_Nonsyn/470
All600$NA_Syn <- All600$NA_Syn/470
All600$M_Nonsyn <- All600$M_Nonsyn/328
All600$M_Syn <- All600$M_Syn/328
All600$NS_Nonsyn <- All600$NS_Nonsyn/280
All600$NS_Syn <- All600$NS_Syn/280


## Estimate the slope by Gene ##
### 1. PB2 NonSyn / Syn #### -- 761 AA
PB2_Nonsyn_year <- lmer(PB2_Nonsyn ~  DecimalYear_0+ (1 | vaccine_code), 
                      data = All600)

summary(PB2_Nonsyn_year)
fixef(PB2_Nonsyn_year)
confint(PB2_Nonsyn_year, level = 0.95)

PB2_Syn_year <- lmer(PB2_Syn ~ DecimalYear_0 + (1 | vaccine_code), 
                        data = All600)

summary(PB2_Syn_year)
fixef(PB2_Syn_year)
confint(PB2_Syn_year, level = 0.95)


### 2. PB1 NonSyn / Syn #### -- 758 AA

PB1_Nonsyn_year <- lmer(PB1_Nonsyn ~ DecimalYear_0 + (1 | vaccine_code), 
                        data = All600)

summary(PB1_Nonsyn_year)
fixef(PB1_Nonsyn_year)
confint(PB1_Nonsyn_year, level = 0.95)

PB1_Syn_year <- lmer(PB1_Syn ~ DecimalYear_0 + (1 | vaccine_code), 
                     data = All600)

summary(PB1_Syn_year)
fixef(PB1_Syn_year)
confint(PB1_Syn_year, level = 0.95)



### 3. PA NonSyn / Syn #### -- 717 AA

PA_Nonsyn_year <- lmer(PA_Nonsyn ~ DecimalYear_0 + (1 | vaccine_code), 
                        data = All600)

summary(PA_Nonsyn_year)
fixef(PA_Nonsyn_year)
confint(PA_Nonsyn_year, level = 0.95)

PA_Syn_year <- lmer(PA_Syn ~ DecimalYear_0 + (1 | vaccine_code), 
                     data = All600)

summary(PA_Syn_year)
fixef(PA_Syn_year)
confint(PA_Syn_year, level = 0.95)




### 4. HA NonSyn / Syn #### -- 567 AA

HA_Nonsyn_year <- lmer( HA_Nonsyn ~ DecimalYear_0 + (1 | vaccine_code), 
                       data = All600)

summary(HA_Nonsyn_year)
fixef(HA_Nonsyn_year)
confint(HA_Nonsyn_year, level = 0.95)

HA_Syn_year <- lmer(HA_Syn ~ DecimalYear_0 +  (1 | vaccine_code), 
                    data = All600)

summary(HA_Syn_year)
fixef(HA_Syn_year)
confint(HA_Syn_year, level = 0.95)



### 5. NP NonSyn / Syn #### -- 499 AA

NP_Nonsyn_year <- lmer(NP_Nonsyn ~ DecimalYear_0 +  (1 | vaccine_code), 
                       data = All600)

summary(NP_Nonsyn_year)
fixef(NP_Nonsyn_year)
confint(NP_Nonsyn_year, level = 0.95)

NP_Syn_year <- lmer(NP_Syn ~ DecimalYear_0 +  (1 | vaccine_code), 
                    data = All600)

summary(NP_Syn_year)
fixef(NP_Syn_year)
confint(NP_Syn_year, level = 0.95)



### 6. NA NonSyn / Syn #### -- 470 AA

NA_Nonsyn_year <- lmer(NA_Nonsyn ~ DecimalYear_0 +  (1 | vaccine_code), 
                       data = All600)

summary(NA_Nonsyn_year)
fixef(NA_Nonsyn_year)
confint(NA_Nonsyn_year, level = 0.95)

NA_Syn_year <- lmer(NA_Syn ~ DecimalYear_0 +  (1 | vaccine_code), 
                    data = All600)

summary(NA_Syn_year)
fixef(NA_Syn_year)
confint(NA_Syn_year, level = 0.95)




### 7. M NonSyn / Syn #### -- 328 AA

M_Nonsyn_year <- lmer(M_Nonsyn ~ DecimalYear_0 +  (1 | vaccine_code), 
                       data = All600)

summary(M_Nonsyn_year)
fixef(M_Nonsyn_year)
confint(M_Nonsyn_year, level = 0.95)

M_Syn_year <- lmer(M_Syn ~ DecimalYear_0 + (1 | vaccine_code), 
                    data = All600)

summary(M_Syn_year)
fixef(M_Syn_year)
confint(M_Syn_year, level = 0.95)


### 8. NS NonSyn / Syn #### -- 280 AA

NS_Nonsyn_year <- lmer(NS_Nonsyn ~ DecimalYear_0 +  (1 | vaccine_code), 
                       data = All600)

summary(NS_Nonsyn_year)
fixef(NS_Nonsyn_year)
confint(NS_Nonsyn_year, level = 0.95)

NS_Syn_year <- lmer(NS_Syn ~ DecimalYear_0 +  (1 | vaccine_code), 
                    data = All600)

summary(NS_Syn_year)
fixef(NS_Syn_year)
confint(NS_Syn_year, level = 0.95)


### Summarize in the Data frame
Slope_Summary <-data.frame()


Slope_Summary <- data.frame(Gene = c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS","PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"),
                            Mutation = c(rep("Nonsyn", 8), rep("Syn", 8)),
                            Slope = c(fixef(PB2_Nonsyn_year)[2], fixef(PB1_Nonsyn_year)[2], fixef(PA_Nonsyn_year)[2], fixef(HA_Nonsyn_year)[2],
                                      fixef(NP_Nonsyn_year)[2], fixef(NA_Nonsyn_year)[2], fixef(M_Nonsyn_year)[2], fixef(NS_Nonsyn_year)[2],
                                      fixef(PB2_Syn_year)[2], fixef(PB1_Syn_year)[2], fixef(PA_Syn_year)[2], fixef(HA_Syn_year)[2],
                                      fixef(NP_Syn_year)[2], fixef(NA_Syn_year)[2], fixef(M_Syn_year)[2], fixef(NS_Syn_year)[2]),
                            Low95 = c(confint(PB2_Nonsyn_year, level = 0.95)[4,1], confint(PB1_Nonsyn_year, level = 0.95)[4,1], confint(PA_Nonsyn_year, level = 0.95)[4,1], confint(HA_Nonsyn_year, level = 0.95)[4,1],
                                      confint(NP_Nonsyn_year, level = 0.95)[4,1], confint(NA_Nonsyn_year, level = 0.95)[4,1], confint(M_Nonsyn_year, level = 0.95)[4,1], confint(NS_Nonsyn_year, level = 0.95)[4,1],
                                      confint(PB2_Syn_year, level = 0.95)[4,1], confint(PB1_Syn_year, level = 0.95)[4,1], confint(PA_Syn_year, level = 0.95)[4,1], confint(HA_Syn_year, level = 0.95)[4,1],
                                      confint(NP_Syn_year, level = 0.95)[4,1], confint(NA_Syn_year, level = 0.95)[4,1], confint(M_Syn_year, level = 0.95)[4,1], confint(NS_Syn_year, level = 0.95)[4,1]),
                            Hgh95 = c(confint(PB2_Nonsyn_year, level = 0.95)[4,2], confint(PB1_Nonsyn_year, level = 0.95)[4,2], confint(PA_Nonsyn_year, level = 0.95)[4,2], confint(HA_Nonsyn_year, level = 0.95)[4,2],
                                      confint(NP_Nonsyn_year, level = 0.95)[4,2], confint(NA_Nonsyn_year, level = 0.95)[4,2], confint(M_Nonsyn_year, level = 0.95)[4,2], confint(NS_Nonsyn_year, level = 0.95)[4,2],
                                      confint(PB2_Syn_year, level = 0.95)[4,2], confint(PB1_Syn_year, level = 0.95)[4,2], confint(PA_Syn_year, level = 0.95)[4,2], confint(HA_Syn_year, level = 0.95)[4,2],
                                      confint(NP_Syn_year, level = 0.95)[4,2], confint(NA_Syn_year, level = 0.95)[4,2], confint(M_Syn_year, level = 0.95)[4,2], confint(NS_Syn_year, level = 0.95)[4,2]))


write.csv(Slope_Summary, "Slope_Summary.csv")



sort(Slope_Summary, Gene)
