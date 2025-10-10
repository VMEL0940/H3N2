## Package ##
library(tidyverse)
library(lme4)
library(MuMIn)
library(pROC)
library(lmtest)
library(corrplot)

## Set working directory ##
setwd("01_GenDisFlu/")
#setwd("01_GenDisFlu/") ## for mac
options(na.action = "na.fail")

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

## 1. Backward model selection for the Train Data 1 ##
### A. Nonsyn Test - Fix HA ##
full_model1 <- glmer(Dom ~ PB2_Nonsyn + PB1_Nonsyn + PA_Nonsyn + HA_Nonsyn + NP_Nonsyn + NA_Nonsyn + M_Nonsyn + NS_Nonsyn + (1 | vaccine_code), 
                    data = All600_Train, family = binomial)

dredge(full_model1, fixed = ~ HA_Nonsyn) ## Fix HA 

## Variables in Top 3 models -- HA, NA, M and NS

### B. Syn Test - Fix HA ##
full_model2 <- glmer(Dom ~ HA_Nonsyn + PB2_Nonsyn + PB1_Syn + PA_Syn + HA_Syn + NP_Syn + NA_Syn + M_Syn + NS_Syn + (1 | vaccine_code), 
                     data = All600_Train, family = binomial)

dredge(full_model2, fixed = ~ HA_Nonsyn)

## Variables in Top 3 models -- HA Nonsyn, HA Syn, M Syn, NP Syn,

### C. Combine NonSyn and Syn -  Fix HA Nonsyn##
full_model3 <- glmer(Dom ~ HA_Nonsyn + NA_Nonsyn + M_Nonsyn + NS_Nonsyn + HA_Syn + NP_Syn + M_Syn + (1 | vaccine_code), 
                     data = All600_Train, family = binomial)

dredge(full_model3, fixed = ~ HA_Nonsyn)

## Final variables - HA Nonsyn, NA Nonsyn, NS Nonsyn, M Syn,


### D. Forward selection ###
#### Candidates of Step 1 ####
model1_1 <- glmer(Dom ~ HA_Nonsyn + NA_Nonsyn + M_Nonsyn + (1 | vaccine_code), 
                  data = All600_Train, family = binomial)

model1_2 <- glmer(Dom ~ HA_Nonsyn + NA_Nonsyn + NS_Nonsyn + M_Nonsyn + (1 | vaccine_code), 
                  data = All600_Train, family = binomial)

model1_3 <- glmer(Dom ~ HA_Nonsyn + NA_Nonsyn + (1 | vaccine_code), 
                  data = All600_Train, family = binomial)

#### Candidates of Step 2 ####
model2_1 <- glmer(Dom ~ HA_Nonsyn + HA_Syn + M_Syn + (1 | vaccine_code), 
                  data = All600_Train, family = binomial)

model2_2 <- glmer(Dom ~ HA_Nonsyn + M_Syn + NP_Syn + (1 | vaccine_code), 
                  data = All600_Train, family = binomial)

model2_3 <- glmer(Dom ~ HA_Nonsyn + HA_Syn + M_Syn + NP_Syn + (1 | vaccine_code), 
                  data = All600_Train, family = binomial)

#### Candidates of Step 3 ####
model3_1 <- glmer(Dom ~ HA_Nonsyn + NA_Nonsyn + NS_Nonsyn + M_Syn + (1 | vaccine_code), 
                  data = All600_Train, family = binomial)

model3_2 <- glmer(Dom ~ HA_Nonsyn + NA_Nonsyn + NS_Nonsyn + M_Syn + M_Nonsyn+ (1 | vaccine_code), 
                  data = All600_Train, family = binomial)

model3_3 <- glmer(Dom ~ HA_Nonsyn + NA_Nonsyn + NS_Nonsyn + M_Syn + HA_Syn + (1 | vaccine_code), 
                  data = All600_Train, family = binomial)

#### Candidates of Step 4 ####
model4_1 <- glmer(Dom ~ HA_Nonsyn * NA_Nonsyn + NS_Nonsyn + M_Syn + (1 | vaccine_code), 
                  data = All600_Train, family = binomial)

dredge(model4_1, fixed = ~ HA_Nonsyn)

model4_2 <- glmer(Dom ~ HA_Nonsyn * NS_Nonsyn + NA_Nonsyn + M_Syn + (1 | vaccine_code), 
                  data = All600_Train, family = binomial)

dredge(model4_2, fixed = ~ HA_Nonsyn)

model4_3 <- glmer(Dom ~ HA_Nonsyn * M_Syn + NA_Nonsyn + NS_Nonsyn + (1 | vaccine_code), 
                  data = All600_Train, family = binomial)

dredge(model4_3, fixed = ~ HA_Nonsyn)

model4_4 <- glmer(Dom ~ HA_Nonsyn + NA_Nonsyn * NS_Nonsyn + M_Syn + (1 | vaccine_code), 
                  data = All600_Train, family = binomial)

dredge(model4_4, fixed = ~ HA_Nonsyn)

model4_5 <- glmer(Dom ~ HA_Nonsyn + NA_Nonsyn * M_Syn + NS_Nonsyn + (1 | vaccine_code), 
                  data = All600_Train, family = binomial)

dredge(model4_5, fixed = ~ HA_Nonsyn)

model4_6 <- glmer(Dom ~ HA_Nonsyn + NA_Nonsyn + NS_Nonsyn * M_Syn + (1 | vaccine_code), 
                  data = All600_Train, family = binomial)

dredge(model4_6, fixed = ~ HA_Nonsyn)


## Three variable models ##

model4_7 <- glmer(Dom ~ HA_Nonsyn * NA_Nonsyn + NS_Nonsyn  + (1 | vaccine_code), 
                  data = All600_Train, family = binomial)

model4_8 <- glmer(Dom ~ HA_Nonsyn + NA_Nonsyn * NS_Nonsyn  + (1 | vaccine_code), 
                  data = All600_Train, family = binomial)

# Extract fixed effects
coef_summary <- summary(model4_8)$coefficients
# Convert coefficients to odds ratios
odds_ratios <- exp(coef_summary[, "Estimate"])


model4_9 <- glmer(Dom ~ HA_Nonsyn * NS_Nonsyn + NA_Nonsyn + (1 | vaccine_code), 
                  data = All600_Train, family = binomial)

model4_10 <- glmer(Dom ~ HA_Nonsyn * NA_Nonsyn + M_Syn  + (1 | vaccine_code), 
                   data = All600_Train, family = binomial)

model4_11 <- glmer(Dom ~ HA_Nonsyn + NA_Nonsyn * M_Syn  + (1 | vaccine_code), 
                  data = All600_Train, family = binomial)

model4_12 <- glmer(Dom ~ HA_Nonsyn * M_Syn + NA_Nonsyn + (1 | vaccine_code), 
                   data = All600_Train, family = binomial)

model4_13 <- glmer(Dom ~ HA_Nonsyn * NS_Nonsyn + M_Syn  + (1 | vaccine_code), 
                   data = All600_Train, family = binomial)

model4_14 <- glmer(Dom ~ HA_Nonsyn + NS_Nonsyn * M_Syn  + (1 | vaccine_code), 
                   data = All600_Train, family = binomial)

model4_15 <- glmer(Dom ~ HA_Nonsyn * M_Syn + NS_Nonsyn + (1 | vaccine_code), 
                   data = All600_Train, family = binomial)



#### A. Predicted Score ####
Final_model1_1_logit_P = predict(model1_1, newdata = All600_Train, type = 'response' )
Final_model1_2_logit_P = predict(model1_2, newdata = All600_Train, type = 'response' )
Final_model1_3_logit_P = predict(model1_3, newdata = All600_Train, type = 'response' )
Final_model2_1_logit_P = predict(model2_1, newdata = All600_Train, type = 'response' )
Final_model2_2_logit_P = predict(model2_2, newdata = All600_Train, type = 'response' )
Final_model2_3_logit_P = predict(model2_3, newdata = All600_Train, type = 'response' )
Final_model3_1_logit_P = predict(model3_1, newdata = All600_Train, type = 'response' )
Final_model3_2_logit_P = predict(model3_2, newdata = All600_Train, type = 'response' )
Final_model3_3_logit_P = predict(model3_3, newdata = All600_Train, type = 'response' )

Final_model4_1_logit_P = predict(model4_1, newdata = All600_Train, type = 'response' )
Final_model4_2_logit_P = predict(model4_2, newdata = All600_Train, type = 'response' )
Final_model4_3_logit_P = predict(model4_3, newdata = All600_Train, type = 'response' )
Final_model4_4_logit_P = predict(model4_4, newdata = All600_Train, type = 'response' )
Final_model4_5_logit_P = predict(model4_5, newdata = All600_Train, type = 'response' )
Final_model4_6_logit_P = predict(model4_6, newdata = All600_Train, type = 'response' )
Final_model4_7_logit_P = predict(model4_7, newdata = All600_Train, type = 'response' )
Final_model4_8_logit_P = predict(model4_8, newdata = All600_Train, type = 'response' )
Final_model4_9_logit_P = predict(model4_9, newdata = All600_Train, type = 'response' )
Final_model4_10_logit_P = predict(model4_10, newdata = All600_Train, type = 'response' )
Final_model4_11_logit_P = predict(model4_11, newdata = All600_Train, type = 'response' )
Final_model4_12_logit_P = predict(model4_12, newdata = All600_Train, type = 'response' )
Final_model4_13_logit_P = predict(model4_13, newdata = All600_Train, type = 'response' )
Final_model4_14_logit_P = predict(model4_14, newdata = All600_Train, type = 'response' )
Final_model4_15_logit_P = predict(model4_15, newdata = All600_Train, type = 'response' )

roc(All600_Train$Dom, Final_model1_1_logit_P)
roc(All600_Train$Dom, Final_model1_2_logit_P)
roc(All600_Train$Dom, Final_model1_3_logit_P)
roc(All600_Train$Dom, Final_model2_1_logit_P)
roc(All600_Train$Dom, Final_model2_2_logit_P)
roc(All600_Train$Dom, Final_model2_3_logit_P)
roc(All600_Train$Dom, Final_model3_1_logit_P)
roc(All600_Train$Dom, Final_model3_2_logit_P)
roc(All600_Train$Dom, Final_model3_3_logit_P)

roc(All600_Train$Dom, Final_model4_1_logit_P)
roc(All600_Train$Dom, Final_model4_2_logit_P)
roc(All600_Train$Dom, Final_model4_3_logit_P)
roc(All600_Train$Dom, Final_model4_4_logit_P)
roc(All600_Train$Dom, Final_model4_5_logit_P)
roc(All600_Train$Dom, Final_model4_6_logit_P)
roc(All600_Train$Dom, Final_model4_7_logit_P)
roc(All600_Train$Dom, Final_model4_8_logit_P)
roc(All600_Train$Dom, Final_model4_9_logit_P)
roc(All600_Train$Dom, Final_model4_10_logit_P)
roc(All600_Train$Dom, Final_model4_11_logit_P)
roc(All600_Train$Dom, Final_model4_12_logit_P)
roc(All600_Train$Dom, Final_model4_13_logit_P)
roc(All600_Train$Dom, Final_model4_14_logit_P)
roc(All600_Train$Dom, Final_model4_15_logit_P)



