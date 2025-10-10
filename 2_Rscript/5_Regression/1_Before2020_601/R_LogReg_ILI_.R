## Package ##
library(tidyverse)
library(lme4)
library(pROC)
## Set working directory ##
setwd("/01_GenDisFlu/")

## Data import ##
AllG_NS <- read.csv("Data/Genetic Distances/GeneCompete/AllG_NS.csv", header = T, na.strings = "")

AllG_NS$compareStrain_V1 <- gsub("_NA_NA_NA_NA", "", AllG_NS$compareStrain)
AllG_NS$compareStrain_V1 <- gsub("_NA_NA_NA", "", AllG_NS$compareStrain_V1)
AllG_NS$compareStrain_V1 <- gsub("_NA_NA", "", AllG_NS$compareStrain_V1)
AllG_NS$compareStrain <- gsub("2009_NA", "2009", AllG_NS$compareStrain_V1)

AllG_NS$Vaccine_code <- factor(AllG_NS$vaccineStrain, 
                                levels = c("EPI103320_A_Moscow_10_1999_NA_NA_NA_NA",
                                           "EPI358781_A_Fujian_411_2002_NA_NA_NA_NA", 
                                           "EPI367109_A_California_7_2004_NA_NA_NA_NA",
                                           "EPI502253_A_Wisconsin_67_2005_NA_NA_NA_NA",
                                           "EPI577980_A_Brisbane_10_2007_NA_NA_NA_NA",
                                           "EPI577969_A_Perth_16_2009_NA_NA_NA_NA", 
                                           "EPI417234_A_Victoria_361_2011_NA_NA_NA_NA",
                                           "EPI614441_A_Switzerland_9715293_2013_NA_NA_NA_NA",
                                           "EPI686117_A_Hong_Kong_15611_2015_NA_NA_NA"),
                                labels = c("1.Mos99",
                                           "2.Fuj02", 
                                           "3.Cal04",
                                           "4.Wis05",
                                           "5.Bris07",
                                           "6.Prth09", 
                                           "7.Vic11",
                                           "8.Swtz13",
                                           "9.HK15"))


Ind <- AllG_NS[,c(3,12,4,6)] %>% 
  spread(key = Gene, value = Nmedian)

colnames(Ind)[1] <- "ID"

## leftjoin * subset ##
Dom <- read.csv("Data/Dom.csv", header = T, na.strings = "")

regdata <- right_join(Ind, Dom, by = "ID")

colnames(regdata)[5] <- "N_A"

Dom2_Prop <- regdata %>% 
  group_by(Vaccine_code, Dom2) %>% 
  summarise(n = n())

DomTrunk_Prop <- regdata %>% 
  group_by(Vaccine_code, DomTrunk) %>% 
  summarise(n = n())


### Subset data - Exclude HK15 ##
train_regdata <- regdata %>% 
  filter(Vaccine_code != "HK15")


### 1. Full model ###

## multilevel multivaraible logistic regression ##
Full_ME_model <- glmer(Dominance ~ PB2+ PB1 + PA + HA + NP + N_A + M + NS + 
                         Reasrt + HA_RBD + HA_15A + CoEvo_PB2 + (1 | Vaccine_code), 
                       data=train_regdata, family=binomial(link="logit"))

# Summarize the model
summary(Full_ME_model)


#### A. Predicted Score ####
Full_logit_P = predict(Full_ME_model, newdata = train_regdata, type = 'response' )
train_regdata$Full_1_logit_P <- Full_logit_P



#### B. Predicted Score ####
train_regdata$Full_1_logit_P <- Full_logit_P



### 2. LM test Final model ###



### 3. Strict Final model ###



### 4. All HA model model ###



### 5. Only HA model model ###




### 6. RBD-15A model model ###



### 7. RBD model model ###




### Final model ##
## multilevel multivaraible logistic regression ##
ME_model <- glmer(Dominance ~ PA + HA + N_A + M + NS +
                    HA_RBD + HA_15A + (1 | Vaccine_code), 
                       data=train_regdata, family=binomial(link="logit"))


# Summarize the model
summary(ME_model)

Rst_logit_P = predict(ME_model, newdata = train_regdata, type = 'response' )
train_regdata$Rst_logit_P <- Rst_logit_P

ranef(Full_ME_model)
ranef(ME_model)



### 1. RBD Only model  ##
## multilevel multivaraible logistic regression ##
RBD_ME_model <- glmer(Dominance ~ HA_RBD + (1 | Vaccine_code), 
                       data=train_regdata, family=binomial(link="logit"))

RBD15A_ME_model <- glmer(Dominance ~ HA_RBD + HA_15A + (1 | Vaccine_code), 
                       data=train_regdata, family=binomial(link="logit"))

HA_ME_model <- glmer(Dominance ~ HA + (1 | Vaccine_code), 
                           data=train_regdata, family=binomial(link="logit"))




# Summarize the model
summary(RBD_ME_model)
### Final model ##
summary(RBD15A_ME_model)
### Final model ##
summary(HA_ME_model)



RBD_logit_P = predict(RBD_ME_model, newdata = train_regdata, type = 'response' )
RBD15A_logit_P = predict(RBD15A_ME_model, newdata = train_regdata, type = 'response' )
HA_logit_P = predict(HA_ME_model, newdata = train_regdata, type = 'response' )

train_regdata$RBD_logit_P <- RBD_logit_P
train_regdata$RBD15A_logit_P <- RBD15A_logit_P
train_regdata$HA_logit_P <- HA_logit_P









