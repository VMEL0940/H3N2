## Package ##
library(tidyverse)
library(lme4)
library(pROC)
library(lmtest)
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



tMRCA_Dom_Prop <- train_regdata %>% 
  group_by(Vaccine_code, tMRCA_Dom) %>% 
  summarise(n = n())

tMRCA_Trunk_Prop <- train_regdata %>% 
  group_by(Vaccine_code, tMRCA_Trunk_Dom) %>% 
  summarise(n = n())

##### 1 Control model - RBD and 15A #####
RBD_15A_model <- glmer(tMRCA_Trunk_Dom ~ HA_RBD + HA_15A + (1 | Vaccine_code), 
                       data=train_regdata, family=binomial(link="logit"))

# Summarize the model
summary(RBD_15A_model)

### Only 15A model ###
HA_15A_model <- glmer(tMRCA_Trunk_Dom ~ HA_15A + (1 | Vaccine_code), 
                       data=train_regdata, family=binomial(link="logit"))

# Summarize the model
summary(HA_15A_model)

#### A. Predicted Score ####
HA_15A_logit_P = predict(HA_15A_model, newdata = train_regdata, type = 'response' )
train_regdata$HA_15A_logit_P <- HA_15A_logit_P

#### B. ROC ####
roc_score=roc(train_regdata$tMRCA_Trunk_Dom, HA_15A_logit_P) #AUC score
plot(roc_score, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE,
     print.thres=TRUE, print.thres.pch=19, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="1. ROC curve of the HA15A Dominance Prediction Model")



##### 2 Genetic distance model - HA + #####
## Model selection
HA_model <- glmer(tMRCA_Trunk_Dom ~ HA + (1 | Vaccine_code), 
                    data=train_regdata, family=binomial(link="logit"))    

summary(GD_HA_model)

#### A. Predicted Score ####
HA_logit_P = predict(HA_model, newdata = train_regdata, type = 'response' )
train_regdata$HA_logit_P <- HA_logit_P


#### B. ROC ####
roc_score=roc(train_regdata$tMRCA_Trunk_Dom, HA_logit_P) #AUC score
plot(roc_score, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE,
     print.thres=TRUE, print.thres.pch=19, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="1. ROC curve of the Full Dominance Prediction Model")


## Model selection
ME_model_1 <- glmer(tMRCA_Dom ~ HA + N_A + M  + (1 | Vaccine_code), 
                       data=train_regdata, family=binomial(link="logit"))                     

summary(ME_model_1)
lrtest(ME_model_1, Full_ME_model)


### 2. LR test Final model ###
Fit_ME_model_1 <- glmer(tMRCA_Dom ~ HA + N_A + M  +  (1 | Vaccine_code), 
                       data=train_regdata, family=binomial(link="logit"))

summary(Fit_ME_model_1)


#### A. Predicted Score ####
Fit_logit_P = predict(Fit_ME_model_1, newdata = train_regdata, type = 'response' )
train_regdata$Fit_2_logit_P <- Fit_logit_P

#### B. ROC ####
roc_score=roc(train_regdata$tMRCA_Dom, Fit_logit_P) #AUC score
plot(roc_score, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE,
     print.thres=TRUE, print.thres.pch=19, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="2. ROC curve of the Fitted Dominance Prediction Model")



### 3. All HA model model ###
ALLHA_ME_model_1 <- glmer(tMRCA_Dom ~HA + HA_RBD + HA_15A +  (1 | Vaccine_code), 
                        data=train_regdata, family=binomial(link="logit"))

summary(ALLHA_ME_model_1)

#### A. Predicted Score ####
AHA_logit_P = predict(ALLHA_ME_model_1, newdata = train_regdata, type = 'response' )
train_regdata$AllHA_3_logit_P <- AHA_logit_P

#### B. ROC ####
roc_score=roc(train_regdata$tMRCA_Dom, AHA_logit_P) #AUC score
plot(roc_score, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE,
     print.thres=TRUE, print.thres.pch=19, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="3. ROC curve of the All HA Dominance Prediction Model")


### 4. Only HA model model ###
OnlyHA_ME_model_1 <- glmer(tMRCA_Dom ~ HA + (1 | Vaccine_code), 
                          data=train_regdata, family=binomial(link="logit"))

summary(OnlyHA_ME_model_1)


#### A. Predicted Score ####
OHA_logit_P = predict(OnlyHA_ME_model_1, newdata = train_regdata, type = 'response' )
train_regdata$OnlyHA_4_logit_P <- OHA_logit_P

#### B. ROC ####
roc_score=roc(train_regdata$tMRCA_Dom, OHA_logit_P) #AUC score
plot(roc_score, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE,
     print.thres=TRUE, print.thres.pch=19, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="4. ROC curve of the HA Gen distance Dominance Prediction Model")

### 5. RBD-15A model model ###
RBD15A_ME_model_1 <- glmer(tMRCA_Dom ~ HA_RBD + HA_15A + (1 | Vaccine_code), 
                           data=train_regdata, family=binomial(link="logit"))

summary(RBD15A_ME_model_1)


#### A. Predicted Score ####
RBD15A_logit_P = predict(RBD15A_ME_model_1, newdata = train_regdata, type = 'response' )
train_regdata$RBD15A_5_logit_P <- RBD15A_logit_P

#### B. ROC ####
roc_score=roc(train_regdata$tMRCA_Dom, RBD15A_logit_P) #AUC score
plot(roc_score, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE,
     print.thres=TRUE, print.thres.pch=19, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="5. ROC curve of the RBS-15A Dominance Prediction Model")



### 6. RBD model  ###
RBD_ME_model_1 <- glmer(tMRCA_Dom ~ HA_RBD  + (1 | Vaccine_code), 
                           data=train_regdata, family=binomial(link="logit"))

summary(RBD_ME_model_1)


#### A. Predicted Score ####
RBD_logit_P = predict(RBD_ME_model_1, newdata = train_regdata, type = 'response' )
train_regdata$RBD_6_logit_P <- RBD_logit_P

#### B. ROC ####
roc_score=roc(train_regdata$tMRCA_Dom, RBD_logit_P) #AUC score
plot(roc_score, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE,
     print.thres=TRUE, print.thres.pch=19, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="6. ROC curve of the RBS Dominance Prediction Model")












# Summarize the model
ranef(Full_ME_model)
ranef(ME_model)












