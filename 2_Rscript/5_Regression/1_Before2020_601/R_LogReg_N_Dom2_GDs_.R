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

##### 1 Base model - RBD-15A #####
RBD_15A_model <- glmer(subtMRCA_1stDesc ~ HA_RBD + HA_15A + (1 | Vaccine_code), 
                    data=train_regdata, family=binomial(link="logit"))

# Summarize the model
summary(RBD_15A_model)

#### A. Predicted Score ####
RBD_15A_logit_P = predict(RBD_15A_model, newdata = train_regdata, type = 'response' )
train_regdata$RBD_15A_logit_P <- RBD_15A_logit_P

#### B. ROC ####
roc_score=roc(train_regdata$subtMRCA_1stDesc, RBD_15A_logit_P) #AUC score

jpeg(filename = "/01_GenDisFlu/Fig/PredModel/1_RBS15A_Model.jpeg", width = 170, height = 150, units = "mm", res = 520)

plot(roc_score, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE,
     print.thres=TRUE, print.thres.pch=19, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="B1. ROC curve of the HA, RBD-15A position Model")

dev.off()


##### 2 Base model - Koel #####
Koel_model <- glmer(subtMRCA_1stDesc ~ HA_Koel + (1 | Vaccine_code), 
                       data=train_regdata, family=binomial(link="logit"))

# Summarize the model
summary(Koel_model)

#### A. Predicted Score ####
Koel_logit_P = predict(Koel_model, newdata = train_regdata, type = 'response' )
train_regdata$Koel_logit_P <- Koel_logit_P

#### B. ROC ####
roc_score=roc(train_regdata$subtMRCA_1stDesc, Koel_logit_P) #AUC score

jpeg(filename = "/01_GenDisFlu/Fig/PredModel/2_Koel_Model.jpeg", width = 170, height = 150, units = "mm", res = 520)

plot(roc_score, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE,
     print.thres=TRUE, print.thres.pch=19, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="B2. ROC curve of the HA, Koel's position Model")

dev.off()


############# Model selection #############
##### HA only model #####

GD_1_HA_model <- glmer(subtMRCA_1stDesc ~ HA + (1 | Vaccine_code), 
                       data=train_regdata, family=binomial(link="logit"))    

summary(GD_1_HA_model)

#### A. Predicted Score ####
HA_logit_P = predict(GD_1_HA_model, newdata = train_regdata, type = 'response' )
train_regdata$HA_logit_P <- HA_logit_P

#### B. ROC ####
roc(train_regdata$subtMRCA_1stDesc, HA_logit_P) #AUC score



GD_2_HA_PA_model <- glmer(subtMRCA_1stDesc ~ HA + PA + (1 | Vaccine_code), 
                       data=train_regdata, family=binomial(link="logit"))    ## No sig

summary(GD_2_HA_PA_model)

#### A. Predicted Score ####
HA_PA_logit_P = predict(GD_2_HA_PA_model, newdata = train_regdata, type = 'response' )
train_regdata$HA_PA_logit_P <- HA_PA_logit_P

#### B. ROC ####
roc(train_regdata$subtMRCA_1stDesc, HA_PA_logit_P) #AUC score


GD_2_1_HAPA_model <- glmer(subtMRCA_1stDesc ~ HA*PA + (1 | Vaccine_code), 
                          data=train_regdata, family=binomial(link="logit"))    ### With interaction 

summary(GD_2_1_HAPA_model)


#### A. Predicted Score ####
HA_PA_logit_P = predict(GD_2_HA_PA_model, newdata = train_regdata, type = 'response' )
train_regdata$HA_PA_logit_P <- HA_PA_logit_P

#### B. ROC ####
roc(train_regdata$subtMRCA_1stDesc, HA_PA_logit_P) #AUC score



lrtest(GD_1_HA_model, GD_2_1_HAPA_model)

GD_2_2_HA_NA_model <- glmer(subtMRCA_1stDesc ~ HA+ N_A + (1 | Vaccine_code), 
                           data=train_regdata, family=binomial(link="logit"))    ### Sig

summary(GD_2_2_HA_NA_model)

lrtest(GD_1_HA_model, GD_2_2_HA_NA_model)

GD_2_3_HANA_model <- glmer(subtMRCA_1stDesc ~ HA*N_A + (1 | Vaccine_code), 
                            data=train_regdata, family=binomial(link="logit"))    ### No interaction!

summary(GD_2_3_HANA_model)

lrtest(GD_1_HA_model, GD_2_3_HANA_model)

GD_2_4_HA_NS_model <- glmer(subtMRCA_1stDesc ~ HA + NS + (1 | Vaccine_code), 
                           data=train_regdata, family=binomial(link="logit"))    ### NO sig

summary(GD_2_4_HA_NS_model)

lrtest(GD_1_HA_model, GD_2_4_HA_NS_model)

GD_2_5_HANS_model <- glmer(subtMRCA_1stDesc ~ HA*NS + (1 | Vaccine_code), 
                            data=train_regdata, family=binomial(link="logit"))    ### With interaction 

summary(GD_2_5_HANS_model)

lrtest(GD_1_HA_model, GD_2_5_HANS_model)





### Negative genes 

## Model selection
GDN_1_1_HA_PB2_model <- glmer(subtMRCA_1stDesc ~ HA + PB2  +  (1 | Vaccine_code), 
                              data=train_regdata, family=binomial(link="logit"))    ### NO sig
summary(GDN_1_1_HA_PB2_model)

GDN_1_2_HA_PB1_model <- glmer(subtMRCA_1stDesc ~ HA + PB1  +  (1 | Vaccine_code), 
                              data=train_regdata, family=binomial(link="logit"))    ### NO sig
summary(GDN_1_2_HA_PB1_model)

GDN_1_3_HA_NP_model <- glmer(subtMRCA_1stDesc ~ HA + NP +  (1 | Vaccine_code), 
                             data=train_regdata, family=binomial(link="logit"))    ### NO sig
summary(GDN_1_3_HA_NP_model)

GDN_1_4_HA_M_model <- glmer(subtMRCA_1stDesc ~ HA + M +  (1 | Vaccine_code), 
                            data=train_regdata, family=binomial(link="logit"))    ### NO sig
summary(GDN_1_4_HA_M_model)

## Model selection
GDN_1_5_HAPB2_model <- glmer(subtMRCA_1stDesc ~ HA * PB2  +  (1 | Vaccine_code), 
                             data=train_regdata, family=binomial(link="logit"))    ### With interaction
summary(GDN_1_5_HAPB2_model)

GDN_1_6_HAPB1_model <- glmer(subtMRCA_1stDesc ~ HA * PB1  +  (1 | Vaccine_code), 
                             data=train_regdata, family=binomial(link="logit"))    ### NO sig
summary(GDN_1_6_HAPB1_model)

GDN_1_7_HANP_model <- glmer(subtMRCA_1stDesc ~ HA * NP +  (1 | Vaccine_code), 
                            data=train_regdata, family=binomial(link="logit"))    ### NO sig
summary(GDN_1_7_HANP_model)

GDN_1_8_HAM_model <- glmer(subtMRCA_1stDesc ~ HA * M +  (1 | Vaccine_code), 
                           data=train_regdata, family=binomial(link="logit"))    ### NO sig
summary(GDN_1_8_HAM_model)



#### 1. HA + NA  -- No Interaction
#### 2. HA * NS 
#### 3. HA * PA 
#### 4. HA * PB2 

### Three variables from 2_2 HA + NA Model 

GD_3_7_HA_NANS_model <- glmer(subtMRCA_1stDesc ~ HA + N_A*NS + (1 | Vaccine_code), 
                               data=train_regdata, family=binomial(link="logit"))    ### Sig! Better
summary(GD_3_7_HA_NANS_model)

lrtest(GD_2_2_HA_NA_model, GD_3_7_HA_NANS_model)

GD_3_7_1_HANS_NA_model <- glmer(subtMRCA_1stDesc ~ HA*NS + N_A + (1 | Vaccine_code), 
                              data=train_regdata, family=binomial(link="logit"))    ### Sig! but worse
summary(GD_3_7_1_HANS_NA_model)

lrtest(GD_2_2_HA_NA_model, GD_3_7_1_HANS_NA_model)


GD_3_8_HA_NAPA_model <- glmer(subtMRCA_1stDesc ~ HA + N_A*PA + (1 | Vaccine_code), 
                              data=train_regdata, family=binomial(link="logit"))    ### No sign
summary(GD_3_8_HA_NAPA_model)

GD_3_9_HA_NAPB2_model <- glmer(subtMRCA_1stDesc ~ HA + N_A*PB2 + (1 | Vaccine_code), 
                              data=train_regdata, family=binomial(link="logit"))    ### No sign
summary(GD_3_9_HA_NAPB2_model)



#### 1. HA + NA  -- No Interaction
#### 2. HA  + NA * NS >> HA * NS + NA
#### 3. HA * PA 
#### 4. HA * PB2 


### Three variables from 2_2 HA + NA Model 

GD_3_1_HA_NA_PA_model <- glmer(subtMRCA_1stDesc ~ HA+ N_A + PA+  (1 | Vaccine_code), 
                            data=train_regdata, family=binomial(link="logit"))    ### Better
summary(GD_3_1_HA_NA_PA_model)

lrtest(GD_2_2_HA_NA_model, GD_3_1_HA_NA_PA_model)

GD_3_2_HA_NA_NS_model <- glmer(subtMRCA_1stDesc ~ HA+ N_A + NS +  (1 | Vaccine_code), 
                               data=train_regdata, family=binomial(link="logit"))    ### Better
summary(GD_3_2_HA_NA_NS_model)

lrtest(GD_2_2_HA_NA_model, GD_3_2_HA_NA_NS_model)


GD_3_3_HAPA_NA_model <- glmer(subtMRCA_1stDesc ~ HA*PA + N_A +  (1 | Vaccine_code), 
                               data=train_regdata, family=binomial(link="logit"))    ### Better
summary(GD_3_3_HAPA_NA_model)

lrtest(GD_2_2_HA_NA_model, GD_3_3_HAPA_NA_model)
lrtest(GD_2_1_HAPA_model, GD_3_3_HAPA_NA_model)

GD_3_4_HANS_NA_model <- glmer(subtMRCA_1stDesc ~ HA*NS + N_A  +  (1 | Vaccine_code), 
                               data=train_regdata, family=binomial(link="logit"))    ### Better
summary(GD_3_4_HANS_NA_model)

lrtest(GD_2_2_HA_NA_model, GD_3_4_HANS_NA_model)
lrtest(GD_2_5_HANS_model, GD_3_4_HANS_NA_model)

GD_3_5_HA_NANS_model <- glmer(subtMRCA_1stDesc ~ HA + N_A*NS  +  (1 | Vaccine_code), 
                              data=train_regdata, family=binomial(link="logit"))    ### Better
summary(GD_3_5_HA_NANS_model)

lrtest(GD_2_2_HA_NA_model, GD_3_5_HA_NANS_model)
lrtest(GD_2_5_HANS_model, GD_3_5_HA_NANS_model)

GD_3_6_HA_NAPA_model <- glmer(subtMRCA_1stDesc ~ HA + N_A*PA  +  (1 | Vaccine_code), 
                              data=train_regdata, family=binomial(link="logit"))    ### NO
summary(GD_3_6_HA_NAPA_model)

lrtest(GD_2_2_HA_NA_model, GD_3_6_HA_NAPA_model)

lrtest(GD_2_5_HANS_model, GD_3_5_HA_NANS_model)


### Three variables from 2_2 HA + NA Model 

GD_4_5_HAPA_NANS_model <- glmer(subtMRCA_1stDesc ~ HA*PA + N_A*NS  +  (1 | Vaccine_code), 
                              data=train_regdata, family=binomial(link="logit"))    ### Better
summary(GD_4_5_HA_PA_NANS_model)

lrtest(GD_3_5_HA_NANS_model, GD_4_5_HA_PA_NANS_model)
lrtest(GD_2_5_HANS_model, GD_3_5_HA_NANS_model)

GD_4_6_HAPA_NA_NS_model <- glmer(subtMRCA_1stDesc ~ HA*PA + N_A+ NS  +  (1 | Vaccine_code), 
                                data=train_regdata, family=binomial(link="logit"))    ### Better
summary(GD_4_6_HAPA_NA_NS_model)

lrtest(GD_3_5_HA_NANS_model, GD_4_6_HAPA_NA_NS_model)
lrtest(GD_2_5_HANS_model, GD_3_5_HA_NANS_model)






## Model selection -- PB2 
GDN_2_1_HAPB2_NA_model <- glmer(subtMRCA_1stDesc ~ HA * PB2  + N_A + (1 | Vaccine_code), 
                                data=train_regdata, family=binomial(link="logit"))    ### Sig
summary(GDN_2_1_HAPB2_NA_model)

GDN_2_2_HAPB2_NS_model <- glmer(subtMRCA_1stDesc ~ HA * PB2  + NS + (1 | Vaccine_code), 
                                data=train_regdata, family=binomial(link="logit"))    ### No improvement 
summary(GDN_2_2_HAPB2_NS_model)

GDN_2_3_HAPB2_PA_model <- glmer(subtMRCA_1stDesc ~ HA * PB2  + PA + (1 | Vaccine_code), 
                                data=train_regdata, family=binomial(link="logit"))    ### No improvement 
summary(GDN_2_3_HAPB2_PA_model)


## Model selection -- PB2 
GDN_3_1_HAPB2_NANS_model <- glmer(subtMRCA_1stDesc ~ HA * PB2  + N_A*NS + (1 | Vaccine_code), 
                                data=train_regdata, family=binomial(link="logit"))    ### Not sign
summary(GDN_3_1_HAPB2_NANS_model)



##### 1 Genetic distance model - HA Only #####
## Model selection
GD_1_HA_model <- glmer(subtMRCA_1stDesc ~ HA + (1 | Vaccine_code), 
                    data=train_regdata, family=binomial(link="logit"))    

summary(GD_1_HA_model)

#### A. Predicted Score ####
HA_logit_P = predict(GD_1_HA_model, newdata = train_regdata, type = 'response' )
train_regdata$HA_logit_P <- HA_logit_P

#### B. ROC ####
roc_score=roc(train_regdata$subtMRCA_1stDesc, HA_logit_P) #AUC score

jpeg(filename = "/01_GenDisFlu/Fig/PredModel/3_HA_Model.jpeg", width = 170, height = 150, units = "mm", res = 520)

plot(roc_score, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE,
     print.thres=TRUE, print.thres.pch=19, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="1. ROC curve of the HA Genetic distance Model")

dev.off()






##### 4 Genetic distance model - HA + PA #####

## Model selection
GD_2_HAPA_model <- glmer(subtMRCA_1stDesc ~ HA *PA + (1 | Vaccine_code), 
                          data=train_regdata, family=binomial(link="logit"))    

summary(GD_2_HAPA_model)

#### A. Predicted Score ####
HAPA_logit_P = predict(GD_2_HAPA_model, newdata = train_regdata, type = 'response' )
train_regdata$HAPA_logit_P <- HAPA_logit_P

#### B. ROC ####
roc_score=roc(train_regdata$subtMRCA_1stDesc, HAPA_logit_P) #AUC score

jpeg(filename = "/01_GenDisFlu/Fig/PredModel/4_HAXPA_Model.jpeg", width = 170, height = 150, units = "mm", res = 520)

plot(roc_score, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE,
     print.thres=TRUE, print.thres.pch=19, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="2. ROC curve of the HA & PA Genetic distance Model")

dev.off()




##### 4 Genetic distance model - HA PA + NA #####

## Model selection
GD_3_HAPA_NA_model <- glmer(subtMRCA_1stDesc ~ HA*PA + N_A + (1 | Vaccine_code), 
                       data=train_regdata, family=binomial(link="logit"))    

summary(GD_3_HAPA_NA_model)

#### A. Predicted Score ####
HAPA_NA_logit_P = predict(GD_3_HAPA_NA_model, newdata = train_regdata, type = 'response' )
train_regdata$HAPA_NA_logit_P <- HAPA_NA_logit_P

#### B. ROC ####
roc_score=roc(train_regdata$subtMRCA_1stDesc, HAPA_NA_logit_P) #AUC score

jpeg(filename = "/01_GenDisFlu/Fig/PredModel/5_HAPAXNA_Model.jpeg", width = 170, height = 150, units = "mm", res = 520)

plot(roc_score, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE,
     print.thres=TRUE, print.thres.pch=19, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="3. ROC curve of the HA & PA + NA Genetic distance Model")

dev.off()



##### 5 Genetic distance model - HA PA + NA NS #####

## Model selection
GD_4_HAPA_NANS_model <- glmer(subtMRCA_1stDesc ~ HA*PA + N_A*NS + (1 | Vaccine_code), 
                            data=train_regdata, family=binomial(link="logit"))    

summary(GD_4_HAPA_NANS_model)

#### A. Predicted Score ####
HAPA_NANS_logit_P = predict(GD_4_HAPA_NANS_model, newdata = train_regdata, type = 'response' )
train_regdata$HAPA_NANS_logit_P <- HAPA_NANS_logit_P

#### B. ROC ####
roc_score=roc(train_regdata$subtMRCA_1stDesc, HAPA_NANS_logit_P) #AUC score

jpeg(filename = "/01_GenDisFlu/Fig/PredModel/6_HAPAXNANS_Model.jpeg", width = 170, height = 150, units = "mm", res = 520)

plot(roc_score, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE,
     print.thres=TRUE, print.thres.pch=19, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="4. ROC curve of the HA & PA + NA & NS Genetic distance Model")

dev.off()




##### 5 Genetic distance model - HA PA + NA NS #####

## Model selection
GD_4_HAPA_NANS_model <- glmer(subtMRCA_1stDesc ~ HA*PA + N_A*NS + (1 | Vaccine_code), 
                            data=train_regdata, family=binomial(link="logit"))    

summary(GD_4_HAPA_NANS_model)

#### A. Predicted Score ####
HAPA_NANS_logit_P = predict(GD_4_HAPA_NANS_model, newdata = train_regdata, type = 'response' )
train_regdata$HAPA_NANS_logit_P <- HAPA_NANS_logit_P

#### B. ROC ####
roc_score=roc(train_regdata$subtMRCA_1stDesc, HAPA_NANS_logit_P) #AUC score

jpeg(filename = "/01_GenDisFlu/Fig/PredModel/6_HAPAXNANS_Model.jpeg", width = 170, height = 150, units = "mm", res = 520)

plot(roc_score, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE,
     print.thres=TRUE, print.thres.pch=19, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="4. ROC curve of the HA & PA + NA & NS Genetic distance Model")

dev.off()


##### 6 Genetic distance model - HA + NA NS #####

## Model selection
GD_3_5_HA_NANS_model <- glmer(subtMRCA_1stDesc ~ HA + N_A*NS  +  (1 | Vaccine_code), 
                              data=train_regdata, family=binomial(link="logit"))    ### Better
summary(GD_3_5_HA_NANS_model)


#### A. Predicted Score ####
HA_NANS_logit_P = predict(GD_3_5_HA_NANS_model, newdata = train_regdata, type = 'response' )
train_regdata$HA_NANS_logit_P <- HA_NANS_logit_P

#### B. ROC ####
roc_score=roc(train_regdata$subtMRCA_1stDesc, HA_NANS_logit_P) #AUC score

jpeg(filename = "/01_GenDisFlu/Fig/PredModel/7_HAXNANS_Model.jpeg", width = 170, height = 150, units = "mm", res = 520)

plot(roc_score, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE,
     print.thres=TRUE, print.thres.pch=19, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="4. ROC curve of the HA + NA & NS Genetic distance Model")

dev.off()




