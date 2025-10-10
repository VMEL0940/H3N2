## Package ##
library(tidyverse)
library(lme4)
library(pROC)
library(lmtest)
library(corrplot)

## Set working directory ##
setwd("01_GenDisFlu/")
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


## Making Test data -- HK15
test_regdata <- regdata %>% 
  filter(vaccineStrain == "EPI686117_A_Hong_Kong_15611_2015_NA_NA_NA")

## set the same intercept as Swtz13
test_regdata$Vaccine_code <- "7.Vic11"




##### 1.HA only model #####
GD_1_HA_model <- glmer(subtMRCA_1stDesc ~ HA + (1 | Vaccine_code), 
                       data=train_regdata, family=binomial(link="logit"))    

summary(GD_1_HA_model)

#### A. Predicted Score ####
HA_logit_P_HK15 = predict(GD_1_HA_model, newdata = test_regdata, type = 'response' )
test_regdata$HA_logit_P <- HA_logit_P_HK15


#### 2.3 HA + NA Model #####
GD_2_3_HA_NA_model <- glmer(subtMRCA_1stDesc ~ HA+ N_A + (1 | Vaccine_code), 
                            data=train_regdata, family=binomial(link="logit"))    ### Sig

summary(GD_2_3_HA_NA_model)

#### A. Predicted Score ####
HA_NA_logit_P_HK15 = predict(GD_2_3_HA_NA_model, newdata = test_regdata, type = 'response' )
test_regdata$HA_NA_logit_P <- HA_NA_logit_P_HK15


#### 2.3.1 HA*PA + NA Model #####
GD_2_3_1_HAPA_NA_model <- glmer(subtMRCA_1stDesc ~ HA*PA + N_A + (1 | Vaccine_code), 
                                data=train_regdata, family=binomial(link="logit"))    ### Sig

summary(GD_2_3_1_HAPA_NA_model)

#### A. Predicted Score ####
HAPA_NA_logit_P_HK15 = predict(GD_2_3_1_HAPA_NA_model, newdata = test_regdata, type = 'response' )
test_regdata$HAPA_NA_logit_P <- HAPA_NA_logit_P_HK15


#### 2.3.5 HA*PA + NA*NS Model #####
GD_2_3_5_HAPA_NANS_model <- glmer(subtMRCA_1stDesc ~ HA*PA + N_A*NS + (1 | Vaccine_code), 
                                data=train_regdata, family=binomial(link="logit"))    ### Sig

summary(GD_2_3_5_HAPA_NANS_model)

#### B. Predicted Score -- Test data ####
HAPA_NANS_logit_P_HK15 = predict(GD_2_3_5_HAPA_NANS_model, newdata = test_regdata, type = 'response' )
test_regdata$HAPA_NANS_logit_P <- HAPA_NANS_logit_P_HK15


#### 2.3.6 HA + NA*NS Model #####
GD_2_3_6_HA_NANS_model <- glmer(subtMRCA_1stDesc ~ HA + N_A*NS + (1 | Vaccine_code), 
                                  data=train_regdata, family=binomial(link="logit"))    ### Sig

summary(GD_2_3_6_HA_NANS_model)

#### B. Predicted Score -- Test data ####
HAPA_NANS_logit_P_HK15 = predict(GD_2_3_6_HA_NANS_model, newdata = test_regdata, type = 'response' )
test_regdata$HAPA_NANS_logit_P <- HAPA_NANS_logit_P_HK15





## Data import --- 702 Data ##
All_HA_702 <- read.csv("Data/Veri_701_2023/H3N2_Index_20231218_701.csv", header = T, na.strings = "")

HK15_702 <- All_HA_702 %>% 
  filter(vaccine_code == "HK15")

colnames(HK15_702)[1] <- "ID"

HK15Pred <- left_join(test_regdata, HK15_702, by = "ID")



###############
### B. ROC ####
###############

### 1. HA only ###

roc_score_1_HA <- roc(HK15Pred$Dom, HK15Pred$HA_logit_P) #AUC score

jpeg(filename = "/01_GenDisFlu/Fig/PredModel/1_HA_Model.jpeg", width = 170, height = 150, units = "mm", res = 520)

plot(roc_score_1_HA, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="1. ROC curve of the HA GD Model")

dev.off()


### 2. HA + NA only ###
roc_score_2_HA_NA <- roc(HK15Pred$Dom, HK15Pred$HA_NA_logit_P) #AUC score

jpeg(filename = "/01_GenDisFlu/Fig/PredModel/2_HA_NA_Model.jpeg", width = 170, height = 150, units = "mm", res = 520)

plot(roc_score_2_HA_NA, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="2. ROC curve of the HA + NA GD Model")

dev.off()



#### B. ROC ####
roc_score_3_HAPA_NA <- roc(HK15Pred$Dom, HK15Pred$HAPA_NA_logit_P) #AUC score

jpeg(filename = "/01_GenDisFlu/Fig/PredModel/3_HAPA_NA_Model.jpeg", width = 170, height = 150, units = "mm", res = 520)

plot(roc_score_3_HAPA_NA, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="3. ROC curve of the HA + NA + PA GD Model")

dev.off()


roc_score_4_HAPA_NANS <- roc(HK15Pred$Dom, HK15Pred$HAPA_NANS_logit_P) #AUC score

jpeg(filename = "/01_GenDisFlu/Fig/PredModel/4_HAPA_NANS_Model.jpeg", width = 170, height = 150, units = "mm", res = 520)

plot(roc_score_4_HAPA_NANS, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="4. ROC curve of the HA + NA + PA + NS GD Model - HK15")

dev.off()


HK15Pred$Dom

## Box plot ##
HK15Pred %>% 
  ggplot() +
  geom_jitter(aes(x= as.factor(Dom), y= HA_logit_P), width = 0.25, size = 0.50, color = "gray40") +
  geom_boxplot(aes(x= as.factor(Dom), y = HA_logit_P, color = as.factor(Dom)), size = 0.8, outlier.shape = NA, alpha = 0) +
  ylab("Prediction score 1")+
  xlab("Dominance")+
  theme_bw()


HK15Pred %>% 
  ggplot() +
  geom_jitter(aes(x= as.factor(Dom), y= HA_NA_logit_P), width = 0.25, size = 0.50, color = "gray40") +
  geom_boxplot(aes(x= as.factor(Dom), y = HA_NA_logit_P, color = as.factor(Dom)), size = 0.8, outlier.shape = NA, alpha = 0) +
  ylab("Prediction score 2")+
  xlab("Dominance")+
  theme_bw()


HK15Pred %>% 
  ggplot() +
  geom_jitter(aes(x= as.factor(Dom), y= HAPA_NA_logit_P), width = 0.25, size = 0.50, color = "gray40") +
  geom_boxplot(aes(x= as.factor(Dom), y = HAPA_NA_logit_P, color = as.factor(Dom)), size = 0.8, outlier.shape = NA, alpha = 0) +
  ylab("Prediction score 3")+
  xlab("Dominance")+
  theme_bw()


HK15Pred %>% 
  ggplot() +
  geom_jitter(aes(x= as.factor(Dom), y= HAPA_NANS_logit_P), width = 0.25, size = 0.50, color = "gray40") +
  geom_boxplot(aes(x= as.factor(Dom), y = HAPA_NANS_logit_P, color = as.factor(Dom)), size = 0.8, outlier.shape = NA, alpha = 0) +
  ylab("Prediction score 4 ")+
  xlab("Dominance")+
  theme_bw()

