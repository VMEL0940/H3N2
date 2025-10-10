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

### 1. Base model ###
Basemodel <- glmer(Dom ~ HA_Nonsyn + (1 | vaccine_code), 
                  data = All600_Train, family = binomial)

### 2. Final model ###
Finalmodel <- glmer(Dom ~ HA_Nonsyn + NA_Nonsyn * NS_Nonsyn  + (1 | vaccine_code), 
                  data = All600_Train, family = binomial)

### 3. Control model 1 ###
Contrl1model <- glmer(Dom ~ HA_RBD + HA_15A + (1 | vaccine_code), 
                   data = All600_Train, family = binomial)

### 4. Control model 2 ###
Contrl2model <- glmer(Dom ~ HA_Koel + (1 | vaccine_code), 
                   data = All600_Train, family = binomial)

#### A. Predicted Score ####
Base_model_1_logit_P = predict(Basemodel, newdata = All600_Train, type = 'response' )
Final_model_2_logit_P = predict(Finalmodel, newdata = All600_Train, type = 'response' )
Control1_model_3_logit_P = predict(Contrl1model, newdata = All600_Train, type = 'response' )
Control2_model_4_logit_P = predict(Contrl2model, newdata = All600_Train, type = 'response' )

Base600 <- roc(All600_Train$Dom, Base_model_1_logit_P)
Final600 <- roc(All600_Train$Dom, Final_model_2_logit_P)
Control1_600 <- roc(All600_Train$Dom, Control1_model_3_logit_P)
Control2_600 <- roc(All600_Train$Dom, Control2_model_4_logit_P)

jpeg(filename = "01_GenDisFlu/Fig/PredModel/0_601Model/1_Base_Model.jpeg", width = 170, height = 150, units = "mm", res = 400)

plot(Base600, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="The Base model with train data")

dev.off()

jpeg(filename = "01_GenDisFlu/Fig/PredModel/0_601Model/2_Final_Model.jpeg", width = 170, height = 150, units = "mm", res = 400)

plot(Final600, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="The Final model with train data")

dev.off()


jpeg(filename = "01_GenDisFlu/Fig/PredModel/0_601Model/3_Ctrl1_Model.jpeg", width = 170, height = 150, units = "mm", res = 400)

plot(Control1_600, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="The Control model 1 with train data")

dev.off()

jpeg(filename = "01_GenDisFlu/Fig/PredModel/0_601Model/4_Ctrl2_Model.jpeg", width = 170, height = 150, units = "mm", res = 400)

plot(Control2_600, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="The Control model 2 with train data")

dev.off()



