## Package ##
library(tidyverse)
library(lme4)
library(pROC)
library(lmtest)
library(corrplot)

####### Preset for regression data ###########

## Set working directory ##
setwd("01_GenDisFlu/")
#setwd("01_GenDisFlu/") ## for mac

## Data import ##
AllG_NS <- read.csv("Data/Genetic Distances/Verification_702/Summary/AllG_NS.csv", header = T, na.strings = "")

## Remove VAccine strain ##
AllG_NS <- AllG_NS %>% 
  filter(compareStrain != vaccineStrain)

## Cleaning compare strain name ## 
AllG_NS$compareStrain_V1 <- gsub("_NA_NA_NA_NA", "", AllG_NS$compareStrain)
AllG_NS$compareStrain_V1 <- gsub("_NA_NA_NA", "", AllG_NS$compareStrain_V1)
AllG_NS$compareStrain_V1 <- gsub("_NA_NA", "", AllG_NS$compareStrain_V1)
AllG_NS$compareStrain <- gsub("2009_NA", "2009", AllG_NS$compareStrain_V1)

## Spread data -- making ind. Variable ## 
Ind <- AllG_NS[,c(3,2,4,6)] %>% 
  spread(key = Gene, value = Nmedian)

## Rename - compareStrain to ID ## 
colnames(Ind)[1] <- "ID"

## Load data of 702 - Dominanat ## 
Dom <- read.csv("Data/Veri_701_2023/H3N2_Index_20231218_701.csv", header = T, na.strings = "")

## Rename - compareStrain to ID ## 
colnames(Dom)[1] <- "ID"

## leftjoin * subset ##
regdata <- right_join(Ind, Dom, by = "ID")

## remove  Moscow Strain ##
regdata <- regdata  %>% 
  filter(ID != "EPI103320_A_Moscow_10_1999")

## Rename NA to N_A ##
colnames(regdata)[5] <- "N_A"

## remove vaccineStrain column ##
regdata <- regdata[c(-2)]

write.csv(regdata, "regdata.csv")

####### Preset for Prediction of 12. Cam20 ###########

### Subset data - Training ##
train_regdata <- regdata %>% 
  filter(vaccine_code != "Dar21" & vaccine_code != "Cam20" 
         & vaccine_code != "Mos99" & vaccine_code != "Fuj02" & vaccine_code != "Cal04")

unique(train_regdata$vaccine_code)

### Subset data - Testing  ##
test_regdata <- regdata %>% 
  filter(vaccine_code == "Cam20")

#unique(test_regdata$vaccine_code) -- Only 1 group 

test_regdata$vaccine_code <- "HK19"

###### Model selection ########### 

##### 1.HA only model #####
GD_1_HA_model <- glmer(Dom ~ HA + (1 | vaccine_code), 
                       data=train_regdata, family=binomial(link="logit"))    

summary(GD_1_HA_model)
AIC(GD_1_HA_model)
### TRAINING ROC ##

#### A. Predicted Score ####
HA_logit_P = predict(GD_1_HA_model, newdata = train_regdata, type = 'response' )
train_regdata$HA_logit_P <- HA_logit_P

#### B. ROC ####
roc_score_1_HA <- roc(train_regdata$Dom, HA_logit_P) #AUC score

jpeg(filename = "01_GenDisFlu/Fig/PredModel/4_Cam20Model/1_HA_Model.jpeg", width = 170, height = 150, units = "mm", res = 520)

plot(roc_score_1_HA, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="1. ROC curve of the HA GD Model for Training Data")

dev.off()

## Cross-validation using 12. Cam20 ##
## Test data ##

#### A. Predicted Score ####
HA_logit_P_CV = predict(GD_1_HA_model, newdata = test_regdata, type = 'response' )
test_regdata$HA_logit_P_CV <- HA_logit_P_CV

#### B. ROC ####
roc_score_1_HA_CV <- roc(test_regdata$Dom, HA_logit_P_CV) #AUC score

jpeg(filename = "01_GenDisFlu/Fig/PredModel/4_Cam20Model/1_HA_Model_CV.jpeg", width = 170, height = 150, units = "mm", res = 520)

plot(roc_score_1_HA_CV, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="1.1. ROC curve of the HA GD Model for Test Data")

dev.off()




##### 2.HA + 1 Model #####

#### 2.1 HA + PA Model #####
GD_2_HA_PA_model <- glmer(Dom ~ HA + PA + (1 | vaccine_code), 
                          data=train_regdata, family=binomial(link="logit"))    ## No sig

summary(GD_2_HA_PA_model)

lrtest(GD_1_HA_model, GD_2_HA_PA_model)  ## Significant -- Yes 

#### 2.1.1 HA X PA Model #####
GD_2_2_1_HAPA_model <- glmer(Dom ~ HA*PA + (1 | vaccine_code), 
                             data=train_regdata, family=binomial(link="logit"))    ### With interaction 

summary(GD_2_2_1_HAPA_model)

lrtest(GD_1_HA_model, GD_2_2_1_HAPA_model) ## Significant but not good predictors


#### 2.2. HA + NA Model #####
GD_2_2_HA_NA_model <- glmer(Dom ~ HA + N_A + (1 | vaccine_code), 
                            data=train_regdata, family=binomial(link="logit"))    

summary(GD_2_2_HA_NA_model)
lrtest(GD_1_HA_model, GD_2_2_HA_NA_model) ## Non-Significant


#### 2.2.1 HA * NA Model ##### ---- Significant
GD_2_2_1_HANA_model <- glmer(Dom ~ HA*N_A + (1 | vaccine_code), 
                             data=train_regdata, family=binomial(link="logit"))    ### No interaction!

summary(GD_2_2_1_HANA_model)

lrtest(GD_1_HA_model, GD_2_2_1_HANA_model) ## Significant -- Yes!! 


#### 2.3 HA + NS Model #####
GD_2_3_HA_NS_model <- glmer(Dom ~ HA + NS + (1 | vaccine_code), 
                            data=train_regdata, family=binomial(link="logit"))    ### No interaction!

summary(GD_2_3_HA_NS_model)

lrtest(GD_1_HA_model, GD_2_3_HA_NS_model) ## Significant


#### 2.3.1 HA * NS Model #####
GD_2_3_1_HANS_model <- glmer(Dom ~ HA * NS + (1 | vaccine_code), 
                             data=train_regdata, family=binomial(link="logit"))    ### No interaction!

summary(GD_2_3_1_HANS_model)

lrtest(GD_1_HA_model, GD_2_3_1_HANS_model) ## Significant but not good predictors


#### 2.2.1. HA*NA Model #####
GD_2_2_1_HANA_model <- glmer(Dom ~ HA*N_A + (1 | vaccine_code), 
                             data=train_regdata, family=binomial(link="logit"))    ### No interaction!
AIC(GD_2_2_1_HANA_model)
### TRAINING ROC ##
#### A. Predicted Score ####
HANA_logit_P = predict(GD_2_2_1_HANA_model, newdata = train_regdata, type = 'response' )
train_regdata$HANA_logit_P <- HANA_logit_P

#### B. ROC ####
roc_score_2_2_1_HANA <- roc(train_regdata$Dom, HANA_logit_P) #AUC score

jpeg(filename = "01_GenDisFlu/Fig/PredModel/4_Cam20Model/2_2_HANA_Model.jpeg", width = 170, height = 150, units = "mm", res = 520)

plot(roc_score_2_2_1_HANA, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="2. ROC curve of the HA * NA GD Model for Training Data")

dev.off()

## Cross-validation using 12. Cam20 ##
## Test data ##

#### A. Predicted Score ####
HANA_logit_P_CV = predict(GD_2_2_1_HANA_model, newdata = test_regdata, type = 'response' )
test_regdata$HANA_logit_P_CV <- HANA_logit_P_CV

#### B. ROC ####
roc_score_2_2_HANA_CV <- roc(test_regdata$Dom, HANA_logit_P_CV) #AUC score

jpeg(filename = "01_GenDisFlu/Fig/PredModel/4_Cam20Model/2_2_HANA_Model_CV.jpeg", width = 170, height = 150, units = "mm", res = 520)

plot(roc_score_2_2_HANA_CV, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="2.1. ROC curve of the HA * NA GD Model for Test Data")

dev.off()



##### 3.HA * NA + 1 Model #####

#### 3.1 HA * NA + PA Model #####
GD_3_1_HANA_PA_model <- glmer(Dom ~ HA * N_A + PA + (1 | vaccine_code), 
                              data=train_regdata, family=binomial(link="logit"))    ### No interaction!

summary(GD_3_1_HANA_PA_model)

lrtest(GD_2_2_1_HANA_model, GD_3_1_HANA_PA_model) ## Significant -- Yes
AIC(GD_3_1_HANA_PA_model)

#### 3.2 HA * NA + NS Model #####
GD_3_2_HANA_NS_model <- glmer(Dom ~ HA * N_A + NS + (1 | vaccine_code), 
                              data=train_regdata, family=binomial(link="logit"))    ### No interaction!

summary(GD_3_2_HANA_NS_model)

lrtest(GD_2_2_1_HANA_model, GD_3_2_HANA_NS_model) ## Non - Significant -- NO



### TRAINING ROC ##
#### A. Predicted Score ####
HANA_PA_logit_P = predict(GD_3_1_HANA_PA_model, newdata = train_regdata, type = 'response' )
train_regdata$HANA_PA_logit_P <- HANA_PA_logit_P

#### B. ROC ####
roc_score_3_1_HANA_PA <- roc(train_regdata$Dom, HANA_PA_logit_P) #AUC score

jpeg(filename = "01_GenDisFlu/Fig/PredModel/4_Cam20Model/3_1_HANA_PA_Model.jpeg", width = 170, height = 150, units = "mm", res = 520)

plot(roc_score_3_1_HANA_PA, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="3. ROC curve of the HA * NA + PA GD Model for Training Data")

dev.off()

## Cross-validation using 12. Cam20 ##
## Test data ##

#### A. Predicted Score ####
HANA_PA_logit_P_CV = predict(GD_3_1_HANA_PA_model, newdata = test_regdata, type = 'response' )
test_regdata$HANA_PA_logit_P_CV <- HANA_PA_logit_P_CV

#### B. ROC ####
roc_score_3_1_HANA_PA_CV <- roc(test_regdata$Dom, HANA_PA_logit_P_CV) #AUC score

jpeg(filename = "01_GenDisFlu/Fig/PredModel/4_Cam20Model/3_1_HANA_PA_Model_CV.jpeg", width = 170, height = 150, units = "mm", res = 520)

plot(roc_score_3_1_HANA_PA_CV, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="3.1 ROC curve of the HA * NA + PA GD Model for Test Data")

dev.off()



##### 4.HA * NA + PA + 1 Model #####
#### 4.1 HA * NA + PA + NS Model #####
GD_4_1_HANA_PA_NS_model <- glmer(Dom ~ HA * N_A + PA + NS + (1 | vaccine_code), 
                               data=train_regdata, family=binomial(link="logit"))    ### No interaction!

summary(GD_4_1_HANA_PA_NS_model)

lrtest(GD_3_1_HANA_PA_model, GD_4_1_HANA_PA_NS_model) ## Non - Significant -- NO


#### 4.2 HA * NA + PA*NS Model #####
GD_4_2_HANA_PANS_model <- glmer(Dom ~ HA*N_A + PA*NS + (1 | vaccine_code), 
                                data=train_regdata, family=binomial(link="logit"))    ### No interaction!

summary(GD_4_2_HANA_PANS_model)

lrtest(GD_3_1_HANA_PA_model, GD_4_2_HANA_PANS_model) ## Non - Significant -- NO


