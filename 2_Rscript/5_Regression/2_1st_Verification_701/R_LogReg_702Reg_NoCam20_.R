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

### Subset data - Exclude Dar21 ##
train_regdata <- regdata %>% 
  filter(vaccine_code != "Dar21" & vaccine_code != "Cam20")


## Model selection ## 
##### 1.HA only model #####
GD_1_HA_model <- glmer(Dom ~ HA + (1 | vaccine_code), 
                       data=train_regdata, family=binomial(link="logit"))    

summary(GD_1_HA_model)

#### A. Predicted Score ####
HA_logit_P = predict(GD_1_HA_model, newdata = train_regdata, type = 'response' )
train_regdata$HA_logit_P <- HA_logit_P

#### B. ROC ####
roc_score_1_HA <- roc(train_regdata$Dom, HA_logit_P) #AUC score

jpeg(filename = "/01_GenDisFlu/Fig/PredModel/1_HA_Model.jpeg", width = 170, height = 150, units = "mm", res = 520)

plot(roc_score_1_HA, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="1. ROC curve of the HA GD Model")

dev.off()


##### 2.HA + 1 Model #####
#### 2.1 HA + PA Model #####
GD_2_HA_PA_model <- glmer(Dom ~ HA + PA + (1 | vaccine_code), 
                          data=train_regdata, family=binomial(link="logit"))    ## No sig

summary(GD_2_HA_PA_model)

lrtest(GD_1_HA_model, GD_2_HA_PA_model)


#### 2.1.1 HA X PA Model #####
GD_2_2_1_HAPA_model <- glmer(Dom ~ HA*PA + (1 | vaccine_code), 
                           data=train_regdata, family=binomial(link="logit"))    ### With interaction 

summary(GD_2_2_1_HAPA_model)

lrtest(GD_1_HA_model, GD_2_2_1_HAPA_model)



#### 2.2. HA + NA Model #####
GD_2_2_HA_NA_model <- glmer(Dom ~ HA + N_A + (1 | vaccine_code), 
                            data=train_regdata, family=binomial(link="logit"))    ### Sig

summary(GD_2_2_HA_NA_model)

lrtest(GD_1_HA_model, GD_2_2_HA_NA_model)

#### 2.2.1 HA * NA Model ##### ---- Significant
GD_2_2_1_HANA_model <- glmer(Dom ~ HA*N_A + (1 | vaccine_code), 
                           data=train_regdata, family=binomial(link="logit"))    ### No interaction!

summary(GD_2_2_1_HANA_model)

lrtest(GD_1_HA_model, GD_2_2_1_HANA_model)


#### A. Predicted Score ####
HANA_logit_P = predict(GD_2_2_1_HANA_model, newdata = train_regdata, type = 'response' )
train_regdata$HANA_logit_P <- HANA_logit_P

#### B. ROC ####
roc_score_2_HANA <- roc(train_regdata$Dom, HANA_logit_P) #AUC score

jpeg(filename = "/01_GenDisFlu/Fig/PredModel/702_2_HAPA_NA_Model.jpeg", width = 170, height = 150, units = "mm", res = 520)

plot(roc_score_2_HANA, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="2. ROC curve of the HA X NA GD Model")

dev.off()


#### 2.3 HA + NS Model #####
GD_2_3_HA_NS_model <- glmer(Dom ~ HA + NS + (1 | vaccine_code), 
                            data=train_regdata, family=binomial(link="logit"))    ### No interaction!

summary(GD_2_3_HA_NS_model)

lrtest(GD_1_HA_model, GD_2_3_HA_NS_model)

#### 2.3.1 HA * NS Model #####
GD_2_3_1_HANS_model <- glmer(Dom ~ HA * NS + (1 | vaccine_code), 
                           data=train_regdata, family=binomial(link="logit"))    ### No interaction!

summary(GD_2_3_1_HANS_model)

lrtest(GD_1_HA_model, GD_2_3_1_HANS_model)


#### 2.4 HA + PB2 Model #####
GD_2_4_HA_PB2_model <- glmer(Dom ~ HA + PB2 + (1 | vaccine_code), 
                           data=train_regdata, family=binomial(link="logit"))    ### No interaction!

summary(GD_2_4_HA_PB2_model)

lrtest(GD_1_HA_model, GD_2_4_HA_PB2_model)


#### 2.4.1 HA * PB2 Model #####
GD_2_4_1_HAPB2_model <- glmer(Dom ~ HA * PB2 + (1 | vaccine_code), 
                             data=train_regdata, family=binomial(link="logit"))    ### No interaction!

summary(GD_2_4_1_HAPB2_model)

lrtest(GD_1_HA_model, GD_2_4_1_HAPB2_model)


#### 2.5 HA + PB1 Model #####
GD_2_5_HA_PB1_model <- glmer(Dom ~ HA + PB1 + (1 | vaccine_code), 
                              data=train_regdata, family=binomial(link="logit"))    ### No interaction!

summary(GD_2_5_HA_PB1_model)

lrtest(GD_1_HA_model, GD_2_5_HA_PB1_model)


#### 2.5.1 HA * PB1 Model #####
GD_2_5_1_HAPB1_model <- glmer(Dom ~ HA * PB1 + (1 | vaccine_code), 
                             data=train_regdata, family=binomial(link="logit"))    ### No interaction!

summary(GD_2_5_1_HAPB1_model)

lrtest(GD_1_HA_model, GD_2_5_1_HAPB1_model)



#### 2.6 HA + NP Model #####
GD_2_6_HA_NP_model <- glmer(Dom ~ HA + NP + (1 | vaccine_code), 
                              data=train_regdata, family=binomial(link="logit"))    ### No interaction!

summary(GD_2_6_HA_NP_model)

lrtest(GD_1_HA_model, GD_2_6_HA_NP_model)



#### 2.6.1 HA * NP Model #####
GD_2_6_1_HANP_model <- glmer(Dom ~ HA * NP + (1 | vaccine_code), 
                            data=train_regdata, family=binomial(link="logit"))    ### No interaction!

summary(GD_2_6_1_HANP_model)

lrtest(GD_1_HA_model, GD_2_6_1_HANP_model)



#### 2.7 HA + M Model #####
GD_2_7_HA_M_model <- glmer(Dom ~ HA + M + (1 | vaccine_code), 
                            data=train_regdata, family=binomial(link="logit"))    ### No interaction!

summary(GD_2_7_HA_M_model)

lrtest(GD_1_HA_model, GD_2_7_HA_M_model)


#### 2.7.1 HA X M Model ##### -- Signficant but unstable
GD_2_7_1_HA_M_model <- glmer(Dom ~ HA * M + (1 | vaccine_code), 
                           data=train_regdata, family=binomial(link="logit"))    ### No interaction!

summary(GD_2_7_1_HA_M_model)

lrtest(GD_1_HA_model, GD_2_7_1_HA_M_model)

#### A. Predicted Score ####
HAM_logit_P = predict(GD_2_7_1_HA_M_model, newdata = train_regdata, type = 'response' )
train_regdata$HAM_logit_P <- HAM_logit_P

#### B. ROC ####
roc_score_3_HAM <- roc(train_regdata$Dom, HAM_logit_P) #AUC score

jpeg(filename = "/01_GenDisFlu/Fig/PredModel/702_2_HAM_Model.jpeg", width = 170, height = 150, units = "mm", res = 520)

plot(roc_score_3_HAM, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="3. ROC curve of the HA X M GD Model")

dev.off()

