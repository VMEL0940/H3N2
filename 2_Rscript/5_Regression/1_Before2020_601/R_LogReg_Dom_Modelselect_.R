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
     max.auc.polygon=TRUE,print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=19, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="Ctrl1. ROC curve of the HA, RBD-15A Prediction Model")

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
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="Ctrl2. ROC curve of the HA, Koel's position Model")

dev.off()


############# Model selection #############
##### 1.HA only model #####

GD_1_HA_model <- glmer(subtMRCA_1stDesc ~ HA + (1 | Vaccine_code), 
                       data=train_regdata, family=binomial(link="logit"))    

summary(GD_1_HA_model)

#### A. Predicted Score ####
HA_logit_P = predict(GD_1_HA_model, newdata = train_regdata, type = 'response' )
train_regdata$HA_logit_P <- HA_logit_P

#### B. ROC ####
roc_score_1_HA <- roc(train_regdata$subtMRCA_1stDesc, HA_logit_P) #AUC score

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
GD_2_HA_PA_model <- glmer(subtMRCA_1stDesc ~ HA + PA + (1 | Vaccine_code), 
                       data=train_regdata, family=binomial(link="logit"))    ## No sig

summary(GD_2_HA_PA_model)

#### A. Predicted Score ####
HA_PA_logit_P = predict(GD_2_HA_PA_model, newdata = train_regdata, type = 'response' )
train_regdata$HA_PA_logit_P <- HA_PA_logit_P

#### B. ROC ####
roc(train_regdata$subtMRCA_1stDesc, HA_PA_logit_P) #AUC score

lrtest(GD_1_HA_model, GD_2_HA_PA_model)


#### 2.2 HA X PA Model #####
GD_2_2_HAPA_model <- glmer(subtMRCA_1stDesc ~ HA*PA + (1 | Vaccine_code), 
                          data=train_regdata, family=binomial(link="logit"))    ### With interaction 

summary(GD_2_2_HAPA_model)

#### A. Predicted Score ####
HAPA_logit_P = predict(GD_2_2_HAPA_model, newdata = train_regdata, type = 'response' )
train_regdata$HAPA_logit_P <- HAPA_logit_P

#### B. ROC ####
roc(train_regdata$subtMRCA_1stDesc, HAPA_logit_P) #AUC score

lrtest(GD_1_HA_model, GD_2_2_HAPA_model)

#### 2.3 HA + NA Model #####
GD_2_3_HA_NA_model <- glmer(subtMRCA_1stDesc ~ HA+ N_A + (1 | Vaccine_code), 
                           data=train_regdata, family=binomial(link="logit"))    ### Sig

summary(GD_2_3_HA_NA_model)

#### A. Predicted Score ####
HA_NA_logit_P = predict(GD_2_3_HA_NA_model, newdata = train_regdata, type = 'response' )
train_regdata$HA_NA_logit_P <- HA_NA_logit_P

#### B. ROC ####
roc_score_2_HA_NA <- roc(train_regdata$subtMRCA_1stDesc, HA_NA_logit_P) #AUC score

jpeg(filename = "/01_GenDisFlu/Fig/PredModel/2_HA_NA_Model.jpeg", width = 170, height = 150, units = "mm", res = 520)

plot(roc_score_2_HA_NA, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="2. ROC curve of the HA + NA GD Model")

dev.off()




lrtest(GD_1_HA_model, GD_2_3_HA_NA_model)

#### 2.4 HA * NA Model #####
GD_2_4_HANA_model <- glmer(subtMRCA_1stDesc ~ HA*N_A + (1 | Vaccine_code), 
                            data=train_regdata, family=binomial(link="logit"))    ### No interaction!

summary(GD_2_4_HANA_model)

#### A. Predicted Score ####
HANA_logit_P = predict(GD_2_4_HANA_model, newdata = train_regdata, type = 'response' )
train_regdata$HANA_logit_P <- HANA_logit_P

#### B. ROC ####
roc(train_regdata$subtMRCA_1stDesc, HANA_logit_P) #AUC score

lrtest(GD_1_HA_model, GD_2_4_HANA_model)


#### 2.5 HA + NS Model #####
GD_2_5_HA_NS_model <- glmer(subtMRCA_1stDesc ~ HA + NS + (1 | Vaccine_code), 
                           data=train_regdata, family=binomial(link="logit"))    ### No interaction!

summary(GD_2_5_HA_NS_model)

#### A. Predicted Score ####
HA_NS_logit_P = predict(GD_2_5_HA_NS_model, newdata = train_regdata, type = 'response' )
train_regdata$HA_NS_logit_P <- HA_NS_logit_P

#### B. ROC ####
roc(train_regdata$subtMRCA_1stDesc, HA_NS_logit_P) #AUC score

lrtest(GD_1_HA_model, GD_2_5_HA_NS_model)

#### 2.6 HA * NS Model #####
GD_2_6_HANS_model <- glmer(subtMRCA_1stDesc ~ HA * NS + (1 | Vaccine_code), 
                            data=train_regdata, family=binomial(link="logit"))    ### No interaction!

summary(GD_2_6_HANS_model)

#### A. Predicted Score ####
HANS_logit_P = predict(GD_2_6_HANS_model, newdata = train_regdata, type = 'response' )
train_regdata$HANS_logit_P <- HANS_logit_P

#### B. ROC ####
roc(train_regdata$subtMRCA_1stDesc, HANS_logit_P) #AUC score

lrtest(GD_1_HA_model, GD_2_6_HANS_model)


#### 2.3.1 HA*PA + NA Model #####
GD_2_3_1_HAPA_NA_model <- glmer(subtMRCA_1stDesc ~ HA*PA + N_A + (1 | Vaccine_code), 
                            data=train_regdata, family=binomial(link="logit"))    ### Sig

summary(GD_2_3_1_HAPA_NA_model)

#### A. Predicted Score ####
HAPA_NA_logit_P = predict(GD_2_3_1_HAPA_NA_model, newdata = train_regdata, type = 'response' )
train_regdata$HAPA_NA_logit_P <- HAPA_NA_logit_P

#### B. ROC ####
roc_score_3_HAPA_NA <- roc(train_regdata$subtMRCA_1stDesc, HAPA_NA_logit_P) #AUC score

jpeg(filename = "/01_GenDisFlu/Fig/PredModel/3_HAPA_NA_Model.jpeg", width = 170, height = 150, units = "mm", res = 520)

plot(roc_score_3_HAPA_NA, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="3. ROC curve of the HA + NA + PA GD Model")

dev.off()



lrtest(GD_2_3_HA_NA_model, GD_2_3_1_HAPA_NA_model)

#### 2.3.2 HA + NA*PA Model #####
GD_2_3_2_HA_NAPA_model <- glmer(subtMRCA_1stDesc ~ HA + N_A*PA + (1 | Vaccine_code), 
                                data=train_regdata, family=binomial(link="logit"))    ### Sig

summary(GD_2_3_2_HA_NAPA_model)

#### A. Predicted Score ####
HA_NAPA_logit_P = predict(GD_2_3_2_HA_NAPA_model, newdata = train_regdata, type = 'response' )
train_regdata$HA_NAPA_logit_P <- HA_NAPA_logit_P

#### B. ROC ####
roc(train_regdata$subtMRCA_1stDesc, HA_NAPA_logit_P) #AUC score

lrtest(GD_2_3_HA_NA_model, GD_2_3_2_HA_NAPA_model)


#### 2.3.3 HA*NS + NA Model #####
GD_2_3_3_HANS_NA_model <- glmer(subtMRCA_1stDesc ~ HA*NS + N_A + (1 | Vaccine_code), 
                                data=train_regdata, family=binomial(link="logit"))    ### Sig

summary(GD_2_3_3_HANS_NA_model)

#### A. Predicted Score ####
HANS_NA_logit_P = predict(GD_2_3_3_HANS_NA_model, newdata = train_regdata, type = 'response' )
train_regdata$HANS_NA_logit_P <- HANS_NA_logit_P

#### B. ROC ####
roc(train_regdata$subtMRCA_1stDesc, HANS_NA_logit_P) #AUC score

lrtest(GD_2_3_HA_NA_model, GD_2_3_3_HANS_NA_model)


#### 2.3.4 HA + NA*NS Model #####
GD_2_3_4_HA_NANS_model <- glmer(subtMRCA_1stDesc ~ HA + N_A*NS + (1 | Vaccine_code), 
                                data=train_regdata, family=binomial(link="logit"))    ### Sig

summary(GD_2_3_4_HA_NANS_model)

#### A. Predicted Score ####
HA_NANS_logit_P = predict(GD_2_3_4_HA_NANS_model, newdata = train_regdata, type = 'response' )
train_regdata$HA_NANS_logit_P <- HA_NANS_logit_P

#### B. ROC ####
roc(train_regdata$subtMRCA_1stDesc, HA_NANS_logit_P) #AUC score

lrtest(GD_2_3_HA_NA_model, GD_2_3_4_HA_NANS_model)

train_regdata %>% 
  ggplot() +
  geom_jitter(aes(x=factor(subtMRCA_1stDesc), y = HA_NANS_logit_P) , width = 0.25) +
  facet_wrap( ~ vaccineStrain, nrow = 4, ncol = 2) +
  theme_bw() 


#### 2.3.5 HA*PA + NA*NS Model #####
GD_2_3_5_HAPA_NANS_model <- glmer(subtMRCA_1stDesc ~ HA*PA + N_A*NS + (1 | Vaccine_code), 
                                data=train_regdata, family=binomial(link="logit"))    ### Sig

summary(GD_2_3_5_HAPA_NANS_model)

#### A. Predicted Score ####
HAPA_NANS_logit_P = predict(GD_2_3_5_HAPA_NANS_model, newdata = train_regdata, type = 'response' )
train_regdata$HAPA_NANS_logit_P <- HAPA_NANS_logit_P

roc_score_4_HAPA_NANS <- roc(train_regdata$subtMRCA_1stDesc, HAPA_NANS_logit_P) #AUC score

jpeg(filename = "/01_GenDisFlu/Fig/PredModel/4_HAPA_NANS_Model.jpeg", width = 170, height = 150, units = "mm", res = 520)

plot(roc_score_4_HAPA_NANS, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="4. ROC curve of the HA + NA + PA + NS GD Model")

dev.off()




train_regdata %>% 
  ggplot() +
  geom_jitter(aes(x=factor(subtMRCA_1stDesc), y = HAPA_NANS_logit_P) , width = 0.25) +
  facet_wrap( ~ vaccineStrain, nrow = 4, ncol = 2) +
  theme_bw() 



plot(train_regdata$HAPA_NANS_logit_P)

#### B. ROC ####
roc(train_regdata$subtMRCA_1stDesc, HAPA_NANS_logit_P) #AUC score

lrtest(GD_2_3_HA_NA_model, GD_2_3_5_HAPA_NANS_model)
lrtest(GD_2_3_4_HA_NANS_model, GD_2_3_5_HAPA_NANS_model)

plot(GD_2_3_5_HAPA_NANS_model)


#### 2.3.6 HA + NA +NS Model #####
GD_2_3_6_HA_NA_NS_model <- glmer(subtMRCA_1stDesc ~ HA + N_A + NS + (1 | Vaccine_code), 
                                  data=train_regdata, family=binomial(link="logit"))    ### Sig

summary(GD_2_3_6_HA_NA_NS_model)

#### A. Predicted Score ####
HA_NA_NS_logit_P = predict(GD_2_3_6_HA_NA_NS_model, newdata = train_regdata, type = 'response' )
train_regdata$HA_NA_NS_logit_P <- HA_NA_NS_logit_P

#### B. ROC ####
roc(train_regdata$subtMRCA_1stDesc, HA_NA_NS_logit_P) #AUC score

lrtest(GD_2_3_HA_NA_model, GD_2_3_6_HA_NA_NS_model)


############################ DIAGNOSIS ##################


########## Visualization  ######
windowsFonts(Arl = windowsFont("Arial"))

## Data import ##
Mdl_Select <- read.csv("Data/Premodel_Comparison.csv", header = T, na.strings = "")

Mdl_Select$Models <- factor(Mdl_Select$Names, 
                            levels = c("Ctrl Model 1", "Ctrl Model 2",
                                                "Model 1","Model 2",
                                                "Model 3","Model 4"), 
                            labels = c("RBD & 15A","Koel 2013","HA","+ NA",
                                       "+ PA","+ NS"))

jpeg(filename = "/01_GenDisFlu/Fig/PredModel/Modelfit_AIC.jpeg", width = 140, height = 90, units = "mm", res = 1200)

Mdl_Select %>% 
  ggplot(aes(x = Models, y = AIC, group = 1)) +
  geom_point(color = "#FBA01D", size = 4) +
  geom_line(color = "#FBA01D", linewidth = 1.7) +
  xlab("Models") +
  ylab("Akaike's Information Criterion (AIC)") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10, family = "Arl"),
        axis.text.x = element_text(size = 10, family = "Arl"),
        strip.text = element_text(size=10, family = "Arl")) 

dev.off()


jpeg(filename = "/01_GenDisFlu/Fig/PredModel/Modelfit_BIC.jpeg", width = 140, height = 90, units = "mm", res = 1200)

Mdl_Select %>% 
  ggplot(aes(x = Models, y = BIC, group = 1)) +
  geom_point(color = "#1F5C70", size = 4) +
  geom_line(color = "#1F5C70", linewidth = 1.7) +
  xlab("Models") +
  ylab("Bayesian Information Criterion (BIC)") +
  theme_bw()  + 
  theme(axis.text.y = element_text(size = 10, family = "Arl"),
                     axis.text.x = element_text(size = 10, family = "Arl"),
                     strip.text = element_text(size=10, family = "Arl")) 


dev.off()


jpeg(filename = "/01_GenDisFlu/Fig/PredModel/Modelfit_AUC.jpeg", width = 140, height = 90, units = "mm", res = 520)

Mdl_Select %>% 
  ggplot(aes(x = Models, y = ROC, group = 1)) +
  geom_point(color = "#B29476", size = 4) +
  geom_line(color = "#B29476", linewidth = 1.7) +
  xlab("Models") +
  ylab("Area under curve (AUC) of ROC") +
  theme_bw() + 
  theme(axis.text.y = element_text(size = 10, family = "Arl"),
        axis.text.x = element_text(size = 10, family = "Arl"),
        strip.text = element_text(size=10, family = "Arl")) 


dev.off()


