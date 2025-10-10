## Package ##
library(tidyverse)
library(lme4)
library(pROC)
library(lmtest)
library(corrplot)
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

subtMRCA_Neb_Prop <- train_regdata %>% 
  group_by(Vaccine_code, subtMRCA_Neibhor_Dom) %>% 
  summarise(n = n())



##### 1 Control model - RBD-15A #####
RBD_15A_model <- glmer(subtMRCA_Neibhor_Dom ~ HA_RBD + HA_15A + (1 | Vaccine_code), 
                    data=train_regdata, family=binomial(link="logit"))

# Summarize the model
summary(RBD_15A_model)


#### A. Predicted Score ####
RBD_15A_logit_P = predict(RBD_15A_model, newdata = train_regdata, type = 'response' )
train_regdata$RBD_15A_logit_P <- RBD_15A_logit_P

#### B. ROC ####
roc_score=roc(train_regdata$subtMRCA_Neibhor_Dom, RBD_15A_logit_P) #AUC score

jpeg(filename = "/01_GenDisFlu/Fig/PredModel/1_RBS15A_Model.jpeg", width = 170, height = 150, units = "mm", res = 520)

plot(roc_score, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE,
     print.thres=TRUE, print.thres.pch=19, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="1. ROC curve of the HA, RBS-15A Dominance Prediction Model")

dev.off()



##### 2 Control model - Koel #####
Koel_model <- glmer(subtMRCA_Neibhor_Dom ~ HA_Koel + (1 | Vaccine_code), 
                       data=train_regdata, family=binomial(link="logit"))

# Summarize the model
summary(Koel_model)

#### A. Predicted Score ####
Koel_logit_P = predict(Koel_model, newdata = train_regdata, type = 'response' )
train_regdata$Koel_logit_P <- Koel_logit_P



#### B. ROC ####
roc_score=roc(train_regdata$subtMRCA_Neibhor_Dom, Koel_logit_P) #AUC score

jpeg(filename = "/01_GenDisFlu/Fig/PredModel/2_Koel_Model.jpeg", width = 170, height = 150, units = "mm", res = 520)

plot(roc_score, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE,
     print.thres=TRUE, print.thres.pch=19, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="2. ROC curve of the HA, Koel Dominance Prediction Model")

dev.off()







##### 3 Genetic distance model - HA Only #####
## Model selection
GD_1_HA_model <- glmer(subtMRCA_Neibhor_Dom ~ HA + (1 | Vaccine_code), 
                    data=train_regdata, family=binomial(link="logit"))    

summary(GD_1_HA_model)

#### A. Predicted Score ####
HA_logit_P = predict(GD_1_HA_model, newdata = train_regdata, type = 'response' )
train_regdata$HA_logit_P <- HA_logit_P

jpeg(filename = "/01_GenDisFlu/Fig/PredModel/3_HA_Model.jpeg", width = 170, height = 150, units = "mm", res = 520)

#### B. ROC ####
roc_score=roc(train_regdata$subtMRCA_Neibhor_Dom, HA_logit_P) #AUC score
plot(roc_score, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE,
     print.thres=TRUE, print.thres.pch=19, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="3. ROC curve of the Only HA GD Prediction Model")

dev.off()



##### 4 Genetic distance model - HA + NA #####
## Model selection
GD_2_HA_NA_model <- glmer(subtMRCA_Neibhor_Dom ~ HA * N_A + (1 | Vaccine_code), 
                       data=train_regdata, family=binomial(link="logit"))    

summary(GD_2_HA_NA_model)

#### A. Predicted Score ####
HANA_logit_P = predict(GD_2_HA_NA_model, newdata = train_regdata, type = 'response' )
train_regdata$HANA_logit_P <- HANA_logit_P

#### B. ROC ####
roc_score=roc(train_regdata$subtMRCA_Neibhor_Dom, HANA_logit_P) #AUC score

jpeg(filename = "/01_GenDisFlu/Fig/PredModel/4_HAXNA_Model.jpeg", width = 170, height = 150, units = "mm", res = 520)

plot(roc_score, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE,
     print.thres=TRUE, print.thres.pch=19, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="4. ROC curve of the HA & NA GD Prediction Model")

dev.off()




## 5. Full GD Model selection
GD_Full_model <- glmer(subtMRCA_Neibhor_Dom ~ PB2 + PB1 + PA + HA + NP + N_A + M +  NS + (1 | Vaccine_code), 
                       data=train_regdata, family=binomial(link="logit"))                     

summary(GD_Full_model)


#### A. Predicted Score ####
GD_Full_logit_P = predict(GD_Full_model, newdata = train_regdata, type = 'response' )
train_regdata$GD_Full_logit_P <- GD_Full_logit_P

#### B. ROC ####
roc_score=roc(train_regdata$subtMRCA_Neibhor_Dom, GD_Full_logit_P) #AUC score

jpeg(filename = "/01_GenDisFlu/Fig/PredModel/5_FullGD_Model.jpeg", width = 170, height = 150, units = "mm", res = 520)

plot(roc_score, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE,
     print.thres=TRUE, print.thres.pch=19, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="5. ROC curve of the Full GD Dominance Prediction Model")

dev.off()


### 6. Backward Fit GD Final model ###
GD_BackFit_model <- glmer(subtMRCA_Neibhor_Dom ~ PB2 + HA + NP + N_A + M +  NS + (1 | Vaccine_code), 
                      data=train_regdata, family=binomial(link="logit"))                     

summary(GD_BackFit_model)

lrtest(GD_BackFit_model, GD_Full_model)

#### A. Predicted Score ####
GD_BackFit_logit_P = predict(GD_BackFit_model, newdata = train_regdata, type = 'response' )
train_regdata$GD_BackFit_logit_P <- GD_BackFit_logit_P

#### B. ROC ####
roc_score=roc(train_regdata$subtMRCA_Neibhor_Dom, GD_BackFit_logit_P) #AUC score
plot(roc_score, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE,
     print.thres=TRUE, print.thres.pch=19, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="6. ROC curve of the Backward Fit GD Dominance Prediction Model")


## HA ##
train_regdata %>% 
  ggplot(aes(x = Vaccine_code, y = HA, fill = as.factor(subtMRCA_Neibhor_Dom))) +
  geom_boxplot() +
  theme_bw()

## NA ##
train_regdata %>% 
  ggplot(aes(x = Vaccine_code, y = N_A, fill = as.factor(subtMRCA_Neibhor_Dom))) +
  geom_boxplot()+
  theme_bw()

## PA ##
train_regdata %>% 
  ggplot(aes(x = Vaccine_code, y = PA, fill = as.factor(subtMRCA_Neibhor_Dom))) +
  geom_boxplot()+
  theme_bw()

## M ##
train_regdata %>% 
  ggplot(aes(x = Vaccine_code, y = M, fill = as.factor(subtMRCA_Neibhor_Dom))) +
  geom_boxplot()+
  theme_bw()

## NP ##
train_regdata %>% 
  ggplot(aes(x = Vaccine_code, y = NP, fill = as.factor(subtMRCA_Neibhor_Dom))) +
  geom_boxplot()+
  theme_bw()




## CorMatrix for the selection
M <- cor(train_regdata[,c(3:10,24)])
testRes = cor.mtest(train_regdata[,c(3:10,24)], conf.level = 0.95)

## Correlation plot
corrplot(M, p.mat = testRes$p, method = 'circle', insig='blank', tl.srt = 0, 
         sig.level = 0.05, addCoef.col ='black', number.cex = 0.8)


### 7. Forward Fit GD Final model ###
GD_ForwFit_model <- glmer(subtMRCA_Neibhor_Dom ~ HA * PA +  (1 | Vaccine_code), 
                      data=train_regdata, family=binomial(link="logit"))                     

GD_ForwFit_model_2 <- glmer(subtMRCA_Neibhor_Dom ~ HA + PA + M + (1 | Vaccine_code), 
                          data=train_regdata, family=binomial(link="logit"))                     


summary(GD_ForwFit_model)
summary(GD_ForwFit_model_2)

lrtest(GD_ForwFit_model, GD_ForwFit_model_2)

#### A. Predicted Score ####
GD_ForwFit_logit_P = predict(GD_ForwFit_model, newdata = train_regdata, type = 'response' )
train_regdata$GD_ForwFit_logit_P <- GD_ForwFit_logit_P

#### B. ROC ####
roc_score=roc(train_regdata$subtMRCA_Neibhor_Dom, GD_ForwFit_logit_P) #AUC score
plot(roc_score, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE,
     print.thres=TRUE, print.thres.pch=19, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="7. ROC curve of the Forward Fit GD Dominance Prediction Model")




#### A. Predicted Score ####
GD_ForwFit_2_logit_P = predict(GD_ForwFit_model_2, newdata = train_regdata, type = 'response' )
train_regdata$GD_ForwFit_2_logit_P <- GD_ForwFit_2_logit_P

#### B. ROC ####
roc_score=roc(train_regdata$subtMRCA_Neibhor_Dom, GD_ForwFit_2_logit_P) #AUC score
plot(roc_score, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE,
     print.thres=TRUE, print.thres.pch=19, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="7. ROC curve of the Forward Fit GD Dominance Prediction Model 2 ")




chisq.test(train_regdata$subtMRCA_Neibhor_Dom, train_regdata$Reasrt) ## Deleterious 

#chisq.test(train_regdata$subtMRCA_Neibhor_Dom, train_regdata$CoEvo_PB2)
#chisq.test(train_regdata$subtMRCA_Neibhor_Dom, train_regdata$CoEvo_PB1) ## NO 1 Value 
#chisq.test(train_regdata$subtMRCA_Neibhor_Dom, train_regdata$CoEvo_PA) ## One 1 value
#chisq.test(train_regdata$subtMRCA_Neibhor_Dom, train_regdata$CoEvo_NP) ## One 1 value
#chisq.test(train_regdata$subtMRCA_Neibhor_Dom, train_regdata$CoEvo_NA) ## One 1 value
#chisq.test(train_regdata$subtMRCA_Neibhor_Dom, train_regdata$CoEvo_M) ## One 1 value
#chisq.test(train_regdata$subtMRCA_Neibhor_Dom, train_regdata$CoEvo_NS) ## Twp 1 value



## 5. Full All complex Model selection
All_Full_model <- glmer(subtMRCA_Neibhor_Dom ~ PB2 + PB1 + PA + HA + NP + N_A + M +  NS + Reasrt + (1 | Vaccine_code), 
                       data=train_regdata, family=binomial(link="logit"))                     

summary(All_Full_model)


### 6. Fit+Reassrt GD Final model ###
GD_Fit_Reasrt_model <- glmer(subtMRCA_Neibhor_Dom ~ PB2 + HA + NP + N_A + M +  NS + Reasrt + (1 | Vaccine_code), 
                      data=train_regdata, family=binomial(link="logit"))                     

summary(GD_Fit_Reasrt_model)

lrtest(GD_Fit_model, GD_Fit_Reasrt_model)



#### A. Predicted Score ####
GD_Full_logit_P = predict(GD_Full_model, newdata = train_regdata, type = 'response' )
train_regdata$GD_Full_logit_P <- GD_Full_logit_P

#### B. ROC ####
roc_score=roc(train_regdata$subtMRCA_Neibhor_Dom, GD_Full_logit_P) #AUC score
plot(roc_score, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE,
     print.thres=TRUE, print.thres.pch=19, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="6. ROC curve of the Full GD Dominance Prediction Model")




