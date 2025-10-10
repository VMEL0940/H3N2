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
colnames(regdata)[25] <- "Dom"

### Subset data - Exclude HK15 ##
train_regdata <- regdata %>% 
  filter(vaccineStrain != "EPI686117_A_Hong_Kong_15611_2015_NA_NA_NA")

train_regdata <- train_regdata %>% 
  filter(ID != "EPI103320_A_Moscow_10_1999")
train_regdata$vaccine_code <- factor(train_regdata$vaccineStrain, 
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

###### Model selection ########### 
## Regression model ##


##### 0_1.HA-RBD&15A only model ##### ---  Significant
GD_0_1_RBD15A_model <- glmer(Dom ~ HA_RBD + HA_15A +(1 | vaccine_code), 
                       data=train_regdata, family=binomial(link="logit"))    

summary(GD_0_1_RBD15A_model)
AIC(GD_0_1_RBD15A_model)


#### A. Predicted Score ####
RBD15A_logit_P = predict(GD_0_1_RBD15A_model, newdata = train_regdata, type = 'response' )
train_regdata$RBD15A_logit_P <- RBD15A_logit_P

#### B. ROC ####
roc_score_Train_RBD15A <- roc(train_regdata$Dom, RBD15A_logit_P) #AUC score

jpeg(filename = "01_GenDisFlu/Fig/PredModel/0_601Model/0_601_HA_RBD15A_Train.jpeg", width = 170, height = 150, units = "mm", res = 520)
svg(filename = "01_GenDisFlu/Fig/PredModel/0_601Model/0_601_1_HA_RBD15A_Train.svg", width = 17, height = 15)

plot(roc_score_Train_RBD15A, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="B1. ROC curve of the HA, RBD-15A Prediction model")

dev.off()



##### 0_2.HA-Koel only model ##### ---  Significant
GD_0_2_Koel_model <- glmer(Dom ~ HA_Koel +(1 | vaccine_code), 
                             data=train_regdata, family=binomial(link="logit"))    

summary(GD_0_2_Koel_model)
AIC(GD_0_2_Koel_model)


#### A. Predicted Score ####
Koel_logit_P = predict(GD_0_2_Koel_model, newdata = train_regdata, type = 'response' )
train_regdata$Koel_logit_P <- Koel_logit_P

#### B. ROC ####
roc_score_Train_Koel <- roc(train_regdata$Dom, Koel_logit_P) #AUC score

#jpeg(filename = "01_GenDisFlu/Fig/PredModel/0_601Model/0_601_1_HA_Koel_Train.svg", width = 170, height = 150, units = "mm", res = 520)
svg(filename = "01_GenDisFlu/Fig/PredModel/0_601Model/0_601_2_HA_Koel_Train.svg", width = 17, height = 15)

plot(roc_score_Train_Koel, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="B2. ROC curve of the HA, Koel Prediction model")

dev.off()



##### 1.HA only model ##### ---  Significant
GD_1_HA_model <- glmer(Dom ~ HA + (1 | vaccine_code), 
                       data=train_regdata, family=binomial(link="logit"))    

summary(GD_1_HA_model)
AIC(GD_1_HA_model)


##### 1.1 NA only model ##### ---  Significant
GD_1_1_NA_model <- glmer(Dom ~ N_A + (1 | vaccine_code), 
                         data=train_regdata, family=binomial(link="logit"))    

summary(GD_1_1_NA_model)
AIC(GD_1_1_NA_model)


##### 2. HA + 1  model #####

##### 2.1 HA + NA  model ##### ---  Significant
GD_2_1_HA_NA_model <- glmer(Dom ~ HA + N_A + (1 | vaccine_code), 
                            data=train_regdata, family=binomial(link="logit"))    

summary(GD_2_1_HA_NA_model)
AIC(GD_2_1_HA_NA_model)


##### 2.1.1 HA X NA  model #####  --- Error 
GD_2_1_1_HANA_model <- glmer(Dom ~ HA * N_A + (1 | vaccine_code), 
                             data=train_regdata, family=binomial(link="logit"))    

summary(GD_2_1_1_HANA_model)
AIC(GD_2_1_1_HANA_model)


##### 2.2 HA + PA  model ##### ---  Error
GD_2_2_HA_PA_model <- glmer(Dom ~ HA + PA + (1 | vaccine_code), 
                            data=train_regdata, family=binomial(link="logit"))    

summary(GD_2_2_HA_PA_model)
AIC(GD_2_2_HA_PA_model)


##### 2.3 HA + NS  model ##### ---  Significant
GD_2_3_HA_NS_model <- glmer(Dom ~ HA + NS + (1 | vaccine_code), 
                            data=train_regdata, family=binomial(link="logit"))    

summary(GD_2_3_HA_NS_model)
AIC(GD_2_3_HA_NS_model)


##### 3.1 HA + NA + NS model ##### ---  No Significant
GD_3_1_HA_NA_NS_model <- glmer(Dom ~ HA + N_A + NS + (1 | vaccine_code), 
                               data=train_regdata, family=binomial(link="logit"))    

summary(GD_3_1_HA_NA_NS_model)
AIC(GD_3_1_HA_NA_NS_model)

##### 3.1.1 HA + NA*NS model ##### ---  No Significant
GD_3_1_1_HA_NANS_model <- glmer(Dom ~ HA + N_A * NS + (1 | vaccine_code), 
                               data=train_regdata, family=binomial(link="logit"))    

summary(GD_3_1_1_HA_NANS_model)
AIC(GD_3_1_1_HA_NANS_model)

##### 3.1.2 HA*NS + NA model ##### ---  No Significant
GD_3_1_2_HANS_NA_model <- glmer(Dom ~ HA * NS + N_A  + (1 | vaccine_code), 
                                data=train_regdata, family=binomial(link="logit"))    

summary(GD_3_1_2_HANS_NA_model)
AIC(GD_3_1_2_HANS_NA_model)



##### 3.2 HA + NA + PA model ##### ---  No Significant
GD_3_2_HA_NA_PA_model <- glmer(Dom ~ HA + N_A + PA + (1 | vaccine_code), 
                               data=train_regdata, family=binomial(link="logit"))    

summary(GD_3_2_HA_NA_PA_model)
AIC(GD_3_2_HA_NA_PA_model)


##### 3.2.1 HA *PA + NA  model ##### ---  No Significant
GD_3_2_1_HAPA_NA_model <- glmer(Dom ~ HA * PA + N_A  + (1 | vaccine_code), 
                               data=train_regdata, family=binomial(link="logit"))    

summary(GD_3_2_1_HAPA_NA_model)
AIC(GD_3_2_1_HAPA_NA_model)

##### 3.2.2 HA  + NA*PA  model ##### ---  No Significant
GD_3_2_2_HA_NAPA_model <- glmer(Dom ~ HA  + N_A * PA + (1 | vaccine_code), 
                                data=train_regdata, family=binomial(link="logit"))    

summary(GD_3_2_2_HA_NAPA_model)
AIC(GD_3_2_2_HA_NAPA_model)



##### 4.1 HA + NA*NS + PA model ##### ---  No Significant
GD_4_1_HA_NANS_PA_model <- glmer(Dom ~ HA + N_A * NS + PA + (1 | vaccine_code), 
                                data=train_regdata, family=binomial(link="logit"))    

summary(GD_4_1_HA_NANS_PA_model)
AIC(GD_4_1_HA_NANS_PA_model)

##### 4.2 HA*PA + NA*NS  model ##### ---  No Significant
GD_4_2_HAPA_NANS_model <- glmer(Dom ~ HA * PA + N_A * NS  + (1 | vaccine_code), 
                                 data=train_regdata, family=binomial(link="logit"))    

summary(GD_4_2_HAPA_NANS_model)
AIC(GD_4_2_HAPA_NANS_model)



######## Final model ##########
##### 3.1.1  HA + NA*NS model ##### ---  Significant
GD_3_1_1_HA_NANS_model <- glmer(Dom ~ HA + N_A * NS + (1 | vaccine_code), 
                                data=train_regdata, family=binomial(link="logit"))    

#### A. Predicted Score ####
HA_NANS_logit_P = predict(GD_3_1_1_HA_NANS_model, newdata = train_regdata, type = 'response' )
train_regdata$HA_NANS_logit_P <- HA_NANS_logit_P

#### B. ROC ####
roc_score_Train_HA_NANS <- roc(train_regdata$Dom, HA_NANS_logit_P) #AUC score

jpeg(filename = "01_GenDisFlu/Fig/PredModel/0_601Model/0_601_HA_NANS_Train.jpeg", width = 170, height = 150, units = "mm", res = 520)

plot(roc_score_Train_HA_NANS, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="ROC curve of the Regression Model for Train set")

dev.off()


######## BAse model ##########
##### 1.HA only model ##### ---  Significant
GD_1_HA_model <- glmer(Dom ~ HA + (1 | vaccine_code), 
                       data=train_regdata, family=binomial(link="logit"))    

#### A. Predicted Score ####
HA_logit_P = predict(GD_1_HA_model, newdata = train_regdata, type = 'response' )
train_regdata$HA_logit_P <- HA_logit_P

#### B. ROC ####
roc_score_Train_HA <- roc(train_regdata$Dom, HA_logit_P) #AUC score

jpeg(filename = "01_GenDisFlu/Fig/PredModel/0_601Model/0_601_HA_Train.jpeg", width = 170, height = 150, units = "mm", res = 520)

plot(roc_score_Train_HA, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="ROC curve of the HA Model for Train set")

dev.off()




######## BAse model ##########
##### 1.HA only model ##### ---  Significant
GD_1_1_HA_Koel_model <- glmer(Dom ~ HA_Koel + (1 | vaccine_code), 
                       data=train_regdata, family=binomial(link="logit"))    

#### A. Predicted Score ####
HA_Koel_logit_P = predict(GD_1_1_HA_Koel_model, newdata = train_regdata, type = 'response' )
train_regdata$HA_Koel_logit_P <- HA_Koel_logit_P

#### B. ROC ####
roc_score_Train_HA_Koel <- roc(train_regdata$Dom, HA_Koel_logit_P) #AUC score

jpeg(filename = "01_GenDisFlu/Fig/PredModel/0_601Model/0_601_1_HA_Koel_Train.jpeg", width = 170, height = 150, units = "mm", res = 520)

plot(roc_score_Train_HA_Koel, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="ROC curve of the HA-Koel Model for Train set")

dev.off()



######## BAse model ##########
##### 1.HA only model ##### ---  Significant
GD_1_2_HA_RBD15_model <- glmer(Dom ~ HA_RBD + HA_15A + (1 | vaccine_code), 
                              data=train_regdata, family=binomial(link="logit"))    

#### A. Predicted Score ####
HA_RBD15_logit_P = predict(GD_1_2_HA_RBD15_model, newdata = train_regdata, type = 'response' )
train_regdata$HA_RBD15_logit_P <- HA_RBD15_logit_P

#### B. ROC ####
roc_score_Train_HA_RBD15 <- roc(train_regdata$Dom, HA_RBD15_logit_P) #AUC score

jpeg(filename = "01_GenDisFlu/Fig/PredModel/0_601Model/0_601_2_HA_RBD15_Train.jpeg", width = 170, height = 150, units = "mm", res = 520)

plot(roc_score_Train_HA_RBD15, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="ROC curve of the HA-RBD&15A Model for Train set")

dev.off()



################ Model Comparison --- AIC + Pred distribution #####

##### 0_1.HA-RBD&15A only model ##### ---  Significant
GD_0_1_RBD15A_model <- glmer(Dom ~ HA_RBD + HA_15A +(1 | vaccine_code), 
                             data=train_regdata, family=binomial(link="logit"))    

summary(GD_0_1_RBD15A_model)
AIC(GD_0_1_RBD15A_model)

HA_RBD15_logit_P = predict(GD_0_1_RBD15A_model, newdata = train_regdata, type = 'response' )
train_regdata$HA_RBD15_logit_P <- HA_RBD15_logit_P

roc_score_Train_HA_RBD15 <- roc(train_regdata$Dom, HA_RBD15_logit_P) #AUC score


##### 0_2.HA-Koel only model ##### ---  Significant
GD_0_2_Koel_model <- glmer(Dom ~ HA_Koel +(1 | vaccine_code), 
                           data=train_regdata, family=binomial(link="logit"))    

summary(GD_0_2_Koel_model)
AIC(GD_0_2_Koel_model)

HA_Koel_logit_P = predict(GD_0_2_Koel_model, newdata = train_regdata, type = 'response' )
train_regdata$HA_Koel_logit_P <- HA_Koel_logit_P

roc_score_Train_HA_Koel <- roc(train_regdata$Dom, HA_Koel_logit_P) #AUC score

##### 1.HA only model ##### ---  Significant
GD_1_HA_model <- glmer(Dom ~ HA + (1 | vaccine_code), 
                       data=train_regdata, family=binomial(link="logit"))    

summary(GD_1_HA_model)
AIC(GD_1_HA_model)

HA_logit_P = predict(GD_1_HA_model, newdata = train_regdata, type = 'response' )
train_regdata$HA_logit_P <- HA_logit_P

roc_score_Train_HA <- roc(train_regdata$Dom, HA_logit_P) #AUC score

##### 1.1 NA only model ##### ---  Significant
GD_1_1_NA_model <- glmer(Dom ~ N_A + (1 | vaccine_code), 
                         data=train_regdata, family=binomial(link="logit"))    

summary(GD_1_1_NA_model)
AIC(GD_1_1_NA_model)

NA_logit_P = predict(GD_1_1_NA_model, newdata = train_regdata, type = 'response' )
train_regdata$NA_logit_P <- NA_logit_P

roc_score_Train_NA <- roc(train_regdata$Dom, NA_logit_P) #AUC score


##### 2.1 HA + NA  model ##### ---  Significant
GD_2_1_HA_NA_model <- glmer(Dom ~ HA + N_A + (1 | vaccine_code), 
                            data=train_regdata, family=binomial(link="logit"))    

summary(GD_2_1_HA_NA_model)
AIC(GD_2_1_HA_NA_model)

HA_NA_logit_P = predict(GD_2_1_HA_NA_model, newdata = train_regdata, type = 'response' )
train_regdata$HA_NA_logit_P <- HA_NA_logit_P

roc_score_Train_HA_NA <- roc(train_regdata$Dom, HA_NA_logit_P) #AUC score

######## Final model ##########
##### 3.1.1  HA + NA*NS model ##### ---  Significant
GD_3_1_1_HA_NANS_model <- glmer(Dom ~ HA + N_A * NS + (1 | vaccine_code), 
                                data=train_regdata, family=binomial(link="logit"))    

#### A. Predicted Score ####
HA_NANS_logit_P = predict(GD_3_1_1_HA_NANS_model, newdata = train_regdata, type = 'response' )
train_regdata$HA_NANS_logit_P <- HA_NANS_logit_P

roc_score_Train_HA_NANS <- roc(train_regdata$Dom, HA_NANS_logit_P) #AUC score


##### 4.1 HA + NA*NS + PA model ##### ---  No Significant
GD_4_1_HA_NANS_PA_model <- glmer(Dom ~ HA + N_A * NS + PA + (1 | vaccine_code), 
                                 data=train_regdata, family=binomial(link="logit"))    

summary(GD_4_1_HA_NANS_PA_model)
AIC(GD_4_1_HA_NANS_PA_model)


#### A. Predicted Score ####
HA_NANS_PA_logit_P = predict(GD_4_1_HA_NANS_PA_model, newdata = train_regdata, type = 'response' )
train_regdata$HA_NANS_PA_logit_P <- HA_NANS_PA_logit_P


roc_score_Train_HA_NANS_PA <- roc(train_regdata$Dom, HA_NANS_PA_logit_P) #AUC score


## Subset prediction comparison ##

PredCompare <- train_regdata[c(1, 25:33)]

colnames(PredCompare)[2] <- "Dominance"
colnames(PredCompare)[4] <- "HA-RBD+15A"
colnames(PredCompare)[5] <- "HA-Koel"
colnames(PredCompare)[6] <- "HA-NS"
colnames(PredCompare)[7] <- "NA-NS"
colnames(PredCompare)[8] <- "HA+NA-NS"
colnames(PredCompare)[9] <- "HA+NA+NS-NS"
colnames(PredCompare)[10] <- "HA+NA+NS+PA-NS"

PredCompare_gat <-  PredCompare %>% 
  gather("Model", "Pred_Val", 4:10) 

PredCompare_gat$Model <- factor(PredCompare_gat$Model, levels = c("HA-RBD+15A",
                                                                   "HA-Koel",
                                                                   "HA-NS",
                                                                   "NA-NS",
                                                                   "HA+NA-NS",
                                                                   "HA+NA+NS-NS",
                                                                   "HA+NA+NS+PA-NS"))

PredCompare_gat$Dominance <- factor(PredCompare_gat$Dominance, levels = c(0,1), labels = c("Extinct","Dominant"))


VP3Col <- c("#333333","#96fd29")


jpeg(filename = "01_GenDisFlu/Fig/PredModel/0_601Model/0_601_PredScore_Compare.jpeg", width = 220, height = 120, units = "mm", res = 720)

PredCompare_gat%>% 
  ggplot() + theme_bw() +
  geom_point(position = position_jitterdodge(jitter.width = 0.2,
                                             dodge.width = 0.7), aes(x= Model, y = Pred_Val, color = as.factor(Dominance))) +
  geom_boxplot(aes(x= Model, y = Pred_Val, color = as.factor(Dominance)), size = 0.8, outlier.shape = NA, alpha = 0) +
  scale_color_manual(values = VP3Col, name="Dominance", labels=c('Extinct', 'Dominant')) + 
  xlab("Models") +
  ylab("Prediction score")

dev.off()


### Comparing AIC and ROC of 7 Models  ###


No <- 1:7
Models <- c("HA-RBD+15A","HA-Koel",
            "HA-NS","NA-NS",
            "HA+NA-NS","HA+NA+NS-NS","HA+NA+NS+PA-NS")
AIC <- c(AIC(GD_0_1_RBD15A_model), AIC(GD_0_2_Koel_model),
         AIC(GD_1_HA_model), AIC(GD_1_1_NA_model),
         AIC(GD_2_1_HA_NA_model),
         AIC(GD_3_1_1_HA_NANS_model), AIC(GD_4_1_HA_NANS_PA_model))

AUC_ROC <- c(auc(roc_score_Train_HA_RBD15), auc(roc_score_Train_HA_Koel),
         auc(roc_score_Train_HA), auc(roc_score_Train_NA),
         auc(roc_score_Train_HA_NA),
         auc(roc_score_Train_HA_NANS), auc(roc_score_Train_HA_NANS_PA))

Model_0_gof <- data.frame(No, Models, AIC, AUC_ROC)

Model_0_gof$Models <- factor(Model_0_gof$Models, levels = c("HA-RBD+15A",
                                                                  "HA-Koel",
                                                                  "HA-NS",
                                                                  "NA-NS",
                                                                  "HA+NA-NS",
                                                                  "HA+NA+NS-NS",
                                                                  "HA+NA+NS+PA-NS"))



coeff <- 250
AICcolor <- "#FBA01D"
AUCcolor <- "#1F5C70"

jpeg(filename = "01_GenDisFlu/Fig/PredModel/0_601Model/0_601_AICAUC_Compare.jpeg", width = 220, height = 100, units = "mm", res = 720)

Model_0_gof %>% 
  ggplot(aes(x=Models)) +
  geom_line( aes(y=AIC, group=1), size = 1.1, color = AICcolor) + 
  geom_point( aes(y=AIC), size = 2.3, color = AICcolor) +
  geom_line( aes(y=AUC_ROC*coeff, group=1),size = 1.1,color = AUCcolor) + 
  geom_point( aes(y=AUC_ROC*coeff), size = 2.3, color = AUCcolor) +
  scale_y_continuous(name = "AIC",
        sec.axis = sec_axis(trans=~./coeff, name="AUC of ROC")) + 
  theme_bw() +
  theme(
    axis.title.y = element_text(color = AICcolor, size=12),
    axis.title.y.right = element_text(color = AUCcolor, size=12)
  ) 

dev.off()
