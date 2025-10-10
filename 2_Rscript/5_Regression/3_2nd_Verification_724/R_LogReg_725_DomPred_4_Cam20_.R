## Package ##
library(tidyverse)
library(lme4)
library(pROC)
library(lmtest)
library(corrplot)

## Set working directory ##
setwd("01_GenDisFlu/")
#setwd("01_GenDisFlu/")

## Data import ##
AllG_NS <- read.csv("Data/Veri_725_2024/FullReg_724.csv", header = T, na.strings = "")

regdata <- AllG_NS[,c(1:8, 10, 12, 14, 16, 18, 20, 22)] 

colnames(regdata)[1]  <- c("ID")
colnames(regdata)[8]  <- c("PB2")
colnames(regdata)[9]  <- c("PB1")
colnames(regdata)[10]  <- c("PA")
colnames(regdata)[11]  <- c("HA")
colnames(regdata)[12]  <- c("NP")
colnames(regdata)[13]  <- c("N_A")
colnames(regdata)[14]  <- c("M")
colnames(regdata)[15]  <- c("NS")

### Subset data - Exclude MosCow ##

unique(regdata$vaccine_code)
## Regression model ##

## 4. Cam20 Prediction ##
### Subset data - Training ##
train_regdata <- regdata %>% 
  filter(vaccine_code != "Dar21" & vaccine_code != "Cam20" & vaccine_code != "Cal04" 
         & vaccine_code != "Fuj02" & vaccine_code != "Mos99")

### Subset data - Testing ##
test_regdata <- regdata %>% 
  filter(vaccine_code == "Cam20")

## set the same intercept as Swtz13
test_regdata$vaccine_code <- "HK19"


###### Model selection ########### 
## Regression model ##

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


##### 2.4 NA + PA  model ##### ---  Significant
GD_2_4_NA_PA_model <- glmer(Dom ~ N_A + PA + (1 | vaccine_code), 
                            data=train_regdata, family=binomial(link="logit"))    

summary(GD_2_4_NA_PA_model)
AIC(GD_2_4_NA_PA_model)


##### 2.4.1 NA * PA  model ##### ---  Error
GD_2_4_1_NAPA_model <- glmer(Dom ~ N_A * PA + (1 | vaccine_code), 
                            data=train_regdata, family=binomial(link="logit"))    

summary(GD_2_4_1_NAPA_model)
AIC(GD_2_4_1_NAPA_model)


##### 2.5 NA + NS  model ##### ---  Significant
GD_2_5_NA_NS_model <- glmer(Dom ~ N_A + NS + (1 | vaccine_code), 
                            data=train_regdata, family=binomial(link="logit"))    

summary(GD_2_5_NA_NS_model)
AIC(GD_2_5_NA_NS_model)


##### 2.5.1 NA *  NS  model ##### ---  Significant
GD_2_5_1_NANS_model <- glmer(Dom ~ N_A * NS + (1 | vaccine_code), 
                            data=train_regdata, family=binomial(link="logit"))    

summary(GD_2_5_1_NANS_model)
AIC(GD_2_5_1_NANS_model)



##### 3.1 HA + NA + NS model ##### ---  No Significant
GD_3_1_HA_NA_NS_model <- glmer(Dom ~ HA + N_A + NS + (1 | vaccine_code), 
                              data=train_regdata, family=binomial(link="logit"))    

summary(GD_3_1_HA_NA_NS_model)
AIC(GD_3_1_HA_NA_NS_model)

##### 3.2 HA + NA + PA model ##### ---  No Significant
GD_3_2_HA_NA_PA_model <- glmer(Dom ~ HA + N_A + PA + (1 | vaccine_code), 
                               data=train_regdata, family=binomial(link="logit"))    

summary(GD_3_2_HA_NA_PA_model)
AIC(GD_3_2_HA_NA_PA_model)


##### 3.3 NA + NS + PA model ##### ---  No Significant
GD_3_3_NA_PA_NS_model <- glmer(Dom ~ N_A + PA + NS +(1 | vaccine_code), 
                               data=train_regdata, family=binomial(link="logit"))    

summary(GD_3_3_NA_PA_NS_model)
AIC(GD_3_3_NA_PA_NS_model)

##### 3.3.1 NA + NS + NANS + PA model ##### ---  No Significant
GD_3_3_1_NANS_PA_model <- glmer(Dom ~ N_A * NS + PA  +(1 | vaccine_code), 
                               data=train_regdata, family=binomial(link="logit"))    

summary(GD_3_3_1_NANS_PA_model)
AIC(GD_3_3_1_NANS_PA_model)





######## Final model ##########
##### 3.1.1  NA*NS  + PA model ##### ---  Significant
GD_3_3_1_NANS_PA_model <- glmer(Dom ~ N_A * NS + PA  +(1 | vaccine_code), 
                                data=train_regdata, family=binomial(link="logit"))    

#### A. Predicted Score ####
NANS_PA_logit_P = predict(GD_3_3_1_NANS_PA_model, newdata = train_regdata, type = 'response' )
train_regdata$NANS_PA_logit_P <- NANS_PA_logit_P

#### B. ROC ####
roc_score_Train_NANS_PA <- roc(train_regdata$Dom, NANS_PA_logit_P) #AUC score

jpeg(filename = "01_GenDisFlu/Fig/PredModel/725_4_Cam20/4_Cam20_NANS_PA_Train.jpeg", width = 170, height = 150, units = "mm", res = 520)

plot(roc_score_Train_NANS_PA, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="ROC curve of the 4.Cam20 Model for Train set")

dev.off()

## Cross-validation using 11. Cam 20 ##
## Test data ##

#### A. Predicted Score ####
NANS_PA_logit_P_CV = predict(GD_3_3_1_NANS_PA_model, newdata = test_regdata, type = 'response' )
test_regdata$NANS_PA_logit_P_CV <- NANS_PA_logit_P_CV

#### B. ROC ####
roc_score_Test_NANS_PA <- roc(test_regdata$Dom, NANS_PA_logit_P_CV) #AUC score

jpeg(filename = "01_GenDisFlu/Fig/PredModel/725_4_Cam20/4_Cam20_NANS_PA_Test.jpeg", width = 170, height = 150, units = "mm", res = 520)

plot(roc_score_Test_NANS_PA, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="ROC curve of the 4.Cam20 Model for Test set")

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

jpeg(filename = "01_GenDisFlu/Fig/PredModel/725_4_Cam20/4_Cam20_HA_Train.jpeg", width = 170, height = 150, units = "mm", res = 520)

plot(roc_score_Train_HA, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="ROC curve of the 4.Cam20 Base Model for Train set")

dev.off()

## Cross-validation using 11.Cam20 ##
## Test data ##

#### A. Predicted Score ####
HA_logit_P_CV = predict(GD_1_HA_model, newdata = test_regdata, type = 'response' )
test_regdata$HA_logit_P_CV <- HA_logit_P_CV

#### B. ROC ####
roc_score_Test_HA <- roc(test_regdata$Dom, HA_logit_P_CV) #AUC score

jpeg(filename = "01_GenDisFlu/Fig/PredModel/725_4_Cam20/4_Cam20_HA_Test.jpeg", width = 170, height = 150, units = "mm", res = 520)

plot(roc_score_Test_HA, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="ROC curve of the 4.Cam20 Base Model for Test set")

dev.off()



