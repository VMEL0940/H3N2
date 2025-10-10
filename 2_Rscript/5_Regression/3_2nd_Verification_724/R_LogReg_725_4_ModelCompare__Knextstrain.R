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

regdata <- AllG_NS[,c(24, 2:8, 10, 12, 14, 16, 18, 20, 22)] 

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


############# Model 1 -  Dataset 1 HK 5 ###################
## 1. HK15 Prediction ##
### Subset data - Training ##
train_regdata_1 <- regdata %>% 
  filter(vaccine_code != "Dar21" & vaccine_code != "Cam20" & vaccine_code != "HK19" 
         & vaccine_code != "Kan17" & vaccine_code != "HK15")

### Subset data - Testing ##
test_regdata_1 <- regdata %>% 
  filter(vaccine_code == "HK15")

## set the same intercept as Swtz13
test_regdata_1$vaccine_code <- "Swtz13"


######## BAse model ##########
##### 1.HA only model ##### ---  Significant
GD_1_HA_model <- glmer(Dom ~ HA + (1 | vaccine_code), 
                       data=train_regdata, family=binomial(link="logit"))    

summary(GD_1_HA_model) ## AIC 274.4, BIC 286.6 ##

####### Model selection #####
GD_1_2_HA_NA_model <- glmer(Dom ~ HA + N_A + (1 | vaccine_code), 
                                data=train_regdata, family=binomial(link="logit"))  

summary(GD_1_2_HA_NA_model) ## AIC 264.1, BIC 280.3 ##

GD_1_3_HA_NANS_model <- glmer(Dom ~ HA + N_A*NS + (1 | vaccine_code), 
                            data=train_regdata, family=binomial(link="logit"))  

summary(GD_1_3_HA_NANS_model) ## AIC 260.1, BIC 284.3 ##  --- Selection

GD_1_4_HAPA_NANS_model <- glmer(Dom ~ HA*PA + N_A*NS + (1 | vaccine_code), 
                              data=train_regdata, family=binomial(link="logit"))  

summary(GD_1_4_HAPA_NANS_model) ## AIC 263.2, BIC 295.5 ##

######## Final model ##########
##### 3.1.2 HA + NA*NS model ##### ---  Significant
GD_1_3_HA_NANS_model <- glmer(Dom ~ HA + N_A*NS + (1 | vaccine_code), 
                              data=train_regdata, family=binomial(link="logit"))  

summary(GD_1_3_HA_NANS_model)

#### A. Predicted Score ####
HA_NANS_logit_P = predict(GD_1_3_HA_NANS_model, newdata = train_regdata_1, type = 'response' )
train_regdata_1$HA_NANS_logit_P <- HA_NANS_logit_P

#### B. ROC ####
roc_score_Train_HA_NANS <- roc(train_regdata_1$Dom, HA_NANS_logit_P) #AUC score
## 0.8246 ##

## Cross-validation using 9.HK15 ##
## Test data 1 ##

#### A. Predicted Score ####
HA_NANS_logit_P_CV = predict(GD_1_3_HA_NANS_model, newdata = test_regdata_1, type = 'response' )
test_regdata_1$HA_NANS_logit_P_CV <- HA_NANS_logit_P_CV

#### B. ROC ####
roc_score_Test_1_HA_NANS <- roc(test_regdata_1$Dom, HA_NANS_logit_P_CV) #AUC score
## 0.8946 ##



############# Model 2 -  Dataset 2 Kan 17 ###################

## 2. Kan17 Prediction ##
### Subset data - Training ##
train_regdata_2 <- regdata %>% 
  filter(vaccine_code != "Dar21" & vaccine_code != "Cam20" & vaccine_code != "HK19" 
         & vaccine_code != "Kan17" & vaccine_code != "Mos99")

### Subset data - Testing ##
test_regdata_2 <- regdata %>% 
  filter(vaccine_code == "Kan17")

## set the same intercept as Swtz13
test_regdata_2$vaccine_code <- "HK15"

## Cross-validation using 10.Kan17 ##
## Test data ##


######## BAse model ##########
##### 1.HA only model ##### ---  Significant
GD_2_HA_model <- glmer(Dom ~ HA + (1 | vaccine_code), 
                       data=train_regdata_2, family=binomial(link="logit"))    

summary(GD_2_HA_model) ## AIC 316.9, BIC 329.7 ##

####### Model selection #####
GD_2_2_HA_NA_model <- glmer(Dom ~ HA + N_A + (1 | vaccine_code), 
                            data=train_regdata_2, family=binomial(link="logit"))  

summary(GD_2_2_HA_NA_model) ## AIC 304.9, BIC 322.0 ##

GD_2_3_HA_NANS_model <- glmer(Dom ~ HA + N_A*NS + (1 | vaccine_code), 
                              data=train_regdata_2, family=binomial(link="logit"))  

summary(GD_2_3_HA_NANS_model) ## AIC 298.4, BIC 324.1 ##  --- Selection

GD_2_4_HAPA_NANS_model <- glmer(Dom ~ HA*PA + N_A*NS + (1 | vaccine_code), 
                                data=train_regdata_2, family=binomial(link="logit"))  

summary(GD_2_4_HAPA_NANS_model) ## AIC 297.3, BIC 331.6 ##

######## Final model ##########
##### 3.1.2 HA + NA*NS model ##### ---  Significant
GD_2_3_HA_NANS_model <- glmer(Dom ~ HA + N_A*NS + (1 | vaccine_code), 
                              data=train_regdata_2, family=binomial(link="logit"))  

summary(GD_2_3_HA_NANS_model)

#### A. Predicted Score ####
HA_NANS_logit_P = predict(GD_2_3_HA_NANS_model, newdata = train_regdata_2, type = 'response' )
train_regdata_2$HA_NANS_logit_P <- HA_NANS_logit_P

#### B. ROC ####
roc_score_Train_HA_NANS <- roc(train_regdata_2$Dom, HA_NANS_logit_P) #AUC score
## 0.819 ##

## Cross-validation using 10.Kan17  ##
## Test data 1 ##

#### A. Predicted Score ####
HA_NANS_logit_P_CV = predict(GD_2_3_HA_NANS_model, newdata = test_regdata_2, type = 'response' )
test_regdata_2$HA_NANS_logit_P_CV <- HA_NANS_logit_P_CV

#### B. ROC ####
roc_score_Test_HA_NANS <- roc(test_regdata_2$Dom, HA_NANS_logit_P_CV) #AUC score
## 0.5375 ##


############# Model 3 -  Dataset 3 HK19###################
## 3. HK19 Prediction ##
### Subset data - Training ##
train_regdata_3 <- regdata %>% 
  filter(vaccine_code != "Dar21" & vaccine_code != "Cam20" & vaccine_code != "HK19" 
         & vaccine_code != "Fuj02" & vaccine_code != "Mos99")

### Subset data - Testing ##
test_regdata_3 <- regdata %>% 
  filter(vaccine_code == "HK19")

## set the same intercept as Swtz13
test_regdata_3$vaccine_code <- "Kan17"


######## BAse model ##########
##### 1.HA only model ##### ---  Significant
GD_3_HA_model <- glmer(Dom ~ HA + (1 | vaccine_code), 
                       data=train_regdata_3, family=binomial(link="logit"))    

summary(GD_3_HA_model) ## AIC 311.4, BIC 324.1 ##

####### Model selection #####
GD_3_2_HA_NA_model <- glmer(Dom ~ HA + N_A + (1 | vaccine_code), 
                            data=train_regdata_3, family=binomial(link="logit"))  

summary(GD_3_2_HA_NA_model) ## AIC 300.2, BIC 317.1  ##

GD_3_3_HA_NANS_model <- glmer(Dom ~ HA + N_A*NS + (1 | vaccine_code), 
                              data=train_regdata_3, family=binomial(link="logit"))  

summary(GD_3_3_HA_NANS_model) ## AIC 291.2, BIC 316.6 ##  --- Selection

GD_3_4_HAPA_NANS_model <- glmer(Dom ~ HA*PA + N_A*NS + (1 | vaccine_code), 
                                data=train_regdata_3, family=binomial(link="logit"))  

summary(GD_3_4_HAPA_NANS_model) ## AIC 289.4, BIC 323.3 ##

######## Final model ##########
##### 3.1.2 HA + NA*NS model ##### ---  Significant
GD_3_3_HA_NANS_model <- glmer(Dom ~ HA + N_A*NS + (1 | vaccine_code), 
                              data=train_regdata_3, family=binomial(link="logit"))  

summary(GD_3_3_HA_NANS_model)

#### A. Predicted Score ####
HA_NANS_logit_P = predict(GD_3_3_HA_NANS_model, newdata = train_regdata_3, type = 'response' )
train_regdata_3$HA_NANS_logit_P <- HA_NANS_logit_P

#### B. ROC ####
roc_score_Train_HA_NANS <- roc(train_regdata_3$Dom, HA_NANS_logit_P) #AUC score
## 0.8144 ##

## Cross-validation using 11.HK19  ##
## Test data 1 ##

#### A. Predicted Score ####
HA_NANS_logit_P_CV = predict(GD_3_3_HA_NANS_model, newdata = test_regdata_3, type = 'response' )
test_regdata_3$HA_NANS_logit_P_CV <- HA_NANS_logit_P_CV

#### B. ROC ####
roc_score_Test_HA_NANS <- roc(test_regdata_3$Dom, HA_NANS_logit_P_CV) #AUC score
## 0.6296 ##


############# Model 4 -  Dataset 4 Cam20 ###################

## 4. Cam20 Prediction ##
### Subset data - Training ##
train_regdata_4 <- regdata %>% 
  filter(vaccine_code != "Dar21" & vaccine_code != "Cam20" & vaccine_code != "Cal04" 
         & vaccine_code != "Fuj02" & vaccine_code != "Mos99")

### Subset data - Testing ##
test_regdata_4 <- regdata %>% 
  filter(vaccine_code == "Cam20")

## set the same intercept as Swtz13
test_regdata_4$vaccine_code <- "HK19"


######## BAse model ##########
##### 1.HA only model ##### ---  Significant
GD_4_HA_model <- glmer(Dom ~ HA + (1 | vaccine_code), 
                       data=train_regdata_4, family=binomial(link="logit"))    

summary(GD_4_HA_model) ## AIC 335.9, BIC 348.7 ##

####### Model selection #####
GD_4_2_HA_NA_model <- glmer(Dom ~ HA + N_A + (1 | vaccine_code), 
                            data=train_regdata_4, family=binomial(link="logit"))  

summary(GD_4_2_HA_NA_model) ## AIC 329.0, BIC 346.0  ##

GD_4_3_HA_NANS_model <- glmer(Dom ~ HA + N_A*NS + (1 | vaccine_code), 
                              data=train_regdata_4, family=binomial(link="logit"))  

summary(GD_4_3_HA_NANS_model) ## AIC 316.9, BIC 342.4 ##  --- Selection

GD_4_4_HAPA_NANS_model <- glmer(Dom ~ HA*PA + N_A*NS + (1 | vaccine_code), 
                                data=train_regdata_4, family=binomial(link="logit"))  

summary(GD_4_4_HAPA_NANS_model) ## AIC 312.6 , BIC 346.7 ##

######## Final model ##########
##### 3.1.2 HA + NA*NS model ##### ---  Significant
GD_4_3_HA_NANS_model <- glmer(Dom ~ HA + N_A*NS + (1 | vaccine_code), 
                              data=train_regdata_4, family=binomial(link="logit"))  

summary(GD_4_3_HA_NANS_model)

#### A. Predicted Score ####
HA_NANS_logit_P = predict(GD_4_3_HA_NANS_model, newdata = train_regdata_4, type = 'response' )
train_regdata_4$HA_NANS_logit_P <- HA_NANS_logit_P

#### B. ROC ####
roc_score_Train_HA_NANS <- roc(train_regdata_4$Dom, HA_NANS_logit_P) #AUC score
## 0.8187 ##

## Cross-validation using 11.HK19  ##
## Test data 1 ##

#### A. Predicted Score ####
HA_NANS_logit_P_CV = predict(GD_4_3_HA_NANS_model, newdata = test_regdata_4, type = 'response' )
test_regdata_4$HA_NANS_logit_P_CV <- HA_NANS_logit_P_CV

#### B. ROC ####
roc_score_Test_HA_NANS <- roc(test_regdata_4$Dom, HA_NANS_logit_P_CV) #AUC score
## 0.8095 ##



############# Model 5 -  Dataset 5 Dar21 ###################

## 4. Cam20 Prediction ##
### Subset data - Training ##
train_regdata_5 <- regdata %>% 
  filter(vaccine_code != "Dar21" & vaccine_code != "Wis05" & vaccine_code != "Cal04" 
         & vaccine_code != "Fuj02" & vaccine_code != "Mos99")

### Subset data - Testing ##
test_regdata_5 <- regdata %>% 
  filter(vaccine_code == "Dar21")

## set the same intercept as Cam20
test_regdata_5$vaccine_code <- "Cam20"


######## BAse model ##########
##### 1.HA only model ##### ---  Significant
GD_5_HA_model <- glmer(Dom ~ HA + (1 | vaccine_code), 
                       data=train_regdata_5, family=binomial(link="logit"))    

summary(GD_5_HA_model) ## AIC 345.5, BIC 358.3 ##

####### Model selection #####
GD_5_2_HA_NA_model <- glmer(Dom ~ HA + N_A + (1 | vaccine_code), 
                            data=train_regdata_5, family=binomial(link="logit"))  

summary(GD_5_2_HA_NA_model) ## AIC 336.5, BIC 353.6  ##

GD_5_3_HA_NANS_model <- glmer(Dom ~ HA + N_A*NS + (1 | vaccine_code), 
                              data=train_regdata_5, family=binomial(link="logit"))  

summary(GD_5_3_HA_NANS_model) ## AIC 325.6, BIC 351.3 ## 

GD_5_4_HAPA_NANS_model <- glmer(Dom ~ HA*PA + N_A*NS + (1 | vaccine_code), 
                                data=train_regdata_5, family=binomial(link="logit"))  

summary(GD_5_4_HAPA_NANS_model) ## AIC 315.7, BIC 349.9 ## --- Selection

######## Final model ##########
##### 3.1.2 HA + NA*NS model ##### ---  Significant
GD_5_4_HAPA_NANS_model <- glmer(Dom ~ HA*PA + N_A*NS + (1 | vaccine_code), 
                                data=train_regdata_5, family=binomial(link="logit"))  

summary(GD_5_4_HAPA_NANS_model)

#### A. Predicted Score ####
HAPA_NANS_logit_P = predict(GD_5_4_HAPA_NANS_model, newdata = train_regdata_5, type = 'response' )
train_regdata_5$HAPA_NANS_logit_P <- HAPA_NANS_logit_P

#### B. ROC ####
roc_score_Train_HAPA_NANS <- roc(train_regdata_5$Dom, HAPA_NANS_logit_P) #AUC score
## 0.8428 ##

## Cross-validation using 11.HK19  ##
## Test data 1 ##

#### A. Predicted Score ####
HAPA_NANS_logit_P_CV = predict(GD_5_4_HAPA_NANS_model, newdata = test_regdata_5, type = 'response' )
test_regdata_5$HA_NANS_logit_P_CV <- HAPA_NANS_logit_P_CV



Predict_all <- rbind(train_regdata_1[c(1,3,16)], test_regdata_1[c(1,16)], test_regdata_2[c(1,16)], test_regdata_3[c(1,16)], test_regdata_4[c(1,16)], test_regdata_5[c(1,16)])
write.csv(train_regdata_1[c(1,3,16)], "train_regdata_1.csv")
write.csv(test_regdata_1[c(1,3,16)], "test_regdata_1.csv")
write.csv(test_regdata_2[c(1,3,16)], "test_regdata_2.csv")
write.csv(test_regdata_3[c(1,3,16)], "test_regdata_3.csv")
write.csv(test_regdata_4[c(1,3,16)], "test_regdata_4.csv")
write.csv(test_regdata_5[c(1,3,16)], "test_regdata_5.csv")
