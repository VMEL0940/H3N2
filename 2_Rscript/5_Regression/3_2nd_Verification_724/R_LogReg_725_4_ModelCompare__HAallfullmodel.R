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


############# Model 1 -  Dataset 1 HK 15 ###################
## 1. HK15 Prediction ##
### Subset data - Training ##
train_regdata <- regdata %>% 
  filter(vaccine_code != "Dar21" & vaccine_code != "Cam20" & vaccine_code != "HK19" 
         & vaccine_code != "Kan17" & vaccine_code != "HK15")

### Subset data - Testing ##
test_regdata_1 <- regdata %>% 
  filter(vaccine_code == "HK15")

### Subset data - Testing ##
test_regdata_2 <- regdata %>% 
  filter(vaccine_code == "Kan17")

### Subset data - Testing ##
test_regdata_3 <- regdata %>% 
  filter(vaccine_code == "HK19")

### Subset data - Testing ##
test_regdata_4 <- regdata %>% 
  filter(vaccine_code == "Cam20")

## set the same intercept as Swtz13
test_regdata_1$vaccine_code <- "Swtz13"
test_regdata_2$vaccine_code <- "Swtz13"
test_regdata_3$vaccine_code <- "Swtz13"
test_regdata_4$vaccine_code <- "Swtz13"


######## BAse model ##########
##### 1.HA only model ##### ---  Significant
GD_1_HA_model <- glmer(Dom ~ HA + (1 | vaccine_code), 
                       data=train_regdata, family=binomial(link="logit"))    

#### A. Predicted Score ####
HA_logit_P = predict(GD_1_HA_model, newdata = train_regdata, type = 'response' )
train_regdata$HA_logit_P <- HA_logit_P

#### B. ROC ####
roc_score_Train_HA <- roc(train_regdata$Dom, HA_logit_P) #AUC score
## 0.7649 ##

## Cross-validation using 9.HK15 ##
## Test data 1 ##

#### A. Predicted Score ####
HA_logit_P_CV = predict(GD_1_HA_model, newdata = test_regdata_1, type = 'response' )
test_regdata_1$HA_logit_P_CV <- HA_logit_P_CV

#### B. ROC ####
roc_score_Test_1_HA <- roc(test_regdata_1$Dom, HA_logit_P_CV) #AUC score
## 0.8571 ##

## Test data 2 ##

#### A. Predicted Score ####
HA_logit_P_CV = predict(GD_1_HA_model, newdata = test_regdata_2, type = 'response' )
test_regdata_2$HA_logit_P_CV <- HA_logit_P_CV

#### B. ROC ####
roc_score_Test_2_HA <- roc(test_regdata_2$Dom, HA_logit_P_CV) #AUC score
## 0.8438 ##

## Test data 3 ##

#### A. Predicted Score ####
HA_logit_P_CV = predict(GD_1_HA_model, newdata = test_regdata_3, type = 'response' )
test_regdata_3$HA_logit_P_CV <- HA_logit_P_CV

#### B. ROC ####
roc_score_Test_3_HA <- roc(test_regdata_3$Dom, HA_logit_P_CV) #AUC score
## 0.8519 ##

## Test data 4 ##

#### A. Predicted Score ####
HA_logit_P_CV = predict(GD_1_HA_model, newdata = test_regdata_4, type = 'response' )
test_regdata_4$HA_logit_P_CV <- HA_logit_P_CV

#### B. ROC ####
roc_score_Test_4_HA <- roc(test_regdata_4$Dom, HA_logit_P_CV) #AUC score
## 0.8519 ##


######## Final model ##########
##### 3.1.2 HA + NA*NS model ##### ---  Significant
GD_3_2_2_HA_NANS_model <- glmer(Dom ~ HA + N_A*NS + (1 | vaccine_code), 
                                data=train_regdata, family=binomial(link="logit"))    

summary(GD_3_2_2_HA_NANS_model)
#### A. Predicted Score ####
HA_NANS_logit_P = predict(GD_3_2_2_HA_NANS_model, newdata = train_regdata, type = 'response' )
train_regdata$HA_NANS_logit_P <- HA_NANS_logit_P

#### B. ROC ####
roc_score_Train_HA_NANS <- roc(train_regdata$Dom, HA_NANS_logit_P) #AUC score
## 0.8246 ##

## Cross-validation using 9.HK15 ##
## Test data 1 ##

#### A. Predicted Score ####
HA_NANS_logit_P_CV = predict(GD_3_2_2_HA_NANS_model, newdata = test_regdata_1, type = 'response' )
test_regdata_1$HA_NANS_logit_P_CV <- HA_NANS_logit_P_CV

#### B. ROC ####
roc_score_Test_1_HA_NANS <- roc(test_regdata_1$Dom, HA_NANS_logit_P_CV) #AUC score
## 0.8946 ##



## Test data 2 ##

#### A. Predicted Score ####
HA_NANS_logit_P_CV = predict(GD_3_2_2_HA_NANS_model, newdata = test_regdata_2, type = 'response' )
test_regdata_2$HA_NANS_logit_P_CV <- HA_NANS_logit_P_CV

#### B. ROC ####
roc_score_Test_2_HA_NANS <- roc(test_regdata_2$Dom, HA_NANS_logit_P_CV) #AUC score
## 0.65 ##

## Test data 3 ##

#### A. Predicted Score ####
HA_NANS_logit_P_CV = predict(GD_3_2_2_HA_NANS_model, newdata = test_regdata_3, type = 'response' )
test_regdata_3$HA_NANS_logit_P_CV <- HA_NANS_logit_P_CV

#### B. ROC ####
roc_score_Test_3_HA_NANS <- roc(test_regdata_3$Dom, HA_NANS_logit_P_CV) #AUC score
## 0.7778 ##

## Test data 4 ##
#### A. Predicted Score ####
HA_NANS_logit_P_CV = predict(GD_3_2_2_HA_NANS_model, newdata = test_regdata_4, type = 'response' )
test_regdata_4$HA_NANS_logit_P_CV <- HA_NANS_logit_P_CV

#### B. ROC ####
roc_score_Test_4_HA_NANS <- roc(test_regdata_4$Dom, HA_NANS_logit_P_CV) #AUC score
## 0.8651 ##



#### Pred M1 and Other test datasets ####
Model <- c("M1_HK15","M1_HK15","M1_HK15","M1_HK15")
Dataset <- c("Train","Train","Test","Test")
Variable <- c("HA","Full","HA","Full")
AUC_ROC <- c(auc(roc_score_Train_HA), auc(roc_score_Train_HA_NANS),
             auc(roc_score_Test_1_HA), auc(roc_score_Test_1_HA_NANS))

Model_1_Pred <- data.frame(Model, Dataset, Variable, AUC_ROC)


#### M1 and Other test datasets ####
Model <- c("M1_HK15","M1_HK15","M1_HK15","M1_HK15","M1_HK15","M1_HK15",
           "M1_HK15","M1_HK15","M1_HK15","M1_HK15")
Dataset <- c("Train_M1","Train_M1","Test_M1","Test_M2","Test_M3","Test_M4",
                            "Test_M1","Test_M2","Test_M3","Test_M4")
Variable <- c("HA","Final","HA","HA","HA","HA","Final","Final","Final","Final")
AUC_ROC <- c(auc(roc_score_Train_HA), auc(roc_score_Train_HA_NANS),
             auc(roc_score_Test_1_HA), auc(roc_score_Test_2_HA),auc(roc_score_Test_3_HA),auc(roc_score_Test_4_HA),
             auc(roc_score_Test_1_HA_NANS),auc(roc_score_Test_2_HA_NANS),auc(roc_score_Test_3_HA_NANS),auc(roc_score_Test_4_HA_NANS))

Model_1_Pred_Others <- data.frame(Model, Dataset, Variable, AUC_ROC)

Model_1_Pred_Others$Dataset_F <- factor(Model_1_Pred_Others$Dataset, levels = c("Train_M1","Test_M1","Test_M2","Test_M3","Test_M4"), 
                                        labels = c("Train by M1","Test for M1","Test for M2","Test for M3","Test for M4"))

Model_1_Pred_Others$Variable_F <- factor(Model_1_Pred_Others$Variable, levels = c("HA","Final"))

Model2col4 <- c("#A1ABD7","#515b87")

jpeg(filename = "01_GenDisFlu/Fig/PredModel/FivedataCompariing.jpeg", width = 220, height = 110, units = "mm", res = 720)

Model_1_Pred_Others %>% 
  ggplot(aes(x=Dataset_F)) +
  geom_bar( aes(y=AUC_ROC, fill = Variable_F), stat='identity', width=.5, position = "dodge") + 
#  facet_grid(Variable ~ .) + 
  scale_fill_manual(values = Model2col4, name="Model") + 
  ylab("Area under curve of ROC") +
  xlab("Datasets") +
  ylim(0,1) +
  theme_bw() 

dev.off()


############# Model 2 -  Dataset 2 Kan 17 ###################

## 2. Kan17 Prediction ##
### Subset data - Training ##
train_regdata <- regdata %>% 
  filter(vaccine_code != "Dar21" & vaccine_code != "Cam20" & vaccine_code != "HK19" 
         & vaccine_code != "Kan17" & vaccine_code != "Mos99")

### Subset data - Testing ##
test_regdata <- regdata %>% 
  filter(vaccine_code == "Kan17")


## set the same intercept as Swtz13
test_regdata$vaccine_code <- "HK15"


## Cross-validation using 10.Kan17 ##
## Test data ##

######## BAse model ##########
##### 1.HA only model ##### ---  Significant
GD_1_HA_model <- glmer(Dom ~ HA + (1 | vaccine_code), 
                       data=train_regdata, family=binomial(link="logit"))    

summary(GD_1_HA_model)
#### A. Predicted Score ####
HA_logit_P = predict(GD_1_HA_model, newdata = train_regdata, type = 'response' )
train_regdata$HA_logit_P <- HA_logit_P

#### B. ROC ####
roc_score_Train_HA <- roc(train_regdata$Dom, HA_logit_P) #AUC score


## Cross-validation using 9.HK15 ##
## Test data ##

#### A. Predicted Score ####
HA_logit_P_CV = predict(GD_1_HA_model, newdata = test_regdata, type = 'response' )
test_regdata$HA_logit_P_CV <- HA_logit_P_CV

#### B. ROC ####
roc_score_Test_HA <- roc(test_regdata$Dom, HA_logit_P_CV) #AUC score


##### 3.1.2 HA + NA*NS model ##### ---  Significant
GD_3_2_2_HA_NANS_model <- glmer(Dom ~ HA + N_A*NS + (1 | vaccine_code), 
                                data=train_regdata, family=binomial(link="logit"))    
summary(GD_3_2_2_HA_NANS_model)
#### A. Predicted Score ####
HA_NANS_logit_P = predict(GD_3_2_2_HA_NANS_model, newdata = train_regdata, type = 'response' )
train_regdata$HA_NANS_logit_P <- HA_NANS_logit_P

#### B. ROC ####
roc_score_Train_HA_NANS <- roc(train_regdata$Dom, HA_NANS_logit_P) #AUC score


## Cross-validation using 9.HK15 ##
## Test data ##

#### A. Predicted Score ####
HA_NANS_logit_P_CV = predict(GD_3_2_2_HA_NANS_model, newdata = test_regdata, type = 'response' )
test_regdata$HA_NANS_logit_P_CV <- HA_NANS_logit_P_CV

#### B. ROC ####
roc_score_Test_HA_NANS <- roc(test_regdata$Dom, HA_NANS_logit_P_CV) #AUC score

Model <- c("M2_Kan17","M2_Kan17","M2_Kan17","M2_Kan17")
Dataset <- c("Train","Train","Test","Test")
Variable <- c("HA","Full","HA","Full")
AUC_ROC <- c(auc(roc_score_Train_HA), auc(roc_score_Train_HA_NANS),
             auc(roc_score_Test_HA), auc(roc_score_Test_HA_NANS))

Model_2_Pred <- data.frame(Model, Dataset, Variable, AUC_ROC)


############# Model 3 -  Dataset 3 HK19###################


## 3. HK19 Prediction ##
### Subset data - Training ##
train_regdata <- regdata %>% 
  filter(vaccine_code != "Dar21" & vaccine_code != "Cam20" & vaccine_code != "HK19" 
         & vaccine_code != "Fuj02" & vaccine_code != "Mos99")

### Subset data - Testing ##
test_regdata <- regdata %>% 
  filter(vaccine_code == "HK19")


## set the same intercept as Swtz13
test_regdata$vaccine_code <- "Kan17"

######## BAse model ##########
##### 1.HA only model ##### ---  Significant
GD_1_HA_model <- glmer(Dom ~ HA + (1 | vaccine_code), 
                       data=train_regdata, family=binomial(link="logit"))    
summary(GD_1_HA_model)
#### A. Predicted Score ####
HA_logit_P = predict(GD_1_HA_model, newdata = train_regdata, type = 'response' )
train_regdata$HA_logit_P <- HA_logit_P

#### B. ROC ####
roc_score_Train_HA <- roc(train_regdata$Dom, HA_logit_P) #AUC score

## Test data ##

#### A. Predicted Score ####
HA_logit_P_CV = predict(GD_1_HA_model, newdata = test_regdata, type = 'response' )
test_regdata$HA_logit_P_CV <- HA_logit_P_CV


######## Final model ##########
##### 3.1 NA*NS  + PA model ##### ---  Significant
GD_3_2_NANS_HA_model <- glmer(Dom ~ HA + N_A*NS + (1 | vaccine_code), 
                              data=train_regdata, family=binomial(link="logit"))    
summary(GD_3_2_NANS_HA_model)
#### A. Predicted Score ####
NANS_HA_logit_P = predict(GD_3_2_NANS_HA_model, newdata = train_regdata, type = 'response' )
train_regdata$NANS_HA_logit_P <- NANS_HA_logit_P

#### B. ROC ####
roc_score_Train_NANS_HA <- roc(train_regdata$Dom, NANS_HA_logit_P) #AUC score


## Test data ##

#### A. Predicted Score ####
NANS_HA_logit_P_CV = predict(GD_3_2_NANS_HA_model, newdata = test_regdata, type = 'response' )
test_regdata$NANS_HA_logit_P_CV <- NANS_HA_logit_P_CV

#### B. ROC ####
roc_score_Test_NANS_HA <- roc(test_regdata$Dom, NANS_HA_logit_P_CV) #AUC score


Model <- c("M3_HK19","M3_HK19","M3_HK19","M3_HK19")
Dataset <- c("Train","Train","Test","Test")
Variable <- c("HA","Full","HA","Full")
AUC_ROC <- c(auc(roc_score_Train_HA), auc(roc_score_Train_NANS_HA),
             auc(roc_score_Test_HA), auc(roc_score_Test_NANS_HA))

Model_3_Pred <- data.frame(Model, Dataset, Variable, AUC_ROC)



############# Model 4 -  Dataset 4 Cam20 ###################

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


######## BAse model ##########
##### 1.HA only model ##### ---  Significant
GD_1_HA_model <- glmer(Dom ~ HA + (1 | vaccine_code), 
                       data=train_regdata, family=binomial(link="logit"))    
summary(GD_1_HA_model)
#### A. Predicted Score ####
HA_logit_P = predict(GD_1_HA_model, newdata = train_regdata, type = 'response' )
train_regdata$HA_logit_P <- HA_logit_P

#### B. ROC ####
roc_score_Train_HA <- roc(train_regdata$Dom, HA_logit_P) #AUC score


## Cross-validation using 11.Cam20 ##
## Test data ##

#### A. Predicted Score ####
HA_logit_P_CV = predict(GD_1_HA_model, newdata = test_regdata, type = 'response' )
test_regdata$HA_logit_P_CV <- HA_logit_P_CV

#### B. ROC ####
roc_score_Test_HA <- roc(test_regdata$Dom, HA_logit_P_CV) #AUC score




######## Final model ##########
##### 3.1.1  NA*NS  + PA model ##### ---  Significant
GD_3_2_NANS_HA_model <- glmer(Dom ~ N_A * NS + HA  +(1 | vaccine_code), 
                                data=train_regdata, family=binomial(link="logit"))    
summary(GD_3_2_NANS_HA_model)
#### A. Predicted Score ####
NANS_HA_logit_P = predict(GD_3_2_NANS_HA_model, newdata = train_regdata, type = 'response' )
train_regdata$NANS_HA_logit_P <- NANS_HA_logit_P

#### B. ROC ####
roc_score_Train_NANS_HA <- roc(train_regdata$Dom, NANS_HA_logit_P) #AUC score


## Cross-validation using 11. Cam 20 ##
## Test data ##

#### A. Predicted Score ####
NANS_HA_logit_P_CV = predict(GD_3_2_NANS_HA_model, newdata = test_regdata, type = 'response' )
test_regdata$NANS_HA_logit_P_CV <- NANS_HA_logit_P_CV

#### B. ROC ####
roc_score_Test_NANS_HA <- roc(test_regdata$Dom, NANS_HA_logit_P_CV) #AUC score





Model <- c("M4_Cam20","M4_Cam20","M4_Cam20","M4_Cam20")
Dataset <- c("Train","Train","Test","Test")
Variable <- c("HA","Full","HA","Full")
AUC_ROC <- c(auc(roc_score_Train_HA), auc(roc_score_Train_NANS_HA),
             auc(roc_score_Test_HA), auc(roc_score_Test_NANS_HA))

Model_4_Pred <- data.frame(Model, Dataset, Variable, AUC_ROC)


Model_1234_Pred <- rbind(Model_1_Pred, Model_2_Pred, Model_3_Pred, Model_4_Pred)

Model_1234_Pred$Dataset <- factor(Model_1234_Pred$Dataset, levels = c("Train", "Test")) 
Model_1234_Pred$Variable <- factor(Model_1234_Pred$Variable, levels = c("HA", "Full")) 
Model_1234_Pred$Model <- factor(Model_1234_Pred$Model, levels = c("M1_HK15","M2_Kan17","M3_HK19","M4_Cam20"),
                                labels = c("Model 1 (HK15)","Model 2 (Kan17)","Model 3 (HK19)","Model 4 (Cam20)")) 


Model2col4 <- c("#80CBC4","#00695C")


jpeg(filename = "01_GenDisFlu/Fig/PredModel/FourmodelCompariing_HA_NANS.jpeg", width = 220, height = 110, units = "mm", res = 720)

Model_1234_Pred %>% 
  ggplot(aes(x=Dataset)) +
  geom_bar( aes(y=AUC_ROC, fill = Variable), stat='identity', width=.5, position = "dodge") + 
  facet_grid( ~ Model) + 
  scale_fill_manual(values = Model2col4, name="Model") + 
  ylab("Area under curve of ROC") +
  ylim(0,1) +
  theme_bw() 

dev.off()
