## Package ##
library(tidyverse)
library(lme4)
library(MuMIn)
library(pROC)
library(lmtest)
library(corrplot)
library(ggplot2)

## Set working directory ##
setwd("~/H3N2/")

save_eps_base <- function(file, width_mm, height_mm, expr) {
  grDevices::postscript(
    file,
    width  = width_mm/25.4,      # mm -> inch
    height = height_mm/25.4,
    horizontal = FALSE, onefile = FALSE, paper = "special"
  )
  on.exit(dev.off(), add = TRUE)
  force(expr)  # expr 안에서 plot.roc 같은 base plotting 실행
}

####### Evaluate the 601 Data #######

## Data import ##
Meta <- read.csv("1_Data/10_Final_MetaData_Submission/H3N2_Traindata_600strains_1_Metadata.csv", header = T, na.strings = "")
NonSyn <- read.csv("1_Data/10_Final_MetaData_Submission/H3N2_Traindata_600strains_2_NonsynDist.csv", header = T, na.strings = "")
Syn <- read.csv("1_Data/10_Final_MetaData_Submission/H3N2_Traindata_600strains_3_SynDist.csv", header = T, na.strings = "")
Others <- read.csv("1_Data/10_Final_MetaData_Submission/H3N2_Traindata_600strains_4_OtherPredictors.csv", header = T, na.strings = "")

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

Base600 <- roc(All600_Train$Dom, Base_model_1_logit_P) ## 0.801 ##
Final600 <- roc(All600_Train$Dom, Final_model_2_logit_P) ## 0.8287 ##
Control1_600 <- roc(All600_Train$Dom, Control1_model_3_logit_P) ## 0.7009 ##
Control2_600 <- roc(All600_Train$Dom, Control2_model_4_logit_P) ## 0.7598 ##


#### A. Predicted Score ####
Base_model_1_logit_P = predict(Base_model_1, newdata = train_regdata_1, type = 'response' )
train_regdata_1$Base_model_1_logit_P <- Base_model_1_logit_P


#### B. ROC ####
roc_score_Train_BaseModel_1 <- roc(train_regdata_1$Dom, Base_model_1_logit_P) #AUC score
## 0.7649 ##

roc_base_train <- roc(train_regdata_1$Dom,
                      train_regdata_1$Base_model_1_logit_P,
                      quiet = TRUE, na.rm = TRUE)


save_eps_base("5_SFigure/Sfig3_ROC/Sfig_3_E_Train1_Base_ROC.eps", 700, 300, {
  plot.roc(Control1_600,
           print.auc = TRUE, max.auc.polygon = TRUE, print.auc.cex = 1.5,
           print.thres = TRUE, print.thres.pch = 21, print.thres.col = "red",
           auc.polygon = TRUE, auc.polygon.col = "#f8ce1c",
           main = "The Control model 1 (RBD) with train data")
})


save_eps_base("5_SFigure/Sfig3_ROC/Sfig_3_F_Test1_Base_ROC.eps", 700, 300, {
  plot.roc(Control2_600,
           print.auc = TRUE, max.auc.polygon = TRUE, print.auc.cex = 1.5,
           print.thres = TRUE, print.thres.pch = 21, print.thres.col = "red",
           auc.polygon = TRUE, auc.polygon.col = "#f8ce1c",
           main = "The Control model 2 (Koel) with train data")
})

save_eps_base("5_SFigure/Sfig3_ROC/Sfig_3_A_Ctrl1_ROC.eps", 700, 300, {
  plot.roc(Base600,
           print.auc = TRUE, max.auc.polygon = TRUE, print.auc.cex = 1.5,
           print.thres = TRUE, print.thres.pch = 21, print.thres.col = "red",
           auc.polygon = TRUE, auc.polygon.col = "#f8ce1c",
           main = "The Base model with train data")
})


save_eps_base("5_SFigure/Sfig3_ROC/Sfig_3_B_Ctrl1_ROC.eps", 700, 300, {
  plot.roc(Final600,
           print.auc = TRUE, max.auc.polygon = TRUE, print.auc.cex = 1.5,
           print.thres = TRUE, print.thres.pch = 21, print.thres.col = "red",
           auc.polygon = TRUE, auc.polygon.col = "#f8ce1c",
           main = "The Final model with train data")
})



####### Evaluate the 724 Data #######

## Data import ##
Meta <- read.csv("1_Data/10_Final_MetaData_Submission/H3N2_Testdata_724strains_1_Metadata.csv", header = T, na.strings = "")
NonSyn <- read.csv("1_Data/10_Final_MetaData_Submission/H3N2_Testdata_724strains_2_NonsynDist.csv", header = T, na.strings = "")
Syn <- read.csv("1_Data/10_Final_MetaData_Submission/H3N2_Testdata_724strains_3_SynDist.csv", header = T, na.strings = "")

## Combine Data ##
All725 <- left_join(Meta, NonSyn[,c(-2,-3,-4)], by = "ID")
All725 <- left_join(All725, Syn[,c(-2,-3)], by = "ID")

### Data Subseting ###

############
## Data 1 ##
############

### Subset data - Training ##
train_regdata_1 <- All725 %>% 
  filter(vaccine_code != "Dar21" & vaccine_code != "Cam20" & vaccine_code != "HK19" 
         & vaccine_code != "Kan17" & vaccine_code != "HK15")

### Subset data - Testing ##
test_regdata_1 <- All725 %>% 
  filter(vaccine_code == "HK15")

test_regdata_1$vaccine_code <- "Swtz13"

############
## Data 2##
############

### Subset data - Training ##
train_regdata_2 <- All725 %>% 
  filter(vaccine_code != "Dar21" & vaccine_code != "Cam20" & vaccine_code != "HK19" 
         & vaccine_code != "Kan17" & vaccine_code != "Mos99")

### Subset data - Testing ##
test_regdata_2 <- All725 %>% 
  filter(vaccine_code == "Kan17")

test_regdata_2$vaccine_code <- "HK15"

############
## Data 3 ##
############

### Subset data - Training ##
train_regdata_3 <- All725 %>% 
  filter(vaccine_code != "Dar21" & vaccine_code != "Cam20" & vaccine_code != "HK19" 
         & vaccine_code != "Fuj02" & vaccine_code != "Mos99")

### Subset data - Testing ##
test_regdata_3 <- All725 %>% 
  filter(vaccine_code == "HK19")

test_regdata_3$vaccine_code <- "Kan17"

############
## Data 4 ##
############

### Subset data - Training ##
train_regdata_4 <- All725 %>% 
  filter(vaccine_code != "Dar21" & vaccine_code != "Cam20" & vaccine_code != "Cal04" 
         & vaccine_code != "Fuj02" & vaccine_code != "Mos99")

### Subset data - Testing ##
test_regdata_4 <- All725 %>% 
  filter(vaccine_code == "Cam20")

test_regdata_4$vaccine_code <- "HK19"


### Figures ###

## Evaluate the 725 Data -- Data 1 ###
## 1. Model 1 ##
Base_model_1 <- glmer(Dom ~ HA_Nonsyn + (1 | vaccine_code), 
                      data = train_regdata_1, family = binomial)

Final_model_1 <- glmer(Dom ~ HA_Nonsyn + NA_Nonsyn * NS_Nonsyn + (1 | vaccine_code), 
                       data = train_regdata_1, family = binomial)

## E, F. model 1 base - train, test
#### A. Predicted Score ####
Base_model_1_logit_P = predict(Base_model_1, newdata = train_regdata_1, type = 'response' )
train_regdata_1$Base_model_1_logit_P <- Base_model_1_logit_P


#### B. ROC ####
roc_score_Train_BaseModel_1 <- roc(train_regdata_1$Dom, Base_model_1_logit_P) #AUC score
## 0.7649 ##

roc_base_train <- roc(train_regdata_1$Dom,
                      train_regdata_1$Base_model_1_logit_P,
                      quiet = TRUE, na.rm = TRUE)

save_eps_base("5_SFigure/Sfig3_ROC/Sfig_3_E_Train1_Base_ROC.eps", 700, 300, {
  plot.roc(roc_base_train,
           print.auc = TRUE, max.auc.polygon = TRUE, print.auc.cex = 1.5,
           print.thres = TRUE, print.thres.pch = 21, print.thres.col = "red",
           auc.polygon = TRUE, auc.polygon.col = "#f8ce1c",
           main = "The Base model 1 with train data")
})

#### A. Predicted Score -- CV ####
Base_model_1_logit_P_CV = predict(Base_model_1, newdata = test_regdata_1, type = 'response' )
test_regdata_1$Base_model_1_logit_P_CV <- Base_model_1_logit_P_CV

#### B. ROC ####
roc_score_Test_Basemodel_1_CV <- roc(test_regdata_1$Dom, Base_model_1_logit_P_CV) #AUC score
## 0.8571 ##

roc_base_test  <- roc(test_regdata_1$Dom,
                      test_regdata_1$Base_model_1_logit_P_CV,
                      quiet = TRUE, na.rm = TRUE)

save_eps_base("5_SFigure/Sfig3_ROC/Sfig_3_F_Test1_Base_ROC.eps", 700, 300, {
  plot.roc(roc_base_test,
           print.auc = TRUE, max.auc.polygon = TRUE, print.auc.cex = 1.5,
           print.thres = TRUE, print.thres.pch = 21, print.thres.col = "red",
           auc.polygon = TRUE, auc.polygon.col = "#f8ce1c",
           main = "The Base model 1 with test data")
})




## G, H. model 1 final - train, test
#### A. Predicted Score ####
Final_model_1_logit_P = predict(Final_model_1, newdata = train_regdata_1, type = 'response' )
train_regdata_1$Final_model_1_logit_P <- Final_model_1_logit_P

#### B. ROC ####
roc_score_Train_Finalmodel_1 <- roc(train_regdata_1$Dom, Final_model_1_logit_P) #AUC score
## 0.8246 ##

roc_final_train <- roc(train_regdata_1$Dom,
                       train_regdata_1$Final_model_1_logit_P,
                       quiet = TRUE, na.rm = TRUE)


save_eps_base("5_SFigure/Sfig3_ROC/Sfig_3_G_Train1_Final_ROC.eps", 700, 300, {
  plot.roc(roc_final_train,
           print.auc = TRUE, max.auc.polygon = TRUE, print.auc.cex = 1.5,
           print.thres = TRUE, print.thres.pch = 21, print.thres.col = "red",
           auc.polygon = TRUE, auc.polygon.col = "#f8ce1c",
           main = "The Final model 1 with train data")
})

#### A. Predicted Score -- CV ####
Final_model_1_logit_P_CV = predict(Final_model_1, newdata = test_regdata_1, type = 'response' )
test_regdata_1$Final_model_1_logit_P_CV <- Final_model_1_logit_P_CV

#### B. ROC ####
roc_score_Test_Finalmodel_1_CV <- roc(test_regdata_1$Dom, Final_model_1_logit_P_CV) #AUC score
## 0.8946 ##


save_eps_base("5_SFigure/Sfig3_ROC/Sfig_3_H_Test1_Final_ROC.eps", 700, 300, {
  plot.roc(roc_score_Test_Finalmodel_1_CV,
           print.auc = TRUE, max.auc.polygon = TRUE, print.auc.cex = 1.5,
           print.thres = TRUE, print.thres.pch = 21, print.thres.col = "red",
           auc.polygon = TRUE, auc.polygon.col = "#f8ce1c",
           main = "The Final model 1 with test data")
})


## 2. Model 2 ##
## I, J. model 1 base - train, test
Base_model_2 <- glmer(Dom ~ HA_Nonsyn + (1 | vaccine_code), 
                      data = train_regdata_2, family = binomial)

Final_model_2 <- glmer(Dom ~ HA_Nonsyn + NA_Nonsyn * NS_Nonsyn + (1 | vaccine_code), 
                       data = train_regdata_2, family = binomial)

# Extract fixed effects
coef_summary <- summary(Final_model_2)$coefficients
# Convert coefficients to odds ratios
odds_ratios <- exp(coef_summary[, "Estimate"])


#### A. Predicted Score ####
Base_model_2_logit_P = predict(Base_model_2, newdata = train_regdata_2, type = 'response' )
train_regdata_2$Base_model_2_logit_P <- Base_model_2_logit_P

#### B. ROC ####
roc_score_Train_BaseModel_2 <- roc(train_regdata_2$Dom, Base_model_2_logit_P) #AUC score
## 0.7655 ##


save_eps_base("5_SFigure/Sfig3_ROC/Sfig_3_I_Train2_Base_ROC.eps", 700, 300, {
  plot.roc(roc_score_Train_BaseModel_2,
           print.auc = TRUE, max.auc.polygon = TRUE, print.auc.cex = 1.5,
           print.thres = TRUE, print.thres.pch = 21, print.thres.col = "red",
           auc.polygon = TRUE, auc.polygon.col = "#f8ce1c",
           main = "The Base model 2 with train data")
})


Base_model_2_logit_P_CV = predict(Base_model_2, newdata = test_regdata_2, type = 'response' )
test_regdata_2$Base_model_2_logit_P_CV <- Base_model_2_logit_P_CV

#### B. ROC ####
roc_score_Test_Basemodel_2_CV <- roc(test_regdata_2$Dom, Base_model_2_logit_P_CV) #AUC score
## 0.8438 ##

save_eps_base("5_SFigure/Sfig3_ROC/Sfig_3_J_Test2_Base_ROC.eps", 700, 300, {
  plot.roc(roc_score_Test_Basemodel_2_CV,
           print.auc = TRUE, max.auc.polygon = TRUE, print.auc.cex = 1.5,
           print.thres = TRUE, print.thres.pch = 21, print.thres.col = "red",
           auc.polygon = TRUE, auc.polygon.col = "#f8ce1c",
           main = "The Base model 2 with test data")
})

## K, L model 2 final - train, test
#### A. Predicted Score ####
Final_model_2_logit_P = predict(Final_model_2, newdata = train_regdata_2, type = 'response' )
train_regdata_2$Final_model_2_logit_P <- Final_model_2_logit_P

#### B. ROC ####
roc_score_Train_Finalmodel_2 <- roc(train_regdata_2$Dom, Final_model_2_logit_P) #AUC score
## 0.819 ##

save_eps_base("5_SFigure/Sfig3_ROC/Sfig_3_K_Train2_Final_ROC.eps", 700, 300, {
  plot.roc(roc_score_Train_Finalmodel_2,
           print.auc = TRUE, max.auc.polygon = TRUE, print.auc.cex = 1.5,
           print.thres = TRUE, print.thres.pch = 21, print.thres.col = "red",
           auc.polygon = TRUE, auc.polygon.col = "#f8ce1c",
           main = "The Base model 2 with train data")
})


#### A. Predicted Score -- CV ####
Final_model_2_logit_P_CV = predict(Final_model_2, newdata = test_regdata_2, type = 'response' )
test_regdata_2$Final_model_2_logit_P_CV <- Final_model_2_logit_P_CV

#### B. ROC ####
roc_score_Test_Finalmodel_2_CV <- roc(test_regdata_2$Dom, Final_model_2_logit_P_CV) #AUC score
## 0.5375 ##


save_eps_base("5_SFigure/Sfig3_ROC/Sfig_3_L_Test2_Final_ROC.eps", 700, 300, {
  plot.roc(roc_score_Test_Finalmodel_2_CV,
           print.auc = TRUE, max.auc.polygon = TRUE, print.auc.cex = 1.5,
           print.thres = TRUE, print.thres.pch = 21, print.thres.col = "red",
           auc.polygon = TRUE, auc.polygon.col = "#f8ce1c",
           main = "The Base model 2 with test data")
})



## 3. Model 3 ##
## M, N model 3 base - train, test
Base_model_3 <- glmer(Dom ~ HA_Nonsyn + (1 | vaccine_code), 
                      data = train_regdata_3, family = binomial)

Final_model_3 <- glmer(Dom ~ HA_Nonsyn + NA_Nonsyn * NS_Nonsyn + (1 | vaccine_code), 
                       data = train_regdata_3, family = binomial)

# Extract fixed effects
coef_summary <- summary(Final_model_3)$coefficients
# Convert coefficients to odds ratios
odds_ratios <- exp(coef_summary[, "Estimate"])


#### A. Predicted Score ####
Base_model_3_logit_P = predict(Base_model_3, newdata = train_regdata_3, type = 'response' )
train_regdata_3$Base_model_3_logit_P <- Base_model_3_logit_P

#### B. ROC ####
roc_score_Train_BaseModel_3 <- roc(train_regdata_3$Dom, Base_model_3_logit_P) #AUC score
## 0.7731 ##

save_eps_base("5_SFigure/Sfig3_ROC/Sfig_3_M_Train3_Base_ROC.eps", 700, 300, {
  plot.roc(roc_score_Train_BaseModel_3,
           print.auc = TRUE, max.auc.polygon = TRUE, print.auc.cex = 1.5,
           print.thres = TRUE, print.thres.pch = 21, print.thres.col = "red",
           auc.polygon = TRUE, auc.polygon.col = "#f8ce1c",
           main = "The Base model 3 with train data")
})

#### A. Predicted Score -- CV ####
Base_model_3_logit_P_CV = predict(Base_model_3, newdata = test_regdata_3, type = 'response' )
test_regdata_3$Base_model_3_logit_P_CV <- Base_model_3_logit_P_CV

#### B. ROC ####
roc_score_Test_Basemodel_3_CV <- roc(test_regdata_3$Dom, Base_model_3_logit_P_CV) #AUC score
## 0.8519 ##

save_eps_base("5_SFigure/Sfig3_ROC/Sfig_3_N_Test3_Base_ROC.eps", 700, 300, {
  plot.roc(roc_score_Test_Basemodel_3_CV,
           print.auc = TRUE, max.auc.polygon = TRUE, print.auc.cex = 1.5,
           print.thres = TRUE, print.thres.pch = 21, print.thres.col = "red",
           auc.polygon = TRUE, auc.polygon.col = "#f8ce1c",
           main = "The Base model 3 with test data")
})


## O, P model 3 final - train, test
#### A. Predicted Score ####
Final_model_3_logit_P = predict(Final_model_3, newdata = train_regdata_3, type = 'response' )
train_regdata_3$Final_model_3_logit_P <- Final_model_3_logit_P

#### B. ROC ####
roc_score_Train_Finalmodel_3 <- roc(train_regdata_3$Dom, Final_model_3_logit_P) #AUC score
## 0.8144 ##

save_eps_base("5_SFigure/Sfig3_ROC/Sfig_3_O_Train3_Final_ROC.eps", 700, 300, {
  plot.roc(roc_score_Train_Finalmodel_3,
           print.auc = TRUE, max.auc.polygon = TRUE, print.auc.cex = 1.5,
           print.thres = TRUE, print.thres.pch = 21, print.thres.col = "red",
           auc.polygon = TRUE, auc.polygon.col = "#f8ce1c",
           main = "The Final model 3 with train data")
})



#### A. Predicted Score -- CV ####
Final_model_3_logit_P_CV = predict(Final_model_3, newdata = test_regdata_3, type = 'response' )
test_regdata_3$Final_model_3_logit_P_CV <- Final_model_3_logit_P_CV

#### B. ROC ####
roc_score_Test_Finalmodel_3_CV <- roc(test_regdata_3$Dom, Final_model_3_logit_P_CV) #AUC score
## 0.6296 ##

save_eps_base("5_SFigure/Sfig3_ROC/Sfig_3_P_Test3_Final_ROC.eps", 700, 300, {
  plot.roc(roc_score_Test_Finalmodel_3_CV,
           print.auc = TRUE, max.auc.polygon = TRUE, print.auc.cex = 1.5,
           print.thres = TRUE, print.thres.pch = 21, print.thres.col = "red",
           auc.polygon = TRUE, auc.polygon.col = "#f8ce1c",
           main = "The Final model 3 with test data")
})


## 4. Model 4 ##
## Q, R model 4 base - train, test

Base_model_4 <- glmer(Dom ~ HA_Nonsyn + (1 | vaccine_code), 
                      data = train_regdata_4, family = binomial)

Final_model_4 <- glmer(Dom ~ HA_Nonsyn + NA_Nonsyn * NS_Nonsyn + (1 | vaccine_code), 
                       data = train_regdata_4, family = binomial)

# Extract fixed effects
coef_summary <- summary(Final_model_4)$coefficients
# Convert coefficients to odds ratios
odds_ratios <- exp(coef_summary[, "Estimate"])


#### A. Predicted Score ####
Base_model_4_logit_P = predict(Base_model_4, newdata = train_regdata_4, type = 'response' )
train_regdata_4$Base_model_4_logit_P <- Base_model_4_logit_P

#### B. ROC ####
roc_score_Train_BaseModel_4 <- roc(train_regdata_4$Dom, Base_model_4_logit_P) #AUC score
## 0.795 ##

save_eps_base("5_SFigure/Sfig3_ROC/Sfig_3_Q_Train4_Base_ROC.eps", 700, 300, {
  plot.roc(roc_score_Train_BaseModel_4,
           print.auc = TRUE, max.auc.polygon = TRUE, print.auc.cex = 1.5,
           print.thres = TRUE, print.thres.pch = 21, print.thres.col = "red",
           auc.polygon = TRUE, auc.polygon.col = "#f8ce1c",
           main = "The Base model 4 with train data")
})



#### A. Predicted Score -- CV ####
Base_model_4_logit_P_CV = predict(Base_model_4, newdata = test_regdata_4, type = 'response' )
test_regdata_4$Base_model_4_logit_P_CV <- Base_model_4_logit_P_CV

#### B. ROC ####
roc_score_Test_Basemodel_4_CV <- roc(test_regdata_4$Dom, Base_model_4_logit_P_CV) #AUC score
## 0.6587 ##

save_eps_base("5_SFigure/Sfig3_ROC/Sfig_3_R_Test4_Base_ROC.eps", 700, 300, {
  plot.roc(roc_score_Test_Basemodel_4_CV,
           print.auc = TRUE, max.auc.polygon = TRUE, print.auc.cex = 1.5,
           print.thres = TRUE, print.thres.pch = 21, print.thres.col = "red",
           auc.polygon = TRUE, auc.polygon.col = "#f8ce1c",
           main = "The Base model 4 with test data")
})


## S, T model 4 final - train, test
#### A. Predicted Score ####
Final_model_4_logit_P = predict(Final_model_4, newdata = train_regdata_4, type = 'response' )
train_regdata_4$Final_model_4_logit_P <- Final_model_4_logit_P

#### B. ROC ####
roc_score_Train_Finalmodel_4 <- roc(train_regdata_4$Dom, Final_model_4_logit_P) #AUC score
## 0.8187 ##

save_eps_base("5_SFigure/Sfig3_ROC/Sfig_3_S_Train4_Final_ROC.eps", 700, 300, {
  plot.roc(roc_score_Train_Finalmodel_4,
           print.auc = TRUE, max.auc.polygon = TRUE, print.auc.cex = 1.5,
           print.thres = TRUE, print.thres.pch = 21, print.thres.col = "red",
           auc.polygon = TRUE, auc.polygon.col = "#f8ce1c",
           main = "The Final model 4 with train data")
})


#### A. Predicted Score -- CV ####
Final_model_4_logit_P_CV = predict(Final_model_4, newdata = test_regdata_4, type = 'response' )
test_regdata_4$Final_model_4_logit_P_CV <- Final_model_4_logit_P_CV

#### B. ROC ####
roc_score_Test_Finalmodel_4_CV <- roc(test_regdata_4$Dom, Final_model_4_logit_P_CV) #AUC score
## 0.8095 ##

save_eps_base("5_SFigure/Sfig3_ROC/Sfig_3_T_Test4_Final_ROC.eps", 700, 300, {
  plot.roc(roc_score_Test_Finalmodel_4_CV,
           print.auc = TRUE, max.auc.polygon = TRUE, print.auc.cex = 1.5,
           print.thres = TRUE, print.thres.pch = 21, print.thres.col = "red",
           auc.polygon = TRUE, auc.polygon.col = "#f8ce1c",
           main = "The Final model 4 with test data")
})


#### 5. Summary Figure ####
Model <- c("M1","M1","M1","M1","M2","M2","M2","M2",
           "M3","M3","M3","M3","M4","M4","M4","M4")
Dataset <- c("Train","Test","Train","Test","Train","Test","Train","Test",
             "Train","Test","Train","Test","Train","Test","Train","Test")
Variable <- c("Base","Base","Final","Final","Base","Base","Final","Final",
              "Base","Base","Final","Final","Base","Base","Final","Final")
AUC_ROC <- c(auc(roc_score_Train_BaseModel_1), auc(roc_score_Test_Basemodel_1_CV),
             auc(roc_score_Train_Finalmodel_1), auc(roc_score_Test_Finalmodel_1_CV),
             auc(roc_score_Train_BaseModel_2), auc(roc_score_Test_Basemodel_2_CV),
             auc(roc_score_Train_Finalmodel_2), auc(roc_score_Test_Finalmodel_2_CV),
             auc(roc_score_Train_BaseModel_3), auc(roc_score_Test_Basemodel_3_CV),
             auc(roc_score_Train_Finalmodel_3), auc(roc_score_Test_Finalmodel_3_CV),
             auc(roc_score_Train_BaseModel_4), auc(roc_score_Test_Basemodel_4_CV),
             auc(roc_score_Train_Finalmodel_4), auc(roc_score_Test_Finalmodel_4_CV))

Model_Pred <- data.frame(Model, Dataset, Variable, AUC_ROC)
Model_Pred$Model_F <- factor(Model_Pred$Model, levels = c("M1","M2","M3","M4"), labels = c("Model 1","Model 2","Model 3","Model 4"))
Model_Pred$Dataset_F <- factor(Model_Pred$Dataset, levels = c("Train","Test"))
Model_Pred$Variable_F <- factor(Model_Pred$Variable, levels = c("Base","Final"))

Model2col4 <- c("#A1ABD7","#515b87")

Model_Pred %>% 
  ggplot(aes(x=Dataset_F)) +
  geom_bar( aes(y=AUC_ROC, fill = Variable_F), stat='identity', width=.5, position = "dodge") + 
  facet_grid(. ~ Model_F) + 
  scale_fill_manual(values = Model2col4, name="Variable") + 
  ylab("Area under curve of ROC") +
  xlab("Datasets") +
  ylim(0,1) +
  theme_bw() +  theme(strip.background = element_rect(fill="white"), 
                      strip.text.x = element_text(size = 10, family = "sans"))


dev.off()