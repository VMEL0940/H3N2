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

## Data import ##
Meta <- read.csv("Data/AllData/H3N2_Testdata_724strains_1_Metadata.csv", header = T, na.strings = "")
NonSyn <- read.csv("Data/AllData/H3N2_Testdata_724strains_2_NonsynDist.csv", header = T, na.strings = "")
Syn <- read.csv("Data/AllData/H3N2_Testdata_724strains_3_SynDist.csv", header = T, na.strings = "")

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
         & vaccine_code != "Kan17")

### Subset data - Testing ##
test_regdata_2 <- All725 %>% 
  filter(vaccine_code == "Kan17")

test_regdata_2$vaccine_code <- "HK15"

############
## Data 3 ##
############

### Subset data - Training ##
train_regdata_3 <- All725 %>% 
  filter(vaccine_code != "Dar21" & vaccine_code != "Cam20" & vaccine_code != "HK19" )

### Subset data - Testing ##
test_regdata_3 <- All725 %>% 
  filter(vaccine_code == "HK19")

test_regdata_3$vaccine_code <- "Kan17"

############
## Data 4 ##
############

### Subset data - Training ##
train_regdata_4 <- All725 %>% 
  filter(vaccine_code != "Dar21" & vaccine_code != "Cam20")

### Subset data - Testing ##
test_regdata_4 <- All725 %>% 
  filter(vaccine_code == "Cam20")

test_regdata_4$vaccine_code <- "HK19"


## Evaluate the 725 Data -- Data 1 ###
## 1. Regression model selection -Automated / Test## 
options(na.action = "na.fail")
## No automated function ##

full_model_NS_1 <- glmer(Dom ~ PB2_Nonsyn + PB1_Nonsyn + PA_Nonsyn + HA_Nonsyn + NP_Nonsyn + NA_Nonsyn + M_Nonsyn + NS_Nonsyn + (1 | vaccine_code), 
                    data = train_regdata_1, 
                    family = binomial)

model_set_NS_1  <- dredge(full_model_NS_1, fixed = ~HA_Nonsyn)

print(model_set_NS_1)

# Get the best model
best_model_NS_1 <- get.models(model_set_NS_1 , 1)[[1]]

### SyN MOdel selection ###
full_model_Syn_1 <- glmer(Dom ~ PB2_Nonsyn + PB1_Syn + PA_Syn + HA_Syn + NP_Syn + NA_Syn + M_Syn + NS_Syn + (1 | vaccine_code), 
                     data = train_regdata_1, 
                     family = binomial)

model_set_Syn_1 <- dredge(full_model_Syn_1)

print(model_set_Syn_1)

# Get the best model
best_model <- get.models(model_set_Syn_1, 1)[[1]]

## Final selection -- Combine ###
full_model_1 <- glmer(Dom ~ HA_Nonsyn + NA_Nonsyn + M_Nonsyn + NP_Nonsyn + NS_Nonsyn + PA_Syn + HA_Syn + NS_Syn + M_Syn+ (1 | vaccine_code), 
                     data = train_regdata_1, 
                     family = binomial)

model_set_1 <- dredge(full_model_1, fixed = ~HA_Nonsyn)

print(model_set_1)

# Get the best model
best_model <- get.models(model_set_1, 1)[[1]]

### Best moodel by Automated selection -- HA NS + M Syn + NA Nonsyn + NS
Base_model_1 <- glmer(Dom ~ HA_Nonsyn + (1 | vaccine_code), 
                       data = train_regdata_1, 
                       family = binomial)

Final_model_1 <- glmer(Dom ~ HA_Nonsyn + NA_Nonsyn + NS_Nonsyn + M_Syn + (1 | vaccine_code), 
                     data = train_regdata_1, 
                     family = binomial)

Final_model_1_1 <- glmer(Dom ~ HA_Nonsyn * NA_Nonsyn + NS_Nonsyn + M_Syn + (1 | vaccine_code), 
                       data = train_regdata_1, 
                       family = binomial)

Final_model_1_2 <- glmer(Dom ~ HA_Nonsyn * NS_Nonsyn + NA_Nonsyn + M_Syn + (1 | vaccine_code), 
                         data = train_regdata_1, 
                         family = binomial)

Final_model_1_3 <- glmer(Dom ~ HA_Nonsyn * M_Syn + NS_Nonsyn + NA_Nonsyn +  (1 | vaccine_code), 
                         data = train_regdata_1, 
                         family = binomial) #########

Final_model_1_4 <- glmer(Dom ~ NA_Nonsyn * NS_Nonsyn + HA_Nonsyn + M_Syn + (1 | vaccine_code), 
                         data = train_regdata_1, 
                         family = binomial)

Final_model_1_5 <- glmer(Dom ~ NA_Nonsyn * M_Syn + NS_Nonsyn + HA_Nonsyn +  (1 | vaccine_code), 
                         data = train_regdata_1, 
                         family = binomial)

Final_model_1_6 <- glmer(Dom ~ M_Syn * NS_Nonsyn + HA_Nonsyn + NA_Nonsyn + (1 | vaccine_code), 
                         data = train_regdata_1, 
                         family = binomial) #########

# Compare models using AIC
AIC(Base_model_1, Final_model_1, Final_model_1_1, Final_model_1_2, Final_model_1_3, Final_model_1_4, Final_model_1_5 , Final_model_1_6)

Final_model_1_3 <- glmer(Dom ~ HA_Nonsyn * M_Syn + NS_Nonsyn + NA_Nonsyn +  (1 | vaccine_code), 
                         data = train_regdata_1, 
                         family = binomial) #########

Final_model_1_3_1 <- glmer(Dom ~ HA_Nonsyn * M_Syn + NS_Nonsyn * NA_Nonsyn +  (1 | vaccine_code), 
                         data = train_regdata_1, 
                         family = binomial) #################

Final_model_1_6_1 <- glmer(Dom ~ M_Syn * NS_Nonsyn + HA_Nonsyn * NA_Nonsyn + (1 | vaccine_code), 
                         data = train_regdata_1, 
                         family = binomial) #########

Final_model_1_1_1 <- glmer(Dom ~ HA_Nonsyn + NS_Nonsyn * NA_Nonsyn + M_Syn + (1 | vaccine_code), 
                         data = train_regdata_1, 
                         family = binomial) #########

Final_model_1_1_2 <- glmer(Dom ~ HA_Nonsyn + NS_Nonsyn * NA_Nonsyn + (1 | vaccine_code), 
                           data = train_regdata_1, 
                           family = binomial) #########

AIC(Base_model_1, Final_model_1, Final_model_1_1, Final_model_1_3, Final_model_1_3_1, Final_model_1_6, Final_model_1_6_1, 
    Final_model_1_1_1, Final_model_1_1_2)

#### A. Predicted Score ####
Final_model_1_logit_P = predict(Base_model_1, newdata = train_regdata_1, type = 'response' )
train_regdata_1$Final_model_1_logit_P <- Final_model_1_logit_P

#### B. ROC ####
roc_score_Train_Finalmodel_1 <- roc(train_regdata_1$Dom, Final_model_1_logit_P) #AUC score
## 0.7649 ##

#### A. Predicted Score -- CV ####
Final_model_1_logit_P_CV = predict(Base_model_1, newdata = test_regdata_1, type = 'response' )
test_regdata_1$Final_model_1_logit_P_CV <- Final_model_1_logit_P_CV

#### B. ROC ####
roc_score_Test_Finalmodel_1_CV <- roc(test_regdata_1$Dom, Final_model_1_logit_P_CV) #AUC score
## 0.8571 ##

#### A. Predicted Score ####
Final_model_1_logit_P = predict(Final_model_1_1_2, newdata = train_regdata_1, type = 'response' )
train_regdata_1$Final_model_1_logit_P <- Final_model_1_logit_P

#### B. ROC ####
roc_score_Train_Finalmodel_1 <- roc(train_regdata_1$Dom, Final_model_1_logit_P) #AUC score
## 0.8119 ##

#### A. Predicted Score -- CV ####
Final_model_1_logit_P_CV = predict(Final_model_1_1_2, newdata = test_regdata_1, type = 'response' )
test_regdata_1$Final_model_1_logit_P_CV <- Final_model_1_logit_P_CV

#### B. ROC ####
roc_score_Test_Finalmodel_1_CV <- roc(test_regdata_1$Dom, Final_model_1_logit_P_CV) #AUC score
## 0.7185 ##







## Evaluate the 725 Data -- Data 2 ###
## 1. Regression model selection -Automated / Test## 
full_model_NS_2 <- glmer(Dom ~ PB2_Nonsyn + PB1_Nonsyn + PA_Nonsyn + HA_Nonsyn + NP_Nonsyn + NA_Nonsyn + M_Nonsyn + NS_Nonsyn + (1 | vaccine_code), 
                         data = train_regdata_2, 
                         family = binomial)

model_set_NS_2  <- dredge(full_model_NS_2)

print(model_set_NS_2 )

# Get the best model
best_model_NS_2 <- get.models(model_set_NS_2 , 1)[[1]]

### SyN MOdel selection ###
full_model_Syn_2 <- glmer(Dom ~ PB2_Nonsyn + PB1_Syn + PA_Syn + HA_Syn + NP_Syn + NA_Syn + M_Syn + NS_Syn + (1 | vaccine_code), 
                          data = train_regdata_2, 
                          family = binomial)

model_set_Syn_2 <- dredge(full_model_Syn_2)

print(model_set_Syn_2)

# Get the best model
best_model_2 <- get.models(model_set_Syn_2, 1)[[1]]

## Final selection -- Combine ###
full_model_2 <- glmer(Dom ~ HA_Nonsyn + NA_Nonsyn + PA_Nonsyn + NS_Nonsyn + PB1_Syn + HA_Syn + NA_Syn + NS_Syn+ (1 | vaccine_code), 
                      data = train_regdata_2, 
                      family = binomial)

model_set_2 <- dredge(full_model_2)

print(model_set_2)

# Get the best model
best_model_2 <- get.models(model_set_2, 1)[[1]]

### Best moodel by Automated selection -- HA NS + M Syn + NA Nonsyn + NS
Base_model_2 <- glmer(Dom ~ HA_Nonsyn + (1 | vaccine_code), 
                      data = train_regdata_2, 
                      family = binomial)

Final_model_2 <- glmer(Dom ~ HA_Nonsyn + NA_Nonsyn + NS_Nonsyn + PA_Nonsyn  +(1 | vaccine_code), 
                       data = train_regdata_2, 
                       family = binomial)

Final_model_2_1 <- glmer(Dom ~ HA_Nonsyn * NA_Nonsyn + NS_Nonsyn + PA_Nonsyn + (1 | vaccine_code), 
                         data = train_regdata_2, 
                         family = binomial) #########

Final_model_2_2 <- glmer(Dom ~ HA_Nonsyn * NS_Nonsyn + NA_Nonsyn + PA_Nonsyn + (1 | vaccine_code), 
                         data = train_regdata_2, 
                         family = binomial)

Final_model_2_3 <- glmer(Dom ~ HA_Nonsyn * PA_Nonsyn + NS_Nonsyn + NA_Nonsyn +  (1 | vaccine_code), 
                         data = train_regdata_2, 
                         family = binomial) 

Final_model_2_4 <- glmer(Dom ~ NA_Nonsyn * NS_Nonsyn + HA_Nonsyn + PA_Nonsyn + (1 | vaccine_code), 
                         data = train_regdata_2, 
                         family = binomial) #########

Final_model_2_5 <- glmer(Dom ~ NA_Nonsyn * PA_Nonsyn + NS_Nonsyn + HA_Nonsyn +  (1 | vaccine_code), 
                         data = train_regdata_2, 
                         family = binomial)

Final_model_2_6 <- glmer(Dom ~ PA_Nonsyn * NS_Nonsyn + HA_Nonsyn + NA_Nonsyn + (1 | vaccine_code), 
                         data = train_regdata_2, 
                         family = binomial) 

# Compare models using AIC
AIC(Base_model_2, Final_model_2, Final_model_2_1, Final_model_2_2, Final_model_2_3, Final_model_2_4, Final_model_2_5 , Final_model_2_6)

Final_model_2_4 <- glmer(Dom ~ NA_Nonsyn * NS_Nonsyn + HA_Nonsyn + PA_Nonsyn + (1 | vaccine_code), 
                         data = train_regdata_2, 
                         family = binomial) ##################

Final_model_2_4_1 <- glmer(Dom ~ HA_Nonsyn * PA_Nonsyn + NA_Nonsyn * NS_Nonsyn  + (1 | vaccine_code), 
                           data = train_regdata_2, 
                           family = binomial) #########

AIC(Base_model_2, Final_model_2, Final_model_2_4, Final_model_2_4_1)


#### A. Predicted Score ####
Base_model_2_logit_P = predict(Base_model_2, newdata = train_regdata_2, type = 'response' )
train_regdata_2$Base_model_2_logit_P <- Base_model_2_logit_P

#### B. ROC ####
roc_score_Train_Basemodel_2 <- roc(train_regdata_2$Dom, Base_model_2_logit_P) #AUC score
## 0.7772 ##

#### A. Predicted Score -- CV ####
Base_model_2_logit_P_CV = predict(Base_model_2, newdata = test_regdata_2, type = 'response' )
test_regdata_2$Base_model_2_logit_P_CV <- Base_model_2_logit_P_CV

#### B. ROC ####
roc_score_Test_Basemodel_2_CV <- roc(test_regdata_2$Dom, Base_model_2_logit_P_CV) #AUC score
## 0.8438 ##

#### A. Predicted Score ####
Final_model_2_logit_P = predict(Final_model_2_4, newdata = train_regdata_2, type = 'response' )
train_regdata_2$Final_model_2_logit_P <- Final_model_2_logit_P

#### B. ROC ####
roc_score_Train_Finalmodel_2 <- roc(train_regdata_2$Dom, Final_model_2_logit_P) #AUC score
##  0.8439 ##

#### A. Predicted Score -- CV ####
Final_model_2_logit_P_CV = predict(Final_model_2_4, newdata = test_regdata_2, type = 'response' )
test_regdata_2$Final_model_2_logit_P_CV <- Final_model_2_logit_P_CV

#### B. ROC ####
roc_score_Test_Finalmodel_2_CV <- roc(test_regdata_2$Dom, Final_model_2_logit_P_CV) #AUC score
## 0.5 ##



## Evaluate the 725 Data -- Data 3 ###
## 1. Regression model selection -Automated / Test## 
options(na.action = "na.fail")
## No automated function ##
full_model_NS_3 <- glmer(Dom ~ PB2_Nonsyn + PB1_Nonsyn + PA_Nonsyn + HA_Nonsyn + NP_Nonsyn + NA_Nonsyn + M_Nonsyn + NS_Nonsyn + (1 | vaccine_code), 
                         data = train_regdata_3, 
                         family = binomial)

model_set_NS_3  <- dredge(full_model_NS_3)

print(model_set_NS_3)

# Get the best model
best_model_NS_3 <- get.models(model_set_NS_3, 1)[[1]]

### SyN MOdel selection ###
full_model_Syn_3 <- glmer(Dom ~ PB2_Nonsyn + PB1_Syn + PA_Syn + HA_Syn + NP_Syn + NA_Syn + M_Syn + NS_Syn + (1 | vaccine_code), 
                          data = train_regdata_3, 
                          family = binomial)

model_set_Syn_3 <- dredge(full_model_Syn_3)

print(model_set_Syn_3)

# Get the best model
best_model_3 <- get.models(model_set_Syn_3, 1)[[1]]

## Final selection -- Combine ###
full_model_3 <- glmer(Dom ~ HA_Nonsyn + NA_Nonsyn + M_Nonsyn + PA_Nonsyn + HA_Syn + PB1_Syn + PB2_Syn + (1 | vaccine_code), 
                      data = train_regdata_3, 
                      family = binomial)

model_set_3 <- dredge(full_model_3)

print(model_set_3)

# Get the best model
best_model_3 <- get.models(model_set_3, 1)[[1]]

### Best moodel by Automated selection -- HA NS + M Syn + NA Nonsyn + NS
Base_model_3 <- glmer(Dom ~ HA_Nonsyn + (1 | vaccine_code), 
                      data = train_regdata_3, 
                      family = binomial)

Final_model_3 <- glmer(Dom ~ HA_Nonsyn + NA_Nonsyn + M_Nonsyn + PA_Nonsyn  +(1 | vaccine_code), 
                       data = train_regdata_3, 
                       family = binomial)

Final_model_3_1 <- glmer(Dom ~ HA_Nonsyn * NA_Nonsyn + NS_Nonsyn + PA_Nonsyn + (1 | vaccine_code), 
                         data = train_regdata_3, 
                         family = binomial) 

Final_model_3_2 <- glmer(Dom ~ HA_Nonsyn * NS_Nonsyn + NA_Nonsyn + PA_Nonsyn + (1 | vaccine_code), 
                         data = train_regdata_3, 
                         family = binomial)

Final_model_3_3 <- glmer(Dom ~ HA_Nonsyn * PA_Nonsyn + NS_Nonsyn + NA_Nonsyn +  (1 | vaccine_code), 
                         data = train_regdata_3, 
                         family = binomial) 

Final_model_3_4 <- glmer(Dom ~ NA_Nonsyn * NS_Nonsyn + HA_Nonsyn + PA_Nonsyn + (1 | vaccine_code), 
                         data = train_regdata_3, 
                         family = binomial) #########

Final_model_3_5 <- glmer(Dom ~ NA_Nonsyn * PA_Nonsyn + NS_Nonsyn + HA_Nonsyn +  (1 | vaccine_code), 
                         data = train_regdata_3, 
                         family = binomial)

Final_model_3_6 <- glmer(Dom ~ PA_Nonsyn * NS_Nonsyn + HA_Nonsyn + NA_Nonsyn + (1 | vaccine_code), 
                         data = train_regdata_3, 
                         family = binomial) 


# Compare models using AIC
AIC(Base_model_3, Final_model_3, Final_model_3_1, Final_model_3_2, Final_model_3_3, Final_model_3_4, Final_model_3_5, Final_model_3_6)

#### A. Predicted Score ####
Base_model_3_logit_P = predict(Base_model_3, newdata = train_regdata_3, type = 'response' )
train_regdata_3$Base_model_3_logit_P <- Base_model_3_logit_P

#### B. ROC ####
roc_score_Train_Basemodel_3 <- roc(train_regdata_3$Dom, Base_model_3_logit_P) #AUC score
## 0.7768 ##

#### A. Predicted Score -- CV ####
Base_model_3_logit_P_CV = predict(Base_model_3, newdata = test_regdata_3, type = 'response' )
test_regdata_3$Base_model_3_logit_P_CV <- Base_model_3_logit_P_CV

#### B. ROC ####
roc_score_Test_Basemodel_3_CV <- roc(test_regdata_3$Dom, Base_model_3_logit_P_CV) #AUC score
## 0.8519 ##

#### A. Predicted Score ####
Final_model_3_logit_P = predict(Final_model_3_1, newdata = train_regdata_3, type = 'response' )
train_regdata_3$Final_model_3_logit_P <- Final_model_3_logit_P

#### B. ROC ####
roc_score_Train_Finalmodel_3 <- roc(train_regdata_3$Dom, Final_model_3_logit_P) #AUC score
## 0.8107 ##

#### A. Predicted Score -- CV ####
Final_model_3_logit_P_CV = predict(Final_model_3_1, newdata = test_regdata_3, type = 'response' )
test_regdata_3$Final_model_3_logit_P_CV <- Final_model_3_logit_P_CV

#### B. ROC ####
roc_score_Test_Finalmodel_3_CV <- roc(test_regdata_3$Dom, Final_model_3_logit_P_CV) #AUC score
## 0.7926 ##



## Evaluate the 725 Data -- Data 4 ###
## 1. Regression model selection -Automated / Test## 
## No automated function ##
full_model_NS_4 <- glmer(Dom ~ PB2_Nonsyn + PB1_Nonsyn + PA_Nonsyn + HA_Nonsyn + NP_Nonsyn + NA_Nonsyn + M_Nonsyn + NS_Nonsyn + (1 | vaccine_code), 
                         data = train_regdata_4, 
                         family = binomial)

model_set_NS_4  <- dredge(full_model_NS_4)

print(model_set_NS_4)

# Get the best model
best_model_NS_4 <- get.models(model_set_NS_4, 1)[[1]]

### SyN MOdel selection ###
full_model_Syn_4 <- glmer(Dom ~ PB2_Nonsyn + PB1_Syn + PA_Syn + HA_Syn + NP_Syn + NA_Syn + M_Syn + NS_Syn + (1 | vaccine_code), 
                          data = train_regdata_4, 
                          family = binomial)

model_set_Syn_4 <- dredge(full_model_Syn_4)

print(model_set_Syn_4)

# Get the best model
best_model_4 <- get.models(model_set_Syn_4, 1)[[1]]

## Final selection -- Combine ###
full_model_4 <- glmer(Dom ~ HA_Nonsyn + NA_Nonsyn + NA_Nonsyn + PA_Nonsyn + HA_Syn + M_Syn + PB1_Syn + (1 | vaccine_code), 
                      data = train_regdata_4, 
                      family = binomial)

model_set_4 <- dredge(full_model_4)

print(model_set_4)

# Get the best model
best_model_4 <- get.models(model_set_4, 1)[[1]]

### Best moodel by Automated selection -- HA NS + M Syn + NA Nonsyn + NS
Base_model_4 <- glmer(Dom ~ HA_Nonsyn + (1 | vaccine_code), 
                      data = train_regdata_4, 
                      family = binomial)

Final_model_4 <- glmer(Dom ~ HA_Nonsyn + NA_Nonsyn + PA_Nonsyn + M_Syn +(1 | vaccine_code), 
                       data = train_regdata_4, 
                       family = binomial)

Final_model_4_1 <- glmer(Dom ~ HA_Nonsyn * NA_Nonsyn + PA_Nonsyn + M_Syn + (1 | vaccine_code), 
                         data = train_regdata_4, 
                         family = binomial) #########

Final_model_4_2 <- glmer(Dom ~ HA_Nonsyn * PA_Nonsyn + NA_Nonsyn + M_Syn+  (1 | vaccine_code), 
                         data = train_regdata_4, 
                         family = binomial)

Final_model_4_3 <- glmer(Dom ~ HA_Nonsyn * M_Syn + PA_Nonsyn + NA_Nonsyn +  (1 | vaccine_code), 
                         data = train_regdata_4, 
                         family = binomial) 

Final_model_4_4 <- glmer(Dom ~ NA_Nonsyn * PA_Nonsyn + HA_Nonsyn + M_Syn + (1 | vaccine_code), 
                         data = train_regdata_4, 
                         family = binomial) 

Final_model_4_5 <- glmer(Dom ~ NA_Nonsyn *  M_Syn + PA_Nonsyn + HA_Nonsyn + (1 | vaccine_code), 
                         data = train_regdata_4, 
                         family = binomial)

Final_model_4_6 <- glmer(Dom ~ HA_Nonsyn + NA_Nonsyn + M_Syn * PA_Nonsyn +  (1 | vaccine_code), 
                         data = train_regdata_4, 
                         family = binomial) 

# Compare models using AIC
AIC(Base_model_4, Final_model_4, Final_model_4_1, Final_model_4_2, Final_model_4_3, Final_model_4_4, Final_model_4_5, Final_model_4_6)

#### A. Predicted Score ####
Base_model_4_logit_P = predict(Base_model_4, newdata = train_regdata_4, type = 'response' )
train_regdata_4$Base_model_4_logit_P <- Base_model_4_logit_P

#### B. ROC ####
roc_score_Train_Basemodel_4 <- roc(train_regdata_4$Dom, Base_model_4_logit_P) #AUC score
## 0.7907 ##

#### A. Predicted Score -- CV ####
Base_model_4_logit_P_CV = predict(Base_model_4, newdata = test_regdata_4, type = 'response' )
test_regdata_4$Base_model_4_logit_P_CV <- Base_model_4_logit_P_CV

#### B. ROC ####
roc_score_Test_Basemodel_4_CV <- roc(test_regdata_4$Dom, Base_model_4_logit_P_CV) #AUC score
## 0.6587 ##

#### A. Predicted Score ####
Final_model_4_logit_P = predict(Final_model_4_1, newdata = train_regdata_4, type = 'response' )
train_regdata_4$Final_model_4_logit_P <- Final_model_4_logit_P

#### B. ROC ####
roc_score_Train_Finalmodel_4 <- roc(train_regdata_4$Dom, Final_model_4_logit_P) #AUC score
## 0.8023 ##

#### A. Predicted Score -- CV ####
Final_model_4_logit_P_CV = predict(Final_model_4_1, newdata = test_regdata_4, type = 'response' )
test_regdata_4$Final_model_4_logit_P_CV <- Final_model_4_logit_P_CV

#### B. ROC ####
roc_score_Test_Finalmodel_4_CV <- roc(test_regdata_4$Dom, Final_model_4_logit_P_CV) #AUC score
## 0.8651 ##
