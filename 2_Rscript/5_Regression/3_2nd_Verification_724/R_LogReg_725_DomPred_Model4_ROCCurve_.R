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


## Evaluate the 725 Data -- Data 1 ###
## 1. Model 1 ##
Base_model_1 <- glmer(Dom ~ HA_Nonsyn + (1 | vaccine_code), 
                       data = train_regdata_1, family = binomial)

Final_model_1 <- glmer(Dom ~ HA_Nonsyn + NA_Nonsyn * NS_Nonsyn + (1 | vaccine_code), 
                    data = train_regdata_1, family = binomial)


#### A. Predicted Score ####
Base_model_1_logit_P = predict(Base_model_1, newdata = train_regdata_1, type = 'response' )
train_regdata_1$Base_model_1_logit_P <- Base_model_1_logit_P

#### B. ROC ####
roc_score_Train_BaseModel_1 <- roc(train_regdata_1$Dom, Base_model_1_logit_P) #AUC score
## 0.7649 ##

jpeg(filename = "01_GenDisFlu/Fig/PredModel/1_HK15Model/1_Base_Model_Train.jpeg", width = 170, height = 150, units = "mm", res = 400)

plot(roc_score_Train_BaseModel_1, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="The Base model 1 with train data")

dev.off()


#### A. Predicted Score -- CV ####
Base_model_1_logit_P_CV = predict(Base_model_1, newdata = test_regdata_1, type = 'response' )
test_regdata_1$Base_model_1_logit_P_CV <- Base_model_1_logit_P_CV

#### B. ROC ####
roc_score_Test_Basemodel_1_CV <- roc(test_regdata_1$Dom, Base_model_1_logit_P_CV) #AUC score
## 0.8571 ##

jpeg(filename = "01_GenDisFlu/Fig/PredModel/1_HK15Model/1_Base_Model_Test.jpeg", width = 170, height = 150, units = "mm", res = 400)

plot(roc_score_Test_Basemodel_1_CV, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="The Base model 1 with test data")

dev.off()

#### A. Predicted Score ####
Final_model_1_logit_P = predict(Final_model_1, newdata = train_regdata_1, type = 'response' )
train_regdata_1$Final_model_1_logit_P <- Final_model_1_logit_P

#### B. ROC ####
roc_score_Train_Finalmodel_1 <- roc(train_regdata_1$Dom, Final_model_1_logit_P) #AUC score
## 0.8246 ##

jpeg(filename = "01_GenDisFlu/Fig/PredModel/1_HK15Model/1_Final_Model_train.jpeg", width = 170, height = 150, units = "mm", res = 400)

plot(roc_score_Train_Finalmodel_1, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="The Final model 1 with train data")

dev.off()

#### A. Predicted Score -- CV ####
Final_model_1_logit_P_CV = predict(Final_model_1, newdata = test_regdata_1, type = 'response' )
test_regdata_1$Final_model_1_logit_P_CV <- Final_model_1_logit_P_CV

#### B. ROC ####
roc_score_Test_Finalmodel_1_CV <- roc(test_regdata_1$Dom, Final_model_1_logit_P_CV) #AUC score
## 0.8946 ##

jpeg(filename = "01_GenDisFlu/Fig/PredModel/1_HK15Model/1_Final_Model_Test.jpeg", width = 170, height = 150, units = "mm", res = 400)

plot(roc_score_Test_Finalmodel_1_CV, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="The Final model 1 with test data")

dev.off()


## 2. Model 2 ##
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

jpeg(filename = "01_GenDisFlu/Fig/PredModel/2_Kan17Model/2_Base_Model_train.jpeg", width = 170, height = 150, units = "mm", res = 400)

plot(roc_score_Train_BaseModel_2, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="The Base model 2 with train data")

dev.off()


#### A. Predicted Score -- CV ####
Base_model_2_logit_P_CV = predict(Base_model_2, newdata = test_regdata_2, type = 'response' )
test_regdata_2$Base_model_2_logit_P_CV <- Base_model_2_logit_P_CV

#### B. ROC ####
roc_score_Test_Basemodel_2_CV <- roc(test_regdata_2$Dom, Base_model_2_logit_P_CV) #AUC score
## 0.8438 ##

jpeg(filename = "01_GenDisFlu/Fig/PredModel/2_Kan17Model/2_Base_Model_test.jpeg", width = 170, height = 150, units = "mm", res = 400)

plot(roc_score_Test_Basemodel_2_CV, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="The Base model 2 with test data")

dev.off()


#### A. Predicted Score ####
Final_model_2_logit_P = predict(Final_model_2, newdata = train_regdata_2, type = 'response' )
train_regdata_2$Final_model_2_logit_P <- Final_model_2_logit_P

#### B. ROC ####
roc_score_Train_Finalmodel_2 <- roc(train_regdata_2$Dom, Final_model_2_logit_P) #AUC score
## 0.819 ##

jpeg(filename = "01_GenDisFlu/Fig/PredModel/2_Kan17Model/2_Final_Model_train.jpeg", width = 170, height = 150, units = "mm", res = 400)

plot(roc_score_Train_Finalmodel_2, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="The Final model 2 with train data")

dev.off()


#### A. Predicted Score -- CV ####
Final_model_2_logit_P_CV = predict(Final_model_2, newdata = test_regdata_2, type = 'response' )
test_regdata_2$Final_model_2_logit_P_CV <- Final_model_2_logit_P_CV

#### B. ROC ####
roc_score_Test_Finalmodel_2_CV <- roc(test_regdata_2$Dom, Final_model_2_logit_P_CV) #AUC score
## 0.5375 ##

jpeg(filename = "01_GenDisFlu/Fig/PredModel/2_Kan17Model/2_Final_Model_test.jpeg", width = 170, height = 150, units = "mm", res = 400)

plot(roc_score_Test_Finalmodel_2_CV, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="The Final model 2 with test data")

dev.off()



## 3. Model 3 ##
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

jpeg(filename = "01_GenDisFlu/Fig/PredModel/3_HK19Model/3_Base_Model_train.jpeg", width = 170, height = 150, units = "mm", res = 400)

plot(roc_score_Train_BaseModel_3, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="The Base model 3 with train data")

dev.off()


#### A. Predicted Score -- CV ####
Base_model_3_logit_P_CV = predict(Base_model_3, newdata = test_regdata_3, type = 'response' )
test_regdata_3$Base_model_3_logit_P_CV <- Base_model_3_logit_P_CV

#### B. ROC ####
roc_score_Test_Basemodel_3_CV <- roc(test_regdata_3$Dom, Base_model_3_logit_P_CV) #AUC score
## 0.8519 ##

jpeg(filename = "01_GenDisFlu/Fig/PredModel/3_HK19Model/3_Base_Model_test.jpeg", width = 170, height = 150, units = "mm", res = 400)

plot(roc_score_Test_Basemodel_3_CV, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="The Base model 3 with test data")

dev.off()

#### A. Predicted Score ####
Final_model_3_logit_P = predict(Final_model_3, newdata = train_regdata_3, type = 'response' )
train_regdata_3$Final_model_3_logit_P <- Final_model_3_logit_P

#### B. ROC ####
roc_score_Train_Finalmodel_3 <- roc(train_regdata_3$Dom, Final_model_3_logit_P) #AUC score
## 0.8144 ##

jpeg(filename = "01_GenDisFlu/Fig/PredModel/3_HK19Model/3_Final_Model_train.jpeg", width = 170, height = 150, units = "mm", res = 400)

plot(roc_score_Train_Finalmodel_3, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="The Final model 3 with train data")

dev.off()


#### A. Predicted Score -- CV ####
Final_model_3_logit_P_CV = predict(Final_model_3, newdata = test_regdata_3, type = 'response' )
test_regdata_3$Final_model_3_logit_P_CV <- Final_model_3_logit_P_CV

#### B. ROC ####
roc_score_Test_Finalmodel_3_CV <- roc(test_regdata_3$Dom, Final_model_3_logit_P_CV) #AUC score
## 0.6296 ##

jpeg(filename = "01_GenDisFlu/Fig/PredModel/3_HK19Model/3_Final_Model_test.jpeg", width = 170, height = 150, units = "mm", res = 400)

plot(roc_score_Test_Finalmodel_3_CV, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="The Final model 3 with test data")

dev.off()


## 4. Model 4 ##
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

jpeg(filename = "01_GenDisFlu/Fig/PredModel/4_Cam20Model/4_Base_Model_train.jpeg", width = 170, height = 150, units = "mm", res = 400)

plot(roc_score_Train_BaseModel_4, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="The Base model 4 with train data")

dev.off()



#### A. Predicted Score -- CV ####
Base_model_4_logit_P_CV = predict(Base_model_4, newdata = test_regdata_4, type = 'response' )
test_regdata_4$Base_model_4_logit_P_CV <- Base_model_4_logit_P_CV

#### B. ROC ####
roc_score_Test_Basemodel_4_CV <- roc(test_regdata_4$Dom, Base_model_4_logit_P_CV) #AUC score
## 0.6587 ##

jpeg(filename = "01_GenDisFlu/Fig/PredModel/4_Cam20Model/4_Base_Model_test.jpeg", width = 170, height = 150, units = "mm", res = 400)

plot(roc_score_Test_Basemodel_4_CV, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="The Base model 4 with test data")

dev.off()


#### A. Predicted Score ####
Final_model_4_logit_P = predict(Final_model_4, newdata = train_regdata_4, type = 'response' )
train_regdata_4$Final_model_4_logit_P <- Final_model_4_logit_P

#### B. ROC ####
roc_score_Train_Finalmodel_4 <- roc(train_regdata_4$Dom, Final_model_4_logit_P) #AUC score
## 0.8187 ##

jpeg(filename = "01_GenDisFlu/Fig/PredModel/4_Cam20Model/4_Final_Model_train.jpeg", width = 170, height = 150, units = "mm", res = 400)

plot(roc_score_Train_Finalmodel_4, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="The Final model 4 with train data")

dev.off()

#### A. Predicted Score -- CV ####
Final_model_4_logit_P_CV = predict(Final_model_4, newdata = test_regdata_4, type = 'response' )
test_regdata_4$Final_model_4_logit_P_CV <- Final_model_4_logit_P_CV

#### B. ROC ####
roc_score_Test_Finalmodel_4_CV <- roc(test_regdata_4$Dom, Final_model_4_logit_P_CV) #AUC score
## 0.8095 ##

jpeg(filename = "01_GenDisFlu/Fig/PredModel/4_Cam20Model/4_Final_Model_test.jpeg", width = 170, height = 150, units = "mm", res = 400)

plot(roc_score_Test_Finalmodel_4_CV, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE, print.auc.cex= 1.5,
     print.thres=TRUE, print.thres.pch=21, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="The Final model 4 with test data")

dev.off()


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

jpeg(filename = "01_GenDisFlu/Fig/PredModel/ModelPred.jpeg", width = 220, height = 110, units = "mm", res = 720)

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
