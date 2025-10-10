## Package ##
library(tidyverse)
library(lme4)
library(pROC)
## Set working directory ##
setwd("/01_GenDisFlu/")

## Data import ##
AllG_NS <- read.csv("Data/Genetic Distances/GeneCompete/AllG_NS.csv", header = T, na.strings = "")

AllG_NS$compareStrain_V1 <- gsub("_NA_NA_NA_NA", "", AllG_NS$compareStrain)
AllG_NS$compareStrain_V1 <- gsub("_NA_NA_NA", "", AllG_NS$compareStrain_V1)
AllG_NS$compareStrain_V1 <- gsub("_NA_NA", "", AllG_NS$compareStrain_V1)
AllG_NS$compareStrain <- gsub("2009_NA", "2009", AllG_NS$compareStrain_V1)

Ind <- AllG_NS[,c(1,2,4,6)] %>% 
  spread(key = Gene, value = Nmedian)

colnames(Ind)[2] <- "ID"

## leftjoin * subset ##
Dom <- read.csv("Data/Genetic Distances/GeneCompete/Dom.csv", header = T, na.strings = "")

regdata <- right_join(Ind, Dom, by = "ID")

colnames(regdata)[5] <- "N_A"

Dom_Prop <- regdata %>% 
  group_by(Vaccine_code, Dominance) %>% 
  summarise(n = n())


### 1. FUll data set ##
train_regdata <- regdata %>% 
  filter(Vaccine_code != "HK15")

## multilevel multivaraible logistic regression ##
Full_ME_model <- glmer(Dominance ~ PB2+ PB1 + PA + HA + NP + N_A + M + NS + 
                         Reasrt + HA_RBD + HA_15A + 
                         CoEvo_PB2 + (1 | Vaccine_code), 
                       data=train_regdata, family=binomial(link="logit"))

# Summarize the model
summary(Full_ME_model)

Full_logit_P = predict(Full_ME_model, newdata = train_regdata, type = 'response' )
train_regdata$logit_P <- logit_P

jpeg(filename = "/01_GenDisFlu/Presentation/Fig/PredModel/ROC_Full.jpeg", width = 180, height = 160, units = "mm", res = 520)
#ROC-curve using pROC library
roc_score=roc(train_regdata$Dominance, Full_logit_P) #AUC score
plot(roc_score, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE,
     print.thres=TRUE, print.thres.pch=19, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="1. ROC curve of the Full Dominance Prediction Model")
dev.off()
### Final model ##
## multilevel multivaraible logistic regression ##
ME_model <- glmer(Dominance ~ PA + HA + N_A + M + NS +
                    HA_RBD + HA_15A + (1 | Vaccine_code), 
                       data=train_regdata, family=binomial(link="logit"))

# Summarize the model
summary(ME_model)

Rst_logit_P = predict(ME_model, newdata = train_regdata, type = 'response' )
train_regdata$Rst_logit_P <- Rst_logit_P

jpeg(filename = "/01_GenDisFlu/Presentation/Fig/PredModel/ROC_final.jpeg", width = 180, height = 160, units = "mm", res = 520)
#ROC-curve using pROC library
roc_score=roc(train_regdata$Dominance, Rst_logit_P) #AUC score
plot(roc_score, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE,
     print.thres=TRUE, print.thres.pch=19, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="2. ROC curve of the Final Dominance Prediction Model")
dev.off()


Coef_ME <- fixef(ME_model)
CI_ME <- confint(ME_model)

ranef(Full_ME_model)
ranef(ME_model)

OR <- round(exp(Coef_ME),2)
OR_CI <- round(exp(CI_ME),2)




### 1. RBD Only model  ##
## multilevel multivaraible logistic regression ##
RBD_ME_model <- glmer(Dominance ~ HA_RBD + (1 | Vaccine_code), 
                       data=train_regdata, family=binomial(link="logit"))

RBD15A_ME_model <- glmer(Dominance ~ HA_RBD + HA_15A + (1 | Vaccine_code), 
                       data=train_regdata, family=binomial(link="logit"))

HA_ME_model <- glmer(Dominance ~ HA + (1 | Vaccine_code), 
                           data=train_regdata, family=binomial(link="logit"))




# Summarize the model
summary(RBD_ME_model)
### Final model ##
summary(RBD15A_ME_model)
### Final model ##
summary(HA_ME_model)



RBD_logit_P = predict(RBD_ME_model, newdata = train_regdata, type = 'response' )
RBD15A_logit_P = predict(RBD15A_ME_model, newdata = train_regdata, type = 'response' )
HA_logit_P = predict(HA_ME_model, newdata = train_regdata, type = 'response' )

train_regdata$RBD_logit_P <- RBD_logit_P
train_regdata$RBD15A_logit_P <- RBD15A_logit_P
train_regdata$HA_logit_P <- HA_logit_P

jpeg(filename = "/01_GenDisFlu/Presentation/Fig/PredModel/ROC_RBD.jpeg", width = 220, height = 160, units = "mm", res = 520)

#ROC-curve using pROC library
roc_score=roc(train_regdata$Dominance, RBD_logit_P) #AUC score
roc_score
plot(roc_score, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE,
     print.thres=TRUE, print.thres.pch=19, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="3. The Dominance Prediction Model using RBD substitution")
dev.off()



jpeg(filename = "/01_GenDisFlu/Presentation/Fig/PredModel/ROC_RBD15A.jpeg", width = 220, height = 160, units = "mm", res = 520)
#ROC-curve using pROC library
roc_score=roc(train_regdata$Dominance, RBD15A_logit_P) #AUC score
roc_score
plot(roc_score, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE,
     print.thres=TRUE, print.thres.pch=19, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="4. The Dominance Prediction Model using RBD and 15A substitution")
dev.off()



jpeg(filename = "/01_GenDisFlu/Presentation/Fig/PredModel/ROC_HA.jpeg", width = 220, height = 160, units = "mm", res = 520)
#ROC-curve using pROC library
roc_score=roc(train_regdata$Dominance, HA_logit_P) #AUC score
roc_score
plot(roc_score, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE,
     print.thres=TRUE, print.thres.pch=19, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="5. The Dominance Prediction Model using HA distance")
dev.off()








### 1. FUll data set ##
test_regdata <- regdata %>% 
  filter(Vaccine_code == "HK15")

test_regdata$Vaccine_code <- "Vic11"

logit_P_tst = predict(Full_ME_model, newdata = test_regdata, type = 'response' )

hist(logit_P_tst)

test_regdata$logit_P <- logit_P_tst
test_regdata$Vaccine_code <- "HK15"


prd <- rbind(train_regdata, test_regdata)
hist(prd$logit_P)
write.csv(prd, "prd_Fullmodel.csv")

########## Cross validation #############




train.control <- trainControl(method = "repeatedcv", 
                              number = 5, repeats = 3)
# Train the model
model <- train(Fertility ~., data = swiss, method = "lm",
               trControl = train.control)
# Summarize the results
print(model)

### 2. Reassortment effect  ####

chisq.test(table(train_regdata2$Dominance, train_regdata2$Reasrt))

65/(168+65)
23/(183)

