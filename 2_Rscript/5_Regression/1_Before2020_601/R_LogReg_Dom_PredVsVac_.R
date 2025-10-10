## Package ##
library(tidyverse)
library(lme4)
library(pROC)
library(lmtest)
library(corrplot)

## Set working directory ##
setwd("/01_GenDisFlu/")
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


##### 2 Base model - Koel #####
Koel_model <- glmer(subtMRCA_1stDesc ~ HA_Koel + (1 | Vaccine_code), 
                       data=train_regdata, family=binomial(link="logit"))

# Summarize the model
summary(Koel_model)

#### A. Predicted Score ####
Koel_logit_P = predict(Koel_model, newdata = train_regdata, type = 'response' )
train_regdata$Koel_logit_P <- Koel_logit_P


############# Model selection #############
##### 1.HA only model #####

GD_1_HA_model <- glmer(subtMRCA_1stDesc ~ HA + (1 | Vaccine_code), 
                       data=train_regdata, family=binomial(link="logit"))    

summary(GD_1_HA_model)

#### A. Predicted Score ####
HA_logit_P = predict(GD_1_HA_model, newdata = train_regdata, type = 'response' )
train_regdata$HA_logit_P <- HA_logit_P



##### 2.HA + 1 Model #####
#### 2.3 HA + NA Model #####
GD_2_3_HA_NA_model <- glmer(subtMRCA_1stDesc ~ HA+ N_A + (1 | Vaccine_code), 
                           data=train_regdata, family=binomial(link="logit"))    ### Sig

summary(GD_2_3_HA_NA_model)

#### A. Predicted Score ####
HA_NA_logit_P = predict(GD_2_3_HA_NA_model, newdata = train_regdata, type = 'response' )
train_regdata$HA_NA_logit_P <- HA_NA_logit_P




#### 2.3.1 HA*PA + NA Model #####
GD_2_3_1_HAPA_NA_model <- glmer(subtMRCA_1stDesc ~ HA*PA + N_A + (1 | Vaccine_code), 
                            data=train_regdata, family=binomial(link="logit"))    ### Sig

summary(GD_2_3_1_HAPA_NA_model)

#### A. Predicted Score ####
HAPA_NA_logit_P = predict(GD_2_3_1_HAPA_NA_model, newdata = train_regdata, type = 'response' )
train_regdata$HAPA_NA_logit_P <- HAPA_NA_logit_P


#### 2.3.5 HA*PA + NA*NS Model #####
GD_2_3_5_HAPA_NANS_model <- glmer(subtMRCA_1stDesc ~ HA*PA + N_A*NS + (1 | Vaccine_code), 
                                data=train_regdata, family=binomial(link="logit"))    ### Sig

summary(GD_2_3_5_HAPA_NANS_model)

#### A. Predicted Score ####
HAPA_NANS_logit_P = predict(GD_2_3_5_HAPA_NANS_model, newdata = train_regdata, type = 'response' )
train_regdata$HAPA_NANS_logit_P <- HAPA_NANS_logit_P

#write.csv(train_regdata, "Data/Pred_Result.csv")

### Vaccine strain ###
train_regdata_V2 <- train_regdata %>% 
  mutate(Vaccine = if_else(str_detect(ID, "EPI"), 1, 0))

train_regdata_V2 %>% 
  filter(Vaccine == 1) %>% 
  select(ID, RBD_15A_logit_P, Koel_logit_P, HA_logit_P, HA_NA_logit_P, HAPA_NA_logit_P, HAPA_NANS_logit_P, Vaccine_code)

### 1. Boxplot of Dom by model ###
train_regdata_V3 <- train_regdata_V2 %>% 
  gather("Model", "Score", 27:32)

train_regdata_V3$Pred_Model <- factor(train_regdata_V3$Model, levels = c("RBD_15A_logit_P", "Koel_logit_P", "HA_logit_P", 
                                                                         "HA_NA_logit_P", "HAPA_NA_logit_P", "HAPA_NANS_logit_P"),
                                      labels = c("RBD&15A", "Koel","HA", "HA&NA", "HA&NA&PA", "HA&NA&PA&NS"))
### Visualization ###
### 1. Boxplot of Dom by model ###

VP3Col <- c("#030E4F","#F49F1C")


jpeg(filename = "/01_GenDisFlu/Fig/PredModel/Dom_Pred_Boxplot.jpeg", width = 180, height = 150, units = "mm", res = 720)

train_regdata_V3 %>% 
  ggplot() +  theme_bw() +
  geom_point(position = position_jitterdodge(jitter.width = 0.2,
                                             dodge.width = 0.7), aes(x= Pred_Model, y = Score, color = as.factor(subtMRCA_1stDesc))) +
  geom_boxplot(aes(x= Pred_Model, y = Score, color = as.factor(subtMRCA_1stDesc)), size = 0.8, outlier.shape = NA, alpha = 0) +
  scale_color_manual(values = VP3Col, name="Dominance", labels=c('Extinct', 'Dominant')) + 
  xlab("Models") +
  ylab("Prediction score")

dev.off()


### Time!!! ###
Subsetlist <- read.csv("Data/H3N2_Index_20230331.csv")
names(train_regdata_V2)[1] <- "compareStrain"
train_regdata_V4 <- full_join(train_regdata_V2, Subsetlist[,c(1,4:6)], by = "compareStrain")
train_regdata_V4$ColDate <- paste0(train_regdata_V4$Year, "/", train_regdata_V4$Month, "/", train_regdata_V4$Day)
train_regdata_V4$ColDate <- as.Date(train_regdata_V4$ColDate,"%Y/%m/%d")


## Time and vaccine ##
train_regdata_V4 %>%
  ggplot() +
  geom_point(aes(x = ColDate, y = HAPA_NANS_logit_P, color = as.factor(Vaccine)), size = 1.7) +
  geom_vline(xintercept = as.Date(c("2004-10-1", "2005-10-1",
                                    "2006-10-1", "2008-10-1",
                                    "2010-10-1", "2012-10-1", 
                                    "2015-10-1", "2016-10-1")), 
             color = "red", lty = 2, 
             lwd = 0.4) +
  scale_x_date(name = "Year", minor_breaks = "1 year") +
  ylab("Score of HA&NA&PA&NS") +
  scale_color_manual(values = VP3Col, guide="none") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9))


## Time and Dominance ##
train_regdata_V4 %>%
  ggplot() +
  geom_point(aes(x = ColDate, y = HA_NA_logit_P, color = as.factor(subtMRCA_1stDesc)), size = 1.7) +
  geom_vline(xintercept = as.Date(c("2004-10-1", "2005-10-1",
                                    "2006-10-1", "2008-10-1",
                                    "2010-10-1", "2012-10-1", 
                                    "2015-10-1", "2016-10-1")), 
             color = "red", lty = 2, 
             lwd = 0.4) +
  scale_x_date(name = "Year", minor_breaks = "1 year") +
  ylab("Score of HA&NA&PA&NS") +
  scale_color_manual(values = VP3Col, guide="none") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9))




