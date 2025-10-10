## Package ##
library(tidyverse)
library(lme4)
library(pROC)
## Set working directory ##
setwd("01_GenDisFlu/")
#setwd("01_GenDisFlu/")

## Data import ##
NewHA_NS <- read.csv("Evaluation/1st_Derek/HA242_NSdist.csv", header = T, na.strings = "")

colnames(NewHA_NS)[1] <- "ID"

NewHA_NS_reg <- NewHA_NS %>% 
  filter(Group != "G17")

NewHA_NS_reg$Vaccine_code <- factor(NewHA_NS_reg$Group,levels = c("G01","G02","G03","G04","G05","G06","G07","G08","G09",
                                                        "G10","G11","G12","G13","G14","G15","G16"),
                               labels = c("1.HK68","2.EN72", "3.PC73", "4.VI75","5.TE77","6.PH82","7.LE86","8.SI87","9.SH87",
                                          "10.GU89","11.BE89","12.BE92","13.SD93","14.JO94","15.WU95","16.SY97"))


### 1. FUll data set ##


## multilevel multivaraible logistic regression ##
Full_ME_model <- glmer(Dom ~ Nmedian + (1 | Group), 
                       data=NewHA_NS_reg, family=binomial(link="logit"))

# Summarize the model
summary(Full_ME_model)

Full_logit_P = predict(Full_ME_model, newdata = NewHA_NS_reg, type = 'response' )
NewHA_NS_reg$logit_P <- Full_logit_P

jpeg(filename = "01_GenDisFlu/Fig/PredModel/ROC_Evaldata.jpeg", width = 120, height = 110, units = "mm", res = 520)
#ROC-curve using pROC library
roc_score=roc(NewHA_NS_reg$Dom, Full_logit_P) #AUC score
plot(roc_score, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE,
     print.thres=TRUE, print.thres.pch=19, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="1968-2003 Dominance Prediction - 16 group")

dev.off()

exp(0.17003)

### 2. Balance data set ##
NewHA_NS_reg_bal <- NewHA_NS %>% 
  filter(Group == "G01" | Group == "G02" | Group == "G03" | Group == "G06" | Group == "G07" |
         Group == "G08" | Group == "G09" | Group == "G11" | Group == "G12" | Group == "G14" | Group == "G16" )

## multilevel multivaraible logistic regression ##
BAL_ME_model <- glmer(Dominance ~ Nmedian + (1 | Group), 
                      data=NewHA_NS_reg_bal, family=binomial(link="logit"))
                    

# Summarize the model
summary(BAL_ME_model)

BAL_logit_P = predict(BAL_ME_model, newdata = NewHA_NS_reg_bal, type = 'response' )
NewHA_NS_reg_bal$logit_P <- BAL_logit_P

jpeg(filename = "01_GenDisFlu/Fig/PredModel/ROC_UPDATE_data.jpeg", width = 180, height = 160, units = "mm", res = 520)
#ROC-curve using pROC library
roc_score=roc(NewHA_NS_reg_bal$Dominance, BAL_logit_P) #AUC score
plot(roc_score, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE,
     print.thres=TRUE, print.thres.pch=19, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="2. ROC curve of the Dominance Prediction Model using 11 groups")
dev.off()

exp(0.45613)



### 3. G11 data set ##
NewHA_NS_reg_G11 <- NewHA_NS %>% 
  filter(Group == "G11")

## multilevel multivaraible logistic regression ##
G11_ME_model <- glm(Dominance ~ Nmedian, 
                       data=NewHA_NS_reg_G11, family=binomial(link="logit"))

# Summarize the model
summary(G11_ME_model)

G11_logit_P = predict(G11_ME_model, newdata = NewHA_NS_reg_G11, type = 'response' )
NewHA_NS_reg_G11$logit_P <- G11_logit_P

jpeg(filename = "/01_GenDisFlu/Fig/PredModel/ROC_11onlydata.jpeg", width = 180, height = 160, units = "mm", res = 520)
#ROC-curve using pROC library
roc_score=roc(NewHA_NS_reg_G11$Dominance, G11_logit_P) #AUC score
plot(roc_score, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE,
     print.thres=TRUE, print.thres.pch=19, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="3. ROC curve of the Dominance Prediction Model only 11.BE89 group")
dev.off()

exp(0.6059)



### 1. FUll data set ##
NewHA_NS_reg <- NewHA_NS %>% 
  filter(Group != "G17")

## multilevel multivaraible logistic regression ##
Full_ME_model <- glmer(Dominance ~ Nmedian + (1 | Group), 
                       data=NewHA_NS_reg, family=binomial(link="logit"))

# Summarize the model
summary(Full_ME_model)

Full_logit_P = predict(Full_ME_model, newdata = NewHA_NS_reg, type = 'response' )
NewHA_NS_reg$logit_P <- Full_logit_P

jpeg(filename = "/01_GenDisFlu/Presentation/Fig/PredModel/ROC_Full.jpeg", width = 180, height = 160, units = "mm", res = 520)
#ROC-curve using pROC library
roc_score=roc(NewHA_NS_reg$Dominance, Full_logit_P) #AUC score
plot(roc_score, col="black",  
     print.auc=TRUE,  
     max.auc.polygon=TRUE,
     print.thres=TRUE, print.thres.pch=19, print.thres.col = "red",
     auc.polygon=TRUE, auc.polygon.col="#f8ce1c", 
     main ="1. ROC curve of the Full Dominance Prediction Model")
dev.off()


exp(0.6059)



NewHA_NS_reg %>%
  ggplot() +
  geom_boxplot(aes(x = as.factor(Dominance), y = Nmedian, fill = as.factor(Dominance))) +
  facet_wrap(~ Group, nrow = 3, ncol = 6) +
  theme_bw() 
               
NewHA_NS_reg %>%
  ggplot() +
  geom_jitter(aes(x = as.factor(Dominance), y = logit_P, color = as.factor(Dominance))) +
  facet_wrap(~ Group, nrow = 3, ncol = 6) +
  theme_bw() 

jpeg(filename = "/01_GenDisFlu/Fig/NS_Time/HA243_Eval_N_Dom_Boxplot.jpeg", width = 300, height = 250, units = "mm", res = 520)

NewHA_NS_reg %>%
  ggplot() +
  geom_jitter(aes(x = Dominance_F, y = Nmedian, color = Dominance_F), height = 0, width = 0.3) +
  geom_boxplot(aes(x = Dominance_F, y = Nmedian, color = Dominance_F), alpha = 0) +
  scale_color_discrete(name = "Dominance") +
  xlab("Dominance of the Influenza H3N2 strain") +
  ylab("Non-synonymous RNSC in the HA gene from Vaccine strain") +
  facet_wrap(~ Vaccine_code, nrow = 3, ncol = 6) +
  theme_bw() 

dev.off()
jpeg(filename = "/01_GenDisFlu/Fig/NS_Time/HA243_Eval_S_Dom_Boxplot.jpeg", width = 300, height = 250, units = "mm", res = 520)


NewHA_NS_reg %>%
  ggplot() +
  geom_jitter(aes(x = Dominance_F, y = Smedian, color = Dominance_F), height = 0, width = 0.3) +
  geom_boxplot(aes(x = Dominance_F, y = Smedian, color = Dominance_F), alpha = 0) +
  scale_color_discrete(name = "Dominance") +
  xlab("Dominance of the Influenza H3N2 strain") +
  ylab("Synonymous RNSC in the HA gene from Vaccine strain") +
  facet_wrap(~ Vaccine_code, nrow = 3, ncol = 6) +
  theme_bw() 

dev.off()
