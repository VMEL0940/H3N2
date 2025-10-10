### Seasonal H3N2 Prediction - N and S Genetic distance descriptive ###
## Package ##
library(tidyverse)

## Set working directory ##
setwd("~/H3N2/")

AllG_NS <- read.csv("1_Data/2_RNSC/1_MetaData/3_2nd_Verification_724/FullReg_724.csv", header = T, na.strings = "")

### 1. Filter -- Only after 2020 ## 
AllG_NS_Af2020 <- AllG_NS %>% 
  filter(vaccine_code == "Dar21" | vaccine_code == "Cam20" | vaccine_code == "HK19" | vaccine_code == "Kan17")

Smedian <- AllG_NS_Af2020[,c(24, 2:7, 9,11,13,15,17,19,21,23)] 

colnames(Smedian)[1]  <- c("ID")
colnames(Smedian)[8]  <- c("PB2")
colnames(Smedian)[9]  <- c("PB1")
colnames(Smedian)[10]  <- c("PA")
colnames(Smedian)[11]  <- c("HA")
colnames(Smedian)[12]  <- c("NP")
colnames(Smedian)[13]  <- c("NA")
colnames(Smedian)[14]  <- c("M")
colnames(Smedian)[15]  <- c("NS")

Des_Data_Af2020 <- Smedian %>% 
  gather("Gene", "SMed", 8:15)

NMedian <- AllG_NS_Af2020[,c(24, 2:7, 8,10,12,14,16,18,20,22)] 

colnames(NMedian)[1]  <- c("ID")
colnames(NMedian)[8]  <- c("PB2")
colnames(NMedian)[9]  <- c("PB1")
colnames(NMedian)[10]  <- c("PA")
colnames(NMedian)[11]  <- c("HA")
colnames(NMedian)[12]  <- c("NP")
colnames(NMedian)[13]  <- c("NA")
colnames(NMedian)[14]  <- c("M")
colnames(NMedian)[15]  <- c("NS")

Des_Data_Af2020_N <- NMedian %>% 
  gather("Gene", "NMed", 8:15)

Des_Data_Af2020_All <- cbind(Des_Data_Af2020, Des_Data_Af2020_N[c(9)])

Des_Data_Af2020_All$Gene <- factor(Des_Data_Af2020_All$Gene, levels = c("PB2", "PB1", "PA", "HA", 
                              "NP", "NA", "M", "NS"))

GeneCol <- c("#C780E8", "#9D94FF","#59ADF6","#08CAD1","#42D6A4","#F8F38D","#FFB480","#FF6961")


Des_Data_Af2020_All$vaccine_code_F <- factor(Des_Data_Af2020_All$vaccine_code, 
                                levels = c("Kan17",
                                           "HK19", 
                                           "Cam20",
                                           "Dar21"),
                                labels = c("10.Kan17",
                                           "11.HK19",
                                           "12.Cam20",
                                           "13.Dar21"))

### Only After 2020 ##


## Box plot for S value by Gene ##
p <- Des_Data_Af2020_All %>% 
    ggplot() +
    geom_jitter(aes(x= Gene, y= SMed), width = 0.25, size = 0.6, color = "gray40") +
    geom_boxplot(aes(x= Gene, y = SMed), color = GeneCol, size = 1.3, outlier.shape = NA, alpha = 0) +
    scale_color_manual(values = GeneCol, guide="none") +
    scale_fill_manual(values = GeneCol, guide="none") +
    ylab("Synonymous genetic distance")+
    xlab("Gene segment")+
    ylim(-1, 65) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12)) 


ggsave("Fig4B_SynDist_725.svg", plot = p, device = "svg",
       path = "3_Figures/Fig_4/Fig4B_SynDist_725",
       width = 400, height = 300, units = "mm", dpi = 720)

## Box plot for N value by Gene ##
p <- Des_Data_Af2020_All %>% 
    ggplot() +
    geom_jitter(aes(x= Gene, y= NMed), width = 0.25, size = 0.50, color = "gray40") +
    geom_boxplot(aes(x= Gene, y = NMed), color = GeneCol, size = 1.3, outlier.shape = NA, alpha = 0) +
    scale_color_manual(values = GeneCol, guide="none") +
    scale_fill_manual(values = GeneCol, guide="none") +
    ylab("Nonsynonymous genetic distance")+
    xlab("Gene segment")+
    ylim(-1, 65) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12)) 

ggsave("Fig4C_NonSynDist_725.svg", plot = p, device = "svg",
       path = "3_Figures/Fig_4/Fig4C_NonSynDist_725",
       width = 400, height = 300, units = "mm", dpi = 720)




################ For reference ################

## Descriptive Summary of N & S value by Gene and Vaccine period ##
Desc_all <- Des_Data_Af2020_All %>% 
  group_by(Gene, vaccine_code) %>% 
  summarise(Nmedian = mean(NMed), Nse = 1.96*sd(NMed)/sqrt(n()),
            Smedian = mean(SMed), Sse = 1.96*sd(SMed)/sqrt(n()),
            n = n())

Desc_all_N_All <- Des_Data_Af2020_All %>% 
  group_by(Gene) %>% 
  summarise(DisMean = mean(NMed), DiSMed = median(NMed),
            lowerCrIn = quantile(NMed, na.rm = T, probs = 0.05),
            UpperCrIn = quantile(NMed, na.rm = T, probs = 0.95),
            lower = mean(NMed) - 1.96*sd(NMed)/sqrt(n()),
            upper = mean(NMed) + 1.96*sd(NMed)/sqrt(n()),
            n = n())


Desc_all_S_All <- Des_Data_Af2020_All %>% 
  group_by(Gene) %>% 
  summarise(DisMean = mean(SMed), DiSMed = median(SMed),
            lowerCrIn = quantile(SMed, na.rm = T, probs = 0.05),
            UpperCrIn = quantile(SMed, na.rm = T, probs = 0.95),
            lower = mean(SMed) - 1.96*sd(SMed)/sqrt(n()),
            upper = mean(SMed) + 1.96*sd(SMed)/sqrt(n()),
            n = n())



Desc_all_N <- Des_Data_Af2020_All %>% 
  group_by(Gene, vaccine_code) %>% 
  summarise(DisMean = mean(NMed), DiSMed = median(NMed),
            lowerCrIn = quantile(NMed, na.rm = T, probs = 0.05),
            UpperCrIn = quantile(NMed, na.rm = T, probs = 0.95),
            lower = mean(NMed) - 1.96*sd(NMed)/sqrt(n()),
            upper = mean(NMed) + 1.96*sd(NMed)/sqrt(n()),
            n = n())


Desc_all_S <- Des_Data_Af2020_All %>% 
  group_by(Gene, vaccine_code) %>% 
  summarise(DisMean = mean(SMed), DiSMed = median(SMed),
            lowerCrIn = quantile(SMed, na.rm = T, probs = 0.05),
            UpperCrIn = quantile(SMed, na.rm = T, probs = 0.95),
            lower = mean(SMed) - 1.96*sd(SMed)/sqrt(n()),
            upper = mean(SMed) + 1.96*sd(SMed)/sqrt(n()),
            n = n())

write.csv(Desc_all_N, file="1_Data/2_RNSC/Desc_all_N.csv", row.names = F) 
write.csv(Desc_all_S, file="1_Data/2_RNSC/Desc_all_S.csv", row.names = F) 
write.csv(Desc_all_N_All, file="1_Data/2_RNSC/Desc_all_N_All.csv", row.names = F) 
write.csv(Desc_all_S_All, file="1_Data/2_RNSC/Desc_all_S_All.csv", row.names = F) 
write.csv(Desc_all, file="1_Data/2_RNSC/AllG_NS.csv", row.names = F)


