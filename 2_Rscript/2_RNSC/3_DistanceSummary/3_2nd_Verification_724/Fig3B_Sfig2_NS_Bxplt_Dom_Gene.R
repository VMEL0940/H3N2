## Package ##
library(tidyverse)

## Set working directory ##
setwd("~/H3N2/")

## Evaluate the 601 Data -- Data 1 ###

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


### DATA Manupulation for the Nonsyn ##
All600_Train_NonSyn <- All600_Train[1:16]

Gat_Nonsyn <- All600_Train_NonSyn %>% 
  gather("Gene", "NonSyn", 9:16)

Gat_Nonsyn$vaccine_period <- factor(Gat_Nonsyn$vaccine_code, 
                                       levels = c("Mos99",
                                                  "Fuj02", 
                                                  "Cal04",
                                                  "Wis05",
                                                  "Bris07",
                                                  "Prth09", 
                                                  "Vic11",
                                                  "Swtz13"),
                                       labels = c("1.Mos99",
                                                  "2.Fuj02", 
                                                  "3.Cal04",
                                                  "4.Wis05",
                                                  "5.Bris07",
                                                  "6.Prth09", 
                                                  "7.Vic11",
                                                  "8.Swtz13"))

Gat_Nonsyn$Gene_F <- factor(Gat_Nonsyn$Gene, 
                                    levels = c("PB2_Nonsyn",
                                               "PB1_Nonsyn", 
                                               "PA_Nonsyn",
                                               "HA_Nonsyn",
                                               "NP_Nonsyn",
                                               "NA_Nonsyn", 
                                               "M_Nonsyn",
                                               "NS_Nonsyn"),
                                    labels = c("PB2",
                                               "PB1", 
                                               "PA",
                                               "HA",
                                               "NP",
                                               "NA", 
                                               "M",
                                               "NS"))

Gat_Nonsyn$Dominance <- factor(Gat_Nonsyn$Dom, 
                            levels = c(0,1),
                            labels = c("Extinct", "Dominant"))

VP2Col <- c("#000000","#96fd29")

### Box plot ###
p <- Gat_Nonsyn %>% 
    ggplot() +
    geom_jitter(aes(x= Dominance, y= NonSyn), width = 0.25, size = 0.50, color = "gray40") +
    geom_boxplot(aes(x= Dominance, y = NonSyn, color = Dominance), size = 1.1, outlier.shape = NA, alpha = 0) +
    xlab("Dominance") +
    ylab("Nonsynonymous genetic distance")+
    facet_wrap( ~ Gene_F, nrow = 2, ncol = 4) +
    scale_color_manual(values = VP2Col, guide="none") +
    theme_bw() +
    theme(strip.background = element_rect(fill="white"), 
          strip.text.x = element_text(size = 10, face = "bold", family = "sans"))


ggsave("NonSynonymous_Dom.svg", plot = p, device = "svg",
       path = "3_Figures/Fig_3/Fig3B_NonSynDom",
       width = 350, height = 450, units = "mm", dpi = 720)



### DATA Manupulation for the Syn ##
All600_Train_Syn <- All600_Train[c(1:8, 17:24)]

Gat_Syn <- All600_Train_Syn %>% 
  gather("Gene", "Syn", 9:16)

Gat_Syn$vaccine_period <- factor(Gat_Syn$vaccine_code, 
                                 levels = c("Mos99",
                                            "Fuj02", 
                                            "Cal04",
                                            "Wis05",
                                            "Bris07",
                                            "Prth09", 
                                            "Vic11",
                                            "Swtz13"),
                                 labels = c("1.Mos99",
                                            "2.Fuj02", 
                                            "3.Cal04",
                                            "4.Wis05",
                                            "5.Bris07",
                                            "6.Prth09", 
                                            "7.Vic11",
                                            "8.Swtz13"))

Gat_Syn$Gene_F <- factor(Gat_Syn$Gene, 
                         levels = c("PB2_Syn",
                                    "PB1_Syn", 
                                    "PA_Syn",
                                    "HA_Syn",
                                    "NP_Syn",
                                    "NA_Syn", 
                                    "M_Syn",
                                    "NS_Syn"),
                         labels = c("PB2",
                                    "PB1", 
                                    "PA",
                                    "HA",
                                    "NP",
                                    "NA", 
                                    "M",
                                    "NS"))

Gat_Syn$Dominance <- factor(Gat_Syn$Dom, 
                            levels = c(0,1),
                            labels = c("Extinct", "Dominant"))

### Box plot ###

p <- Gat_Syn %>% 
    ggplot() +
    geom_jitter(aes(x= Dominance, y= Syn), width = 0.25, size = 0.50, color = "gray40") +
    geom_boxplot(aes(x= Dominance, y = Syn, color = Dominance), size = 1.1, outlier.shape = NA, alpha = 0) +
    xlab("Dominance") +
    ylab("Synonymous genetic distance")+
    facet_wrap( ~ Gene_F, nrow = 2, ncol = 4) +
    scale_color_manual(values = VP2Col, guide="none") +
    theme_bw() +
    theme(strip.background = element_rect(fill="white"), 
          strip.text.x = element_text(size = 10, face = "bold", family = "sans"))


ggsave("Synonymous_Dom.svg", plot = p, device = "svg",
       path = "5_SFigure/Sfig2_SynDis",
       width = 450, height = 300, units = "mm", dpi = 720)

################ For reference ################

### DATA Manupulation for HA Mutations##
All600_Train_HA <- All600_Train[c(1:8, 25:27)]

Gat_HA <- All600_Train_HA %>% 
  gather("Position", "Mutation", 9:11)

Gat_HA$Position_F <- factor(Gat_HA$Position, 
                         levels = c("HA_RBD",
                                    "HA_15A", 
                                    "HA_Koel"),
                         labels = c("HA-RBD",
                                    "HA-RBD/15A", 
                                    "HA-Koel"))

Gat_HA$Dominance <- factor(Gat_HA$Dom, 
                            levels = c(0,1),
                            labels = c("Extinct", "Dominant"))

### Box plot ###

Gat_HA %>% 
  ggplot() +
  geom_jitter(aes(x= Dominance, y= Mutation), width = 0.25, size = 0.50, color = "gray40") +
  geom_boxplot(aes(x= Dominance, y = Mutation, color = Dominance), size = 1.1, outlier.shape = NA, alpha = 0) +
  xlab("Dominance") +
  ylab("The number of amino acid mutations")+
  facet_wrap( ~ Position_F, nrow = 1, ncol = 3) +
  scale_color_manual(values = VP2Col, guide="none") +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"), 
        strip.text.x = element_text(size = 10, face = "bold", family = "sans"))




