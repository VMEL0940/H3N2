### Seasonal H3N2 Prediction - N and S Genetic distance descriptive ###
## Package ##
library(tidyverse)

## Set working directory ##
setwd("~/H3N2/")

## Data import ##
AllG_NS <- read.csv('1_Data/2_RNSC/3_DistanceLog/AllG_NS.csv', header = T, na.strings = "")

### All Group ###

AllG_NS$Gene <- factor(AllG_NS$Gene, levels = c("PB2", "PB1", "PA", "HA", 
                                                "NP", "NA", "M", "NS"))

GeneCol <- c("#C780E8", "#9D94FF","#59ADF6","#08CAD1","#42D6A4","#F8F38D","#FFB480","#FF6961")

## Descriptive Summary of N & S value by Gene and Vaccine period ##
Desc_all <- AllG_NS %>% 
  group_by(Gene, vaccine_code) %>% 
  summarise(Nmean = mean(Nmedian), Nse = 1.96*sd(Nmedian)/sqrt(n()),
            Smean = mean(Smedian), Sse = 1.96*sd(Smedian)/sqrt(n()),
            n = n())

Desc_all_N_All <- AllG_NS %>% 
  group_by(Gene) %>% 
  summarise(DisMean = mean(Nmedian), DisMedian = median(Nmedian),
            lowerCrIn = quantile(Nmedian, na.rm = T, probs = 0.05),
            UpperCrIn = quantile(Nmedian, na.rm = T, probs = 0.95),
            lower = mean(Nmedian) - 1.96*sd(Nmedian)/sqrt(n()),
            upper = mean(Nmedian) + 1.96*sd(Nmedian)/sqrt(n()),
            n = n())


Desc_all_S_All <- AllG_NS %>% 
  group_by(Gene) %>% 
  summarise(DisMean = mean(Smedian), DisMedian = median(Smedian),
            lowerCrIn = quantile(Smedian, na.rm = T, probs = 0.05),
            UpperCrIn = quantile(Smedian, na.rm = T, probs = 0.95),
            lower = mean(Smedian) - 1.96*sd(Smedian)/sqrt(n()),
            upper = mean(Smedian) + 1.96*sd(Smedian)/sqrt(n()),
            n = n())



Desc_all_N <- AllG_NS %>% 
  group_by(Gene, vaccine_code) %>% 
  summarise(DisMean = mean(Nmedian), DisMedian = median(Nmedian),
            lowerCrIn = quantile(Nmedian, na.rm = T, probs = 0.05),
            UpperCrIn = quantile(Nmedian, na.rm = T, probs = 0.95),
            lower = mean(Nmedian) - 1.96*sd(Nmedian)/sqrt(n()),
            upper = mean(Nmedian) + 1.96*sd(Nmedian)/sqrt(n()),
            n = n())


Desc_all_S <- AllG_NS %>% 
  group_by(Gene, vaccine_code) %>% 
  summarise(DisMean = mean(Smedian), DisMedian = median(Smedian),
            lowerCrIn = quantile(Smedian, na.rm = T, probs = 0.05),
            UpperCrIn = quantile(Smedian, na.rm = T, probs = 0.95),
            lower = mean(Smedian) - 1.96*sd(Smedian)/sqrt(n()),
            upper = mean(Smedian) + 1.96*sd(Smedian)/sqrt(n()),
            n = n())


## Box plot for S value by Gene ##
p <- AllG_NS %>% 
    ggplot() +
    geom_jitter(aes(x= Gene, y= Smedian), width = 0.25, size = 0.6, color = "gray40") +
    geom_boxplot(aes(x= Gene, y = Smedian), color = GeneCol, size = 1.3, outlier.shape = NA, alpha = 0) +
    scale_color_manual(values = GeneCol, guide="none") +
    scale_fill_manual(values = GeneCol, guide="none") +
    ylab("Synonymous genetic distance")+
    xlab("Gene segment")+
    ylim(-1, 65) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12)) 


ggsave("Synonymous_Distance.svg", plot = p, device = "svg", 
       path = "3_Figures/Fig_2/Fig2A_SynDist",
       width = 700, height = 300, units = "mm", dpi = 720)

## Box plot for N value by Gene ##
p <- AllG_NS %>% 
    ggplot() +
    geom_jitter(aes(x= Gene, y= Nmedian), width = 0.25, size = 0.50, color = "gray40") +
    geom_boxplot(aes(x= Gene, y = Nmedian), color = GeneCol, size = 1.3, outlier.shape = NA, alpha = 0) +
    scale_color_manual(values = GeneCol, guide="none") +
    scale_fill_manual(values = GeneCol, guide="none") +
    ylab("Nonsynonymous genetic distance")+
    xlab("Gene segment")+
    ylim(-1, 60) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12)) 


ggsave("NonSynonymous_Distance.svg", plot = p, device = "svg", 
       path = "3_Figures/Fig_2/Fig2B_NonSynDist",
       width = 700, height = 300, units = "mm", dpi = 720)




################ For reference ################

## Box + Jitter plot for N value by Gene and Vaccine period ##
AllG_NS %>% 
  ggplot() +
  geom_jitter(aes(x= Gene, y= Nmedian), width = 0.25, size = 0.75, color = "gray40") +
  geom_boxplot(aes(x= Gene, y = Nmedian, color = Gene), size = 0.8, outlier.shape = NA, alpha = 0) +
  scale_color_manual(values = GeneCol, guide="none") +
  scale_fill_manual(values = GeneCol, guide="none") +
  facet_wrap( ~ vaccine_code, nrow = 3, ncol = 5) +
  ylab("Nonsynonymous genetic distance")+
  xlab("Gene segment")+
  theme_bw() +
  ylim(-1, 65) +
  theme(strip.background = element_rect(fill="white"), 
        strip.text.x = element_text(size = 10, face = "bold", family = "sans"),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 7))


## Box + Jitter plot for S value by Gene and Vaccine period ##
AllG_NS %>% 
  ggplot() +
  geom_jitter(aes(x= Gene, y= Smedian), width = 0.25, size = 0.50, color = "gray40") +
  geom_boxplot(aes(x= Gene, y = Smedian, color = Gene), size = 0.8, outlier.shape = NA, alpha = 0) +
  scale_color_manual(values = GeneCol, guide="none") +
  scale_fill_manual(values = GeneCol, guide="none") +
  facet_wrap( ~ vaccine_code, nrow = 3, ncol = 5) +
  ylab("Synonymous genetic distance")+
  xlab("Gene segment")+
  theme_bw() +
  ylim(-1, 65) +
  theme(strip.background = element_rect(fill="white"), 
        strip.text.x = element_text(size = 10, face = "bold", family = "sans"),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 7))



## Error Bar plot for N value by Gene and Vaccine period ##
Desc_all %>% 
  ggplot() +
  geom_errorbar(aes(x= Gene, ymin=Nmean-Nse,ymax=Nmean+Nse), color = "gray40",
                width=0.3, size = 0.7) +
  geom_point(aes(x= Gene, y=Nmean, fill = Gene), size = 2, shape = 21, color = "gray40") +
  scale_fill_manual(values = GeneCol, guide="none") +
  facet_wrap(~vaccine_code, nrow = 3, ncol = 5) +
  ylab("Renaissance count of Non-synonymous mutation")+
  ylim(0, 40) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 7),
        strip.text = element_text(size=11)) 


## Error Bar plot for S value by Gene and Vaccine period ##
Desc_all %>% 
  ggplot() +
  geom_errorbar(aes(x= Gene, ymin=Smean-Sse, ymax=Smean+Sse), color = "gray40",
                width=0.3, size = 0.7) +
  geom_point(aes(x= Gene, y=Smean, fill = Gene), size = 2, shape = 21, color = "gray40") +
  scale_fill_manual(values = GeneCol, guide="none") +
  facet_wrap(~vaccine_code, nrow = 2, ncol = 5) +
  ylab("Renaissance count of Synonymous mutation")+
  ylim(0, 40) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 7),
        strip.text = element_text(size=11)) 


## Error Bar plot for N value by Gene ##
AllG_NS %>% 
  group_by(Gene) %>% 
  summarise(Nmean = mean(Nmedian), Nse = 1.96*sd(Nmedian)/sqrt(n()),
            Smean = mean(Smedian), Sse = 1.96*sd(Smedian)/sqrt(n()),
            n = n()) %>% 
  ggplot() +
  geom_errorbar(aes(x= Gene, ymin=Nmean-Nse,ymax=Nmean+Nse), color = "gray40",
                width=0.2, size = 0.7) +
  geom_point(aes(x= Gene, y=Nmean, fill = Gene), size = 3, shape = 21, color = "gray40") +
  scale_fill_manual(values = GeneCol, guide="none") +
  ylab("Renaissance count of Non-synonymous mutation")+
  ylim(0, 15) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12)) 


## Error Bar plot for S value by Gene ##
AllG_NS %>% 
  group_by(Gene) %>% 
  summarise(Nmean = mean(Nmedian), Nse = 1.96*sd(Nmedian)/sqrt(n()),
            Smean = mean(Smedian), Sse = 1.96*sd(Smedian)/sqrt(n()),
            n = n()) %>% 
  ggplot() +
  geom_errorbar(aes(x= Gene, ymin=Smean-Sse, ymax=Smean+Sse), color = "gray40",
                width=0.2, size = 0.7) +
  geom_point(aes(x= Gene, y=Smean, fill = Gene), size = 3, shape = 21, color = "gray40") +
  scale_fill_manual(values = GeneCol, guide="none") +
  ylab("Renaissance count of Synonymous mutation")+
  ylim(0, 15) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12)) 


