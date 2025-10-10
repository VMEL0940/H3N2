### Influ H3N2 GenDist Proj ###
## Package ##
library(tidyverse)
library(corrplot)
library(ggplotify)

## Set working directory ##
setwd("~/H3N2/")

## Data import ##
Nonsyn <- read.csv("1_Data/10_Final_MetaData_Submission/H3N2_Traindata_600strains_2_NonsynDist.csv", header = T, na.strings = "")
Syn <- read.csv("1_Data/10_Final_MetaData_Submission/H3N2_Traindata_600strains_3_SynDist.csv", header = T, na.strings = "")

colnames(Nonsyn)[c(4:11)] <- c("PB2", "PB1","PA","HA","NP", "NA","M","NS")
colnames(Syn)[c(4:11)] <- c("PB2", "PB1","PA","HA","NP", "NA","M","NS")

## Build Correlation matrix ##
M <- cor(Syn[,c(4:11)] , method = "pearson")
testRes = cor.mtest(Syn[,c(4:11)], conf.level = 0.50)

## Correlation plot
p <- as.ggplot(~corrplot(M, p.mat = testRes$p, method = "circle", insig = "blank",
                         tl.srt = 0, type = "upper",
                         sig.level = 0.05, addCoef.col = "black", number.cex = 0.8))

ggsave("Synonymous_Corr.svg", plot = p, device = "svg",
       path = "3_Figures/Fig_2/Fig2C_SynCorr",
       width = 350, height = 300, units = "mm", dpi = 720)


## Build Correlation matrix ##
M <- cor(Nonsyn[,c(4:11)] , method = "pearson")
testRes = cor.mtest(Nonsyn[,c(4:11)], conf.level = 0.50)

## Correlation plot
p <- as.ggplot(~corrplot(M, p.mat = testRes$p, method = "circle", insig = "blank",
                         tl.srt = 0, type = "upper",
                         sig.level = 0.05, addCoef.col = "black", number.cex = 0.8))

ggsave("NonSynonymous_Corr.svg", plot = p, device = "svg",
       path = "3_Figures/Fig_2/Fig2D_NonSynCorr",
       width = 350, height = 300, units = "mm", dpi = 720)

################ For reference ################

########## for G1 ###########

## Filter G1 of Nonsyn
G1_N <- Nonsyn %>% 
  filter(vaccine_code == "Mos99")

## Build Correlation matrix ##
M <- cor(G1_N[,c(4:11)] , method = "pearson")
testRes = cor.mtest(G1_N[,c(4:11)], conf.level = 0.50)


## Correlation plot
corrplot(M, p.mat = testRes$p, method = 'circle', insig='blank', tl.srt = 0, 
         sig.level = 0.05, addCoef.col ='black', number.cex = 0.8)

## Spread the table!! ##
G1_S <- Syn %>% 
  filter(vaccine_code == "Mos99")

## Build Correlation matrix ##
M <- cor(G1_S[,c(4:11)] , method = "pearson")
testRes = cor.mtest(G1_S[,c(4:11)], conf.level = 0.50)


## Correlation plot
corrplot(M, p.mat = testRes$p, method = 'circle', insig='blank', tl.srt = 0, 
         sig.level = 0.05, addCoef.col ='black', number.cex = 0.8)



########## for G2 ###########

## Filter G2 of Nonsyn
G2_N <- Nonsyn %>% 
  filter(vaccine_code == "Fuj02")

## Build Correlation matrix ##
M <- cor(G2_N[,c(4:11)] , method = "pearson")
testRes = cor.mtest(G2_N[,c(4:11)], conf.level = 0.95)


## Correlation plot
corrplot(M, p.mat = testRes$p, method = 'circle', insig='blank', tl.srt = 0, 
         sig.level = 0.05, addCoef.col ='black', number.cex = 0.8)


## Spread the table!! ##
G2_S <- Syn %>% 
  filter(vaccine_code == "Fuj02")

## Build Correlation matrix ##
M <- cor(G2_S[,c(4:11)] , method = "pearson")
testRes = cor.mtest(G2_S[,c(4:11)], conf.level = 0.50)


## Correlation plot
corrplot(M, p.mat = testRes$p, method = 'circle', insig='blank', tl.srt = 0, 
         sig.level = 0.05, addCoef.col ='black', number.cex = 0.8)


########## for G3 ###########

## Filter G3 of Nonsyn
G3_N <- Nonsyn %>% 
  filter(vaccine_code == "Cal04")

## Build Correlation matrix ##
M <- cor(G3_N[,c(4:11)] , method = "pearson")
testRes = cor.mtest(G3_N[,c(4:11)], conf.level = 0.95)


## Correlation plot
corrplot(M, p.mat = testRes$p, method = 'circle', insig='blank', tl.srt = 0, 
         sig.level = 0.05, addCoef.col ='black', number.cex = 0.8)


## Spread the table!! ##
G3_S <- Syn %>% 
  filter(vaccine_code == "Cal04")

## Build Correlation matrix ##
M <- cor(G3_S[,c(4:11)] , method = "pearson")
testRes = cor.mtest(G3_S[,c(4:11)], conf.level = 0.50)



## Correlation plot
corrplot(M, p.mat = testRes$p, method = 'circle', insig='blank', tl.srt = 0, 
         sig.level = 0.05, addCoef.col ='black', number.cex = 0.8)


########## for G4 ###########

## Filter G4 of Nonsyn
G4_N <- Nonsyn %>% 
  filter(vaccine_code == "Wis05")

## Build Correlation matrix ##
M <- cor(G4_N[,c(4:11)] , method = "pearson")
testRes = cor.mtest(G4_N[,c(4:11)], conf.level = 0.95)

## Correlation plot
corrplot(M, p.mat = testRes$p, method = 'circle', insig='blank', tl.srt = 0, 
         sig.level = 0.05, addCoef.col ='black', number.cex = 0.8)


## Spread the table!! ##
G4_S <- Syn %>% 
  filter(vaccine_code == "Wis05")

## Build Correlation matrix ##
M <- cor(G4_S[,c(4:11)] , method = "pearson")
testRes = cor.mtest(G4_S[,c(4:11)], conf.level = 0.50)



## Correlation plot
corrplot(M, p.mat = testRes$p, method = 'circle', insig='blank', tl.srt = 0, 
         sig.level = 0.05, addCoef.col ='black', number.cex = 0.8)





########## for G5 ###########

## Filter G5 of Nonsyn
G5_N <- Nonsyn %>% 
  filter(vaccine_code == "Bris07")

## Build Correlation matrix ##
M <- cor(G5_N[,c(4:11)] , method = "pearson")
testRes = cor.mtest(G5_N[,c(4:11)], conf.level = 0.95)

## Correlation plot
corrplot(M, p.mat = testRes$p, method = 'circle', insig='blank', tl.srt = 0, 
         sig.level = 0.05, addCoef.col ='black', number.cex = 0.8)


## Spread the table!! ##
G5_S <- Syn %>% 
  filter(vaccine_code == "Bris07")

## Build Correlation matrix ##
M <- cor(G5_S[,c(4:11)] , method = "pearson")
testRes = cor.mtest(G5_S[,c(4:11)], conf.level = 0.50)

## Correlation plot
corrplot(M, p.mat = testRes$p, method = 'circle', insig='blank', tl.srt = 0, 
         sig.level = 0.05, addCoef.col ='black', number.cex = 0.8)




########## for G6 ###########

## Filter G6 of Nonsyn
G6_N <- Nonsyn %>% 
  filter(vaccine_code == "Prth09")

## Build Correlation matrix ##
M <- cor(G6_N[,c(4:11)] , method = "pearson")
testRes = cor.mtest(G6_N[,c(4:11)], conf.level = 0.95)

## Correlation plot
corrplot(M, p.mat = testRes$p, method = 'circle', insig='blank', tl.srt = 0, 
         sig.level = 0.05, addCoef.col ='black', number.cex = 0.8)

## Spread the table!! ##
G6_S <- Syn %>% 
  filter(vaccine_code == "Prth09")

## Build Correlation matrix ##
M <- cor(G6_S[,c(4:11)] , method = "pearson")
testRes = cor.mtest(G6_S[,c(4:11)], conf.level = 0.50)

## Correlation plot
corrplot(M, p.mat = testRes$p, method = 'circle', insig='blank', tl.srt = 0, 
         sig.level = 0.05, addCoef.col ='black', number.cex = 0.8)



########## for G7 ###########

## Filter G7 of Nonsyn
G7_N <- Nonsyn %>% 
  filter(vaccine_code == "Vic11")

## Build Correlation matrix ##
M <- cor(G7_N[,c(4:11)] , method = "pearson")
testRes = cor.mtest(G7_N[,c(4:11)], conf.level = 0.95)

## Correlation plot
corrplot(M, p.mat = testRes$p, method = 'circle', insig='blank', tl.srt = 0, 
         sig.level = 0.05, addCoef.col ='black', number.cex = 0.8)

## Spread the table!! ##
G7_S <- Syn %>% 
  filter(vaccine_code == "Vic11")

## Build Correlation matrix ##
M <- cor(G7_S[,c(4:11)] , method = "pearson")
testRes = cor.mtest(G7_S[,c(4:11)], conf.level = 0.50)


## Correlation plot
corrplot(M, p.mat = testRes$p, method = 'circle', insig='blank', tl.srt = 0, 
         sig.level = 0.05, addCoef.col ='black', number.cex = 0.8)



########## for G8 ###########

## Filter G8 of Nonsyn
G8_N <- Nonsyn %>% 
  filter(vaccine_code == "Swtz13")

## Build Correlation matrix ##
M <- cor(G8_N[,c(4:11)] , method = "pearson")
testRes = cor.mtest(G8_N[,c(4:11)], conf.level = 0.95)

## Correlation plot
corrplot(M, p.mat = testRes$p, method = 'circle', insig='blank', tl.srt = 0, 
         sig.level = 0.05, addCoef.col ='black', number.cex = 0.8)


## Spread the table!! ##
G8_S <- Syn %>% 
  filter(vaccine_code == "Swtz13")

## Build Correlation matrix ##
M <- cor(G8_S[,c(4:11)] , method = "pearson")
testRes = cor.mtest(G8_S[,c(4:11)], conf.level = 0.50)



## Correlation plot
corrplot(M, p.mat = testRes$p, method = 'circle', insig='blank', tl.srt = 0, 
         sig.level = 0.05, addCoef.col ='black', number.cex = 0.8)

########## for G9 ###########

## Filter G9 of Nonsyn
G9_N <- Nonsyn %>% 
  filter(vaccine_code == "HK15")

## Build Correlation matrix ##
M <- cor(G9_N[,c(4:11)] , method = "pearson")
testRes = cor.mtest(G9_N[,c(4:11)], conf.level = 0.95)

## Correlation plot
corrplot(M, p.mat = testRes$p, method = 'circle', insig='blank', tl.srt = 0, 
         sig.level = 0.05, addCoef.col ='black', number.cex = 0.8)


## Spread the table!! ##
G9_S <- Syn %>% 
  filter(vaccine_code == "HK15")

## Build Correlation matrix ##
M <- cor(G9_S[,c(4:11)] , method = "pearson")
testRes = cor.mtest(G9_S[,c(4:11)], conf.level = 0.50)


## Correlation plot
corrplot(M, p.mat = testRes$p, method = 'circle', insig='blank', tl.srt = 0, 
         sig.level = 0.05, addCoef.col ='black', number.cex = 0.8)

