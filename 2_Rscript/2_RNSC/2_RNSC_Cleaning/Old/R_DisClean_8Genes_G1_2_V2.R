### Influ H3N2 GenDist Proj ###
## Package ##
library(tidyverse)

## Set working directory ##
setwd("/01_GenDisFlu/")

########### HA !!!!!! #################

## Data import ##
HA_G1 <- read.table("Data/Genetic Distances/HA/HA_G1_distances_2.txt", header = T)
NA_G1 <- read.table("Data/Genetic Distances/NA/NA_G1_distances_2.txt", header = T)
NP_G1 <- read.table("Data/Genetic Distances/NP/NP_G1_distances_2.txt", header = T)
M_G1 <- read.table("Data/Genetic Distances/M/M_G1_distances_2.txt", header = T)
PB1_G1 <- read.table("Data/Genetic Distances/PB1/PB1_G1_distances_2.txt", header = T)
PB2_G1 <- read.table("Data/Genetic Distances/PB2/PB2_G1_distances_2.txt", header = T)
PA_G1 <- read.table("Data/Genetic Distances/PA/PA_G1_distances_2.txt", header = T)
NS_G1 <- read.table("Data/Genetic Distances/NS/NS_G1_distances_2.txt", header = T)

## Data Sort ##
HA_G1 <- HA_G1[order(HA_G1$tree.nr, HA_G1$compareStrain),]
NA_G1 <- NA_G1[order(NA_G1$tree.nr, NA_G1$compareStrain),]
NP_G1 <- NP_G1[order(NP_G1$tree.nr, NP_G1$compareStrain),]
M_G1 <- M_G1[order(M_G1$tree.nr, M_G1$compareStrain),]
PB1_G1 <- PB1_G1[order(PB1_G1$tree.nr, PB1_G1$compareStrain),]
PB2_G1 <- PB2_G1[order(PB2_G1$tree.nr, PB2_G1$compareStrain),]
NS_G1 <- NS_G1[order(NS_G1$tree.nr, NS_G1$compareStrain),]
PA_G1 <- PA_G1[order(PA_G1$tree.nr, PA_G1$compareStrain),]

## Combine all data with N
N_G1 <- HA_G1[,c(1,4:6)]
colnames(N_G1)[4] <- "HA"
N_G1$"NA" <- NA_G1[,6]
N_G1$"NP" <- NP_G1[,6]
N_G1$"M" <- M_G1[,6]
N_G1$"PB1" <- PB1_G1[,6]
N_G1$"PB2" <- PB2_G1[,6]
N_G1$"NS" <- NS_G1[,6]
N_G1$"PA" <- PA_G1[,6]


## Summarize the 901 tree values
N_G1_Ndist <- N_G1 %>% 
  gather("Gene", "N", 4:11) %>% 
  group_by(vaccineStrain, compareStrain, Gene) %>% 
  summarise(Nmean = mean(N), Nmedian = median(N), Nerror = mean(N)-median(N))

## Combine all data with S
S_G1 <- HA_G1[,c(1,4,5,7)]
colnames(S_G1)[4] <- "HA"
S_G1$"NA" <- NA_G1[,7]
S_G1$"NP" <- NP_G1[,7]
S_G1$"M" <- M_G1[,7]
S_G1$"PB1" <- PB1_G1[,7]
S_G1$"PB2" <- PB2_G1[,7]
S_G1$"NS" <- NS_G1[,7]
S_G1$"PA" <- PA_G1[,7]

## Summarize the 901 tree values
S_G1_Ndist <- S_G1 %>% 
  gather("Gene", "S", 4:11) %>% 
  group_by(vaccineStrain, compareStrain, Gene) %>% 
  summarise(Smean = mean(S), Smedian = median(S), Serror = mean(S)-median(S))

## Export the intermediate results
G1_NS <- cbind(N_G1_Ndist, S_G1_Ndist[,4:6])

## Make year column
G1_NS$Year <- as.numeric(str_sub(gsub("_NA", "", G1_NS$compareStrain),-4,-1))

write.table(G1_NS,file="Data/Genetic Distances/GeneCompete/G1_NS_2.txt")




## Data Sort ##
HA_G2 <- HA_G2[order(HA_G2$tree.nr, HA_G2$compareStrain),]
NA_G2 <- NA_G2[order(NA_G2$tree.nr, NA_G2$compareStrain),]
NP_G2 <- NP_G2[order(NP_G2$tree.nr, NP_G2$compareStrain),]
M_G2 <- M_G2[order(M_G2$tree.nr, M_G2$compareStrain),]
PB1_G2 <- PB1_G2[order(PB1_G2$tree.nr, PB1_G2$compareStrain),]
PB2_G2 <- PB2_G2[order(PB2_G2$tree.nr, PB2_G2$compareStrain),]
NS_G2 <- NS_G2[order(NS_G2$tree.nr, NS_G2$compareStrain),]
PA_G2 <- PA_G2[order(PA_G2$tree.nr, PA_G2$compareStrain),]

## Combine all data with N
N_G2 <- HA_G2[,c(1,4,5,8)]
colnames(N_G2)[4] <- "HA"
N_G2$"NA" <- NA_G2[,8]
N_G2$"NP" <- NP_G2[,8]
N_G2$"M" <- M_G2[,8]
N_G2$"PB1" <- PB1_G2[,8]
N_G2$"PB2" <- PB2_G2[,8]
N_G2$"NS" <- NS_G2[,8]
N_G2$"PA" <- PA_G2[,8]


## Summarize the 901 tree values
N_G2_Ndist <- N_G2 %>% 
  gather("Gene", "N", 4:11) %>% 
  group_by(vaccineStrain, compareStrain, Gene) %>% 
  summarise(Nmean = mean(N), Nmedian = median(N), Nerror = mean(N)-median(N))

########## Realm of Syn ########

## Combine all data with S
S_G2 <- HA_G2[,c(1,4,5,9)]
colnames(S_G2)[4] <- "HA"
S_G2$"NA" <- NA_G2[,9]
S_G2$"NP" <- NP_G2[,9]
S_G2$"M" <- M_G2[,9]
S_G2$"PB1" <- PB1_G2[,9]
S_G2$"PB2" <- PB2_G2[,9]
S_G2$"NS" <- NS_G2[,9]
S_G2$"PA" <- PA_G2[,9]

## Summarize the 901 tree values
S_G2_Ndist <- S_G2 %>% 
  gather("Gene", "S", 4:11) %>% 
  group_by(vaccineStrain, compareStrain, Gene) %>% 
  summarise(Smean = mean(S), Smedian = median(S), Serror = mean(S)-median(S))

## Export the intermediate results
G2_NS <- cbind(N_G2_Ndist, S_G2_Ndist[,4:6])

## Make year column
G2_NS$Year <- as.numeric(str_sub(gsub("_NA", "", G2_NS$compareStrain),-4,-1))

write.table(G2_NS,file="Data/Genetic Distances/GeneCompete/G2_NS_2.txt")



