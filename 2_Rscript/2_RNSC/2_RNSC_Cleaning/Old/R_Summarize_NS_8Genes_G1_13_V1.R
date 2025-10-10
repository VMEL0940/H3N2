### Influ H3N2 GenDist Proj ###
## Package ##
library(tidyverse)

## Set working directory ##
setwd("01_GenDisFlu/")

########### G1 #################

## Data import ##
HA_G1 <- read.table("Data/Genetic Distances/Verification_Final_720/02_summary/04_HA/HA_summary_G1.txt", header = T)
NA_G1 <- read.table("Data/Genetic Distances/Verification_Final_720/02_summary/06_NA/NA_summary_G1.txt", header = T)
NP_G1 <- read.table("Data/Genetic Distances/Verification_Final_720/02_summary/05_NP/NP_summary_G1.txt", header = T)
M_G1 <- read.table("Data/Genetic Distances/Verification_Final_720/02_summary/07_M/M_summary_G1.txt", header = T)
PB1_G1 <- read.table("Data/Genetic Distances/Verification_Final_720/02_summary/02_PB1/PB1_G1_distances.txt", header = T)
PB2_G1 <- read.table("Data/Genetic Distances/Verification_702/PB2/PB2_G1_distances_2.txt", header = T)
PA_G1 <- read.table("Data/Genetic Distances/Verification_702/PA/PA_G1_distances_2.txt", header = T)
NS_G1 <- read.table("Data/Genetic Distances/Verification_702/NS/NS_G1_distances_2.txt", header = T)

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
  summarise(Nmean = mean(abs(N)), Nmedian = median(abs(N)), Nerror = mean(abs(N))-median(abs(N)))

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
  summarise(Smean = mean(abs(S)), Smedian = median(abs(S)), Serror = mean(abs(S))-median(abs(S)))

## Export the intermediate results
G1_NS <- cbind(N_G1_Ndist, S_G1_Ndist[,4:6])

## Make year column
G1_NS$Year <- as.numeric(str_sub(gsub("_NA", "", G1_NS$compareStrain),-4,-1))


########### G2 #################

## Data import ##
HA_G2 <- read.table("Data/Genetic Distances/Verification_702/HA/HA_G2_distances_2.txt", header = T)
NA_G2 <- read.table("Data/Genetic Distances/Verification_702/NA/NA_G2_distances_2.txt", header = T)
NP_G2 <- read.table("Data/Genetic Distances/Verification_702/NP/NP_G2_distances_2.txt", header = T)
M_G2 <- read.table("Data/Genetic Distances/Verification_702/M/M_G2_distances_2.txt", header = T)
PB1_G2 <- read.table("Data/Genetic Distances/Verification_702/PB1/PB1_G2_distances_2.txt", header = T)
PB2_G2 <- read.table("Data/Genetic Distances/Verification_702/PB2/PB2_G2_distances_2.txt", header = T)
PA_G2 <- read.table("Data/Genetic Distances/Verification_702/PA/PA_G2_distances_2.txt", header = T)
NS_G2 <- read.table("Data/Genetic Distances/Verification_702/NS/NS_G2_distances_2.txt", header = T)




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
N_G2 <- HA_G2[,c(1,4:6)]
colnames(N_G2)[4] <- "HA"
N_G2$"NA" <- NA_G2[,6]
N_G2$"NP" <- NP_G2[,6]
N_G2$"M" <- M_G2[,6]
N_G2$"PB1" <- PB1_G2[,6]
N_G2$"PB2" <- PB2_G2[,6]
N_G2$"NS" <- NS_G2[,6]
N_G2$"PA" <- PA_G2[,6]


## Summarize the 901 tree values
N_G2_Ndist <- N_G2 %>% 
  gather("Gene", "N", 4:11) %>% 
  group_by(vaccineStrain, compareStrain, Gene) %>% 
  summarise(Nmean = mean(abs(N)), Nmedian = median(abs(N)), Nerror = mean(abs(N))-median(abs(N)))

########## Realm of Syn ########

## Combine all data with S
S_G2 <- HA_G2[,c(1,4,5,7)]
colnames(S_G2)[4] <- "HA"
S_G2$"NA" <- NA_G2[,7]
S_G2$"NP" <- NP_G2[,7]
S_G2$"M" <- M_G2[,7]
S_G2$"PB1" <- PB1_G2[,7]
S_G2$"PB2" <- PB2_G2[,7]
S_G2$"NS" <- NS_G2[,7]
S_G2$"PA" <- PA_G2[,7]

## Summarize the 901 tree values
S_G2_Ndist <- S_G2 %>% 
  gather("Gene", "S", 4:11) %>% 
  group_by(vaccineStrain, compareStrain, Gene) %>% 
  summarise(Smean = mean(abs(S)), Smedian = median(abs(S)), Serror = mean(abs(S))-median(abs(S)))

## Export the intermediate results
G2_NS <- cbind(N_G2_Ndist, S_G2_Ndist[,4:6])

## Make year column
G2_NS$Year <- as.numeric(str_sub(gsub("_NA", "", G2_NS$compareStrain),-4,-1))


########### G3 #################

## Data import ##
HA_G3 <- read.table("Data/Genetic Distances/Verification_702/HA/HA_G3_distances_2.txt", header = T)
NA_G3 <- read.table("Data/Genetic Distances/Verification_702/NA/NA_G3_distances_2.txt", header = T)
NP_G3 <- read.table("Data/Genetic Distances/Verification_702/NP/NP_G3_distances_2.txt", header = T)
M_G3 <- read.table("Data/Genetic Distances/Verification_702/M/M_G3_distances_2.txt", header = T)
PB1_G3 <- read.table("Data/Genetic Distances/Verification_702/PB1/PB1_G3_distances_2.txt", header = T)
PB2_G3 <- read.table("Data/Genetic Distances/Verification_702/PB2/PB2_G3_distances_2.txt", header = T)
PA_G3 <- read.table("Data/Genetic Distances/Verification_702/PA/PA_G3_distances_2.txt", header = T)
NS_G3 <- read.table("Data/Genetic Distances/Verification_702/NS/NS_G3_distances_2.txt", header = T)




## Data Sort ##
HA_G3 <- HA_G3[order(HA_G3$tree.nr, HA_G3$compareStrain),]
NA_G3 <- NA_G3[order(NA_G3$tree.nr, NA_G3$compareStrain),]
NP_G3 <- NP_G3[order(NP_G3$tree.nr, NP_G3$compareStrain),]
M_G3 <- M_G3[order(M_G3$tree.nr, M_G3$compareStrain),]
PB1_G3 <- PB1_G3[order(PB1_G3$tree.nr, PB1_G3$compareStrain),]
PB2_G3 <- PB2_G3[order(PB2_G3$tree.nr, PB2_G3$compareStrain),]
NS_G3 <- NS_G3[order(NS_G3$tree.nr, NS_G3$compareStrain),]
PA_G3 <- PA_G3[order(PA_G3$tree.nr, PA_G3$compareStrain),]

## Combine all data with N
N_G3 <- HA_G3[,c(1,4:6)]
colnames(N_G3)[4] <- "HA"
N_G3$"NA" <- NA_G3[,6]
N_G3$"NP" <- NP_G3[,6]
N_G3$"M" <- M_G3[,6]
N_G3$"PB1" <- PB1_G3[,6]
N_G3$"PB2" <- PB2_G3[,6]
N_G3$"NS" <- NS_G3[,6]
N_G3$"PA" <- PA_G3[,6]


## Summarize the 901 tree values
N_G3_Ndist <- N_G3 %>% 
  gather("Gene", "N", 4:11) %>% 
  group_by(vaccineStrain, compareStrain, Gene) %>% 
  summarise(Nmean = mean(abs(N)), Nmedian = median(abs(N)), Nerror = mean(abs(N))-median(abs(N)))

########## Realm of Syn ########

## Combine all data with S
S_G3 <- HA_G3[,c(1,4,5,7)]
colnames(S_G3)[4] <- "HA"
S_G3$"NA" <- NA_G3[,7]
S_G3$"NP" <- NP_G3[,7]
S_G3$"M" <- M_G3[,7]
S_G3$"PB1" <- PB1_G3[,7]
S_G3$"PB2" <- PB2_G3[,7]
S_G3$"NS" <- NS_G3[,7]
S_G3$"PA" <- PA_G3[,7]

## Summarize the 901 tree values
S_G3_Ndist <- S_G3 %>% 
  gather("Gene", "S", 4:11) %>% 
  group_by(vaccineStrain, compareStrain, Gene) %>% 
  summarise(Smean = mean(abs(S)), Smedian = median(abs(S)), Serror = mean(abs(S))-median(abs(S)))

## Export the intermediate results
G3_NS <- cbind(N_G3_Ndist, S_G3_Ndist[,4:6])

## Make year column
G3_NS$Year <- as.numeric(str_sub(gsub("_NA", "", G3_NS$compareStrain),-4,-1))


########### G4 #################

## Data import ##
HA_G4 <- read.table("Data/Genetic Distances/Verification_702/HA/HA_G4_distances_2.txt", header = T)
NA_G4 <- read.table("Data/Genetic Distances/Verification_702/NA/NA_G4_distances_2.txt", header = T)
NP_G4 <- read.table("Data/Genetic Distances/Verification_702/NP/NP_G4_distances_2.txt", header = T)
M_G4 <- read.table("Data/Genetic Distances/Verification_702/M/M_G4_distances_2.txt", header = T)
PB1_G4 <- read.table("Data/Genetic Distances/Verification_702/PB1/PB1_G4_distances_2.txt", header = T)
PB2_G4 <- read.table("Data/Genetic Distances/Verification_702/PB2/PB2_G4_distances_2.txt", header = T)
PA_G4 <- read.table("Data/Genetic Distances/Verification_702/PA/PA_G4_distances_2.txt", header = T)
NS_G4 <- read.table("Data/Genetic Distances/Verification_702/NS/NS_G4_distances_2.txt", header = T)




## Data Sort ##
HA_G4 <- HA_G4[order(HA_G4$tree.nr, HA_G4$compareStrain),]
NA_G4 <- NA_G4[order(NA_G4$tree.nr, NA_G4$compareStrain),]
NP_G4 <- NP_G4[order(NP_G4$tree.nr, NP_G4$compareStrain),]
M_G4 <- M_G4[order(M_G4$tree.nr, M_G4$compareStrain),]
PB1_G4 <- PB1_G4[order(PB1_G4$tree.nr, PB1_G4$compareStrain),]
PB2_G4 <- PB2_G4[order(PB2_G4$tree.nr, PB2_G4$compareStrain),]
NS_G4 <- NS_G4[order(NS_G4$tree.nr, NS_G4$compareStrain),]
PA_G4 <- PA_G4[order(PA_G4$tree.nr, PA_G4$compareStrain),]

## Combine all data with N
N_G4 <- HA_G4[,c(1,4:6)]
colnames(N_G4)[4] <- "HA"
N_G4$"NA" <- NA_G4[,6]
N_G4$"NP" <- NP_G4[,6]
N_G4$"M" <- M_G4[,6]
N_G4$"PB1" <- PB1_G4[,6]
N_G4$"PB2" <- PB2_G4[,6]
N_G4$"NS" <- NS_G4[,6]
N_G4$"PA" <- PA_G4[,6]


## Summarize the 901 tree values
N_G4_Ndist <- N_G4 %>% 
  gather("Gene", "N", 4:11) %>% 
  group_by(vaccineStrain, compareStrain, Gene) %>% 
  summarise(Nmean = mean(abs(N)), Nmedian = median(abs(N)), Nerror = mean(abs(N))-median(abs(N)))

########## Realm of Syn ########

## Combine all data with S
S_G4 <- HA_G4[,c(1,4,5,7)]
colnames(S_G4)[4] <- "HA"
S_G4$"NA" <- NA_G4[,7]
S_G4$"NP" <- NP_G4[,7]
S_G4$"M" <- M_G4[,7]
S_G4$"PB1" <- PB1_G4[,7]
S_G4$"PB2" <- PB2_G4[,7]
S_G4$"NS" <- NS_G4[,7]
S_G4$"PA" <- PA_G4[,7]

## Summarize the 901 tree values
S_G4_Ndist <- S_G4 %>% 
  gather("Gene", "S", 4:11) %>% 
  group_by(vaccineStrain, compareStrain, Gene) %>% 
  summarise(Smean = mean(abs(S)), Smedian = median(abs(S)), Serror = mean(abs(S))-median(abs(S)))

## Export the intermediate results
G4_NS <- cbind(N_G4_Ndist, S_G4_Ndist[,4:6])

## Make year column
G4_NS$Year <- as.numeric(str_sub(gsub("_NA", "", G4_NS$compareStrain),-4,-1))


########### G5 #################

## Data import ##
HA_G5 <- read.table("Data/Genetic Distances/Verification_702/HA/HA_G5_distances_2.txt", header = T)
NA_G5 <- read.table("Data/Genetic Distances/Verification_702/NA/NA_G5_distances_2.txt", header = T)
NP_G5 <- read.table("Data/Genetic Distances/Verification_702/NP/NP_G5_distances_2.txt", header = T)
M_G5 <- read.table("Data/Genetic Distances/Verification_702/M/M_G5_distances_2.txt", header = T)
PB1_G5 <- read.table("Data/Genetic Distances/Verification_702/PB1/PB1_G5_distances_2.txt", header = T)
PB2_G5 <- read.table("Data/Genetic Distances/Verification_702/PB2/PB2_G5_distances_2.txt", header = T)
PA_G5 <- read.table("Data/Genetic Distances/Verification_702/PA/PA_G5_distances_2.txt", header = T)
NS_G5 <- read.table("Data/Genetic Distances/Verification_702/NS/NS_G5_distances_2.txt", header = T)




## Data Sort ##
HA_G5 <- HA_G5[order(HA_G5$tree.nr, HA_G5$compareStrain),]
NA_G5 <- NA_G5[order(NA_G5$tree.nr, NA_G5$compareStrain),]
NP_G5 <- NP_G5[order(NP_G5$tree.nr, NP_G5$compareStrain),]
M_G5 <- M_G5[order(M_G5$tree.nr, M_G5$compareStrain),]
PB1_G5 <- PB1_G5[order(PB1_G5$tree.nr, PB1_G5$compareStrain),]
PB2_G5 <- PB2_G5[order(PB2_G5$tree.nr, PB2_G5$compareStrain),]
NS_G5 <- NS_G5[order(NS_G5$tree.nr, NS_G5$compareStrain),]
PA_G5 <- PA_G5[order(PA_G5$tree.nr, PA_G5$compareStrain),]

## Combine all data with N
N_G5 <- HA_G5[,c(1,4:6)]
colnames(N_G5)[4] <- "HA"
N_G5$"NA" <- NA_G5[,6]
N_G5$"NP" <- NP_G5[,6]
N_G5$"M" <- M_G5[,6]
N_G5$"PB1" <- PB1_G5[,6]
N_G5$"PB2" <- PB2_G5[,6]
N_G5$"NS" <- NS_G5[,6]
N_G5$"PA" <- PA_G5[,6]


## Summarize the 901 tree values
N_G5_Ndist <- N_G5 %>% 
  gather("Gene", "N", 4:11) %>% 
  group_by(vaccineStrain, compareStrain, Gene) %>% 
  summarise(Nmean = mean(abs(N)), Nmedian = median(abs(N)), Nerror = mean(abs(N))-median(abs(N)))

########## Realm of Syn ########

## Combine all data with S
S_G5 <- HA_G5[,c(1,4,5,7)]
colnames(S_G5)[4] <- "HA"
S_G5$"NA" <- NA_G5[,7]
S_G5$"NP" <- NP_G5[,7]
S_G5$"M" <- M_G5[,7]
S_G5$"PB1" <- PB1_G5[,7]
S_G5$"PB2" <- PB2_G5[,7]
S_G5$"NS" <- NS_G5[,7]
S_G5$"PA" <- PA_G5[,7]

## Summarize the 901 tree values
S_G5_Ndist <- S_G5 %>% 
  gather("Gene", "S", 4:11) %>% 
  group_by(vaccineStrain, compareStrain, Gene) %>% 
  summarise(Smean = mean(abs(S)), Smedian = median(abs(S)), Serror = mean(abs(S))-median(abs(S)))

## Export the intermediate results
G5_NS <- cbind(N_G5_Ndist, S_G5_Ndist[,4:6])

## Make year column
G5_NS$Year <- as.numeric(str_sub(gsub("_NA", "", G5_NS$compareStrain),-4,-1))


########### G6 #################

## Data import ##
HA_G6 <- read.table("Data/Genetic Distances/Verification_702/HA/HA_G6_distances_2.txt", header = T)
NA_G6 <- read.table("Data/Genetic Distances/Verification_702/NA/NA_G6_distances_2.txt", header = T)
NP_G6 <- read.table("Data/Genetic Distances/Verification_702/NP/NP_G6_distances_2.txt", header = T)
M_G6 <- read.table("Data/Genetic Distances/Verification_702/M/M_G6_distances_2.txt", header = T)
PB1_G6 <- read.table("Data/Genetic Distances/Verification_702/PB1/PB1_G6_distances_2.txt", header = T)
PB2_G6 <- read.table("Data/Genetic Distances/Verification_702/PB2/PB2_G6_distances_2.txt", header = T)
PA_G6 <- read.table("Data/Genetic Distances/Verification_702/PA/PA_G6_distances_2.txt", header = T)
NS_G6 <- read.table("Data/Genetic Distances/Verification_702/NS/NS_G6_distances_2.txt", header = T)




## Data Sort ##
HA_G6 <- HA_G6[order(HA_G6$tree.nr, HA_G6$compareStrain),]
NA_G6 <- NA_G6[order(NA_G6$tree.nr, NA_G6$compareStrain),]
NP_G6 <- NP_G6[order(NP_G6$tree.nr, NP_G6$compareStrain),]
M_G6 <- M_G6[order(M_G6$tree.nr, M_G6$compareStrain),]
PB1_G6 <- PB1_G6[order(PB1_G6$tree.nr, PB1_G6$compareStrain),]
PB2_G6 <- PB2_G6[order(PB2_G6$tree.nr, PB2_G6$compareStrain),]
NS_G6 <- NS_G6[order(NS_G6$tree.nr, NS_G6$compareStrain),]
PA_G6 <- PA_G6[order(PA_G6$tree.nr, PA_G6$compareStrain),]

## Combine all data with N
N_G6 <- HA_G6[,c(1,4:6)]
colnames(N_G6)[4] <- "HA"
N_G6$"NA" <- NA_G6[,6]
N_G6$"NP" <- NP_G6[,6]
N_G6$"M" <- M_G6[,6]
N_G6$"PB1" <- PB1_G6[,6]
N_G6$"PB2" <- PB2_G6[,6]
N_G6$"NS" <- NS_G6[,6]
N_G6$"PA" <- PA_G6[,6]


## Summarize the 901 tree values
N_G6_Ndist <- N_G6 %>% 
  gather("Gene", "N", 4:11) %>% 
  group_by(vaccineStrain, compareStrain, Gene) %>% 
  summarise(Nmean = mean(abs(N)), Nmedian = median(abs(N)), Nerror = mean(abs(N))-median(abs(N)))

########## Realm of Syn ########

## Combine all data with S
S_G6 <- HA_G6[,c(1,4,5,7)]
colnames(S_G6)[4] <- "HA"
S_G6$"NA" <- NA_G6[,7]
S_G6$"NP" <- NP_G6[,7]
S_G6$"M" <- M_G6[,7]
S_G6$"PB1" <- PB1_G6[,7]
S_G6$"PB2" <- PB2_G6[,7]
S_G6$"NS" <- NS_G6[,7]
S_G6$"PA" <- PA_G6[,7]

## Summarize the 901 tree values
S_G6_Ndist <- S_G6 %>% 
  gather("Gene", "S", 4:11) %>% 
  group_by(vaccineStrain, compareStrain, Gene) %>% 
  summarise(Smean = mean(abs(S)), Smedian = median(abs(S)), Serror = mean(abs(S))-median(abs(S)))

## Export the intermediate results
G6_NS <- cbind(N_G6_Ndist, S_G6_Ndist[,4:6])

## Make year column
G6_NS$Year <- as.numeric(str_sub(gsub("_NA", "", G6_NS$compareStrain),-4,-1))


########### G7 #################

## Data import ##
HA_G7 <- read.table("Data/Genetic Distances/Verification_702/HA/HA_G7_distances_2.txt", header = T)
NA_G7 <- read.table("Data/Genetic Distances/Verification_702/NA/NA_G7_distances_2.txt", header = T)
NP_G7 <- read.table("Data/Genetic Distances/Verification_702/NP/NP_G7_distances_2.txt", header = T)
M_G7 <- read.table("Data/Genetic Distances/Verification_702/M/M_G7_distances_2.txt", header = T)
PB1_G7 <- read.table("Data/Genetic Distances/Verification_702/PB1/PB1_G7_distances_2.txt", header = T)
PB2_G7 <- read.table("Data/Genetic Distances/Verification_702/PB2/PB2_G7_distances_2.txt", header = T)
PA_G7 <- read.table("Data/Genetic Distances/Verification_702/PA/PA_G7_distances_2.txt", header = T)
NS_G7 <- read.table("Data/Genetic Distances/Verification_702/NS/NS_G7_distances_2.txt", header = T)




## Data Sort ##
HA_G7 <- HA_G7[order(HA_G7$tree.nr, HA_G7$compareStrain),]
NA_G7 <- NA_G7[order(NA_G7$tree.nr, NA_G7$compareStrain),]
NP_G7 <- NP_G7[order(NP_G7$tree.nr, NP_G7$compareStrain),]
M_G7 <- M_G7[order(M_G7$tree.nr, M_G7$compareStrain),]
PB1_G7 <- PB1_G7[order(PB1_G7$tree.nr, PB1_G7$compareStrain),]
PB2_G7 <- PB2_G7[order(PB2_G7$tree.nr, PB2_G7$compareStrain),]
NS_G7 <- NS_G7[order(NS_G7$tree.nr, NS_G7$compareStrain),]
PA_G7 <- PA_G7[order(PA_G7$tree.nr, PA_G7$compareStrain),]

## Combine all data with N
N_G7 <- HA_G7[,c(1,4:6)]
colnames(N_G7)[4] <- "HA"
N_G7$"NA" <- NA_G7[,6]
N_G7$"NP" <- NP_G7[,6]
N_G7$"M" <- M_G7[,6]
N_G7$"PB1" <- PB1_G7[,6]
N_G7$"PB2" <- PB2_G7[,6]
N_G7$"NS" <- NS_G7[,6]
N_G7$"PA" <- PA_G7[,6]


## Summarize the 901 tree values
N_G7_Ndist <- N_G7 %>% 
  gather("Gene", "N", 4:11) %>% 
  group_by(vaccineStrain, compareStrain, Gene) %>% 
  summarise(Nmean = mean(abs(N)), Nmedian = median(abs(N)), Nerror = mean(abs(N))-median(abs(N)))

########## Realm of Syn ########

## Combine all data with S
S_G7 <- HA_G7[,c(1,4,5,7)]
colnames(S_G7)[4] <- "HA"
S_G7$"NA" <- NA_G7[,7]
S_G7$"NP" <- NP_G7[,7]
S_G7$"M" <- M_G7[,7]
S_G7$"PB1" <- PB1_G7[,7]
S_G7$"PB2" <- PB2_G7[,7]
S_G7$"NS" <- NS_G7[,7]
S_G7$"PA" <- PA_G7[,7]

## Summarize the 901 tree values
S_G7_Ndist <- S_G7 %>% 
  gather("Gene", "S", 4:11) %>% 
  group_by(vaccineStrain, compareStrain, Gene) %>% 
  summarise(Smean = mean(abs(S)), Smedian = median(abs(S)), Serror = mean(abs(S))-median(abs(S)))

## Export the intermediate results
G7_NS <- cbind(N_G7_Ndist, S_G7_Ndist[,4:6])

## Make year column
G7_NS$Year <- as.numeric(str_sub(gsub("_NA", "", G7_NS$compareStrain),-4,-1))


########### G8 #################

## Data import ##
HA_G8 <- read.table("Data/Genetic Distances/Verification_702/HA/HA_G8_distances_2.txt", header = T)
NA_G8 <- read.table("Data/Genetic Distances/Verification_702/NA/NA_G8_distances_2.txt", header = T)
NP_G8 <- read.table("Data/Genetic Distances/Verification_702/NP/NP_G8_distances_2.txt", header = T)
M_G8 <- read.table("Data/Genetic Distances/Verification_702/M/M_G8_distances_2.txt", header = T)
PB1_G8 <- read.table("Data/Genetic Distances/Verification_702/PB1/PB1_G8_distances_2.txt", header = T)
PB2_G8 <- read.table("Data/Genetic Distances/Verification_702/PB2/PB2_G8_distances_2.txt", header = T)
PA_G8 <- read.table("Data/Genetic Distances/Verification_702/PA/PA_G8_distances_2.txt", header = T)
NS_G8 <- read.table("Data/Genetic Distances/Verification_702/NS/NS_G8_distances_2.txt", header = T)




## Data Sort ##
HA_G8 <- HA_G8[order(HA_G8$tree.nr, HA_G8$compareStrain),]
NA_G8 <- NA_G8[order(NA_G8$tree.nr, NA_G8$compareStrain),]
NP_G8 <- NP_G8[order(NP_G8$tree.nr, NP_G8$compareStrain),]
M_G8 <- M_G8[order(M_G8$tree.nr, M_G8$compareStrain),]
PB1_G8 <- PB1_G8[order(PB1_G8$tree.nr, PB1_G8$compareStrain),]
PB2_G8 <- PB2_G8[order(PB2_G8$tree.nr, PB2_G8$compareStrain),]
NS_G8 <- NS_G8[order(NS_G8$tree.nr, NS_G8$compareStrain),]
PA_G8 <- PA_G8[order(PA_G8$tree.nr, PA_G8$compareStrain),]

## Combine all data with N
N_G8 <- HA_G8[,c(1,4:6)]
colnames(N_G8)[4] <- "HA"
N_G8$"NA" <- NA_G8[,6]
N_G8$"NP" <- NP_G8[,6]
N_G8$"M" <- M_G8[,6]
N_G8$"PB1" <- PB1_G8[,6]
N_G8$"PB2" <- PB2_G8[,6]
N_G8$"NS" <- NS_G8[,6]
N_G8$"PA" <- PA_G8[,6]


## Summarize the 901 tree values
N_G8_Ndist <- N_G8 %>% 
  gather("Gene", "N", 4:11) %>% 
  group_by(vaccineStrain, compareStrain, Gene) %>% 
  summarise(Nmean = mean(abs(N)), Nmedian = median(abs(N)), Nerror = mean(abs(N))-median(abs(N)))

########## Realm of Syn ########

## Combine all data with S
S_G8 <- HA_G8[,c(1,4,5,7)]
colnames(S_G8)[4] <- "HA"
S_G8$"NA" <- NA_G8[,7]
S_G8$"NP" <- NP_G8[,7]
S_G8$"M" <- M_G8[,7]
S_G8$"PB1" <- PB1_G8[,7]
S_G8$"PB2" <- PB2_G8[,7]
S_G8$"NS" <- NS_G8[,7]
S_G8$"PA" <- PA_G8[,7]

## Summarize the 901 tree values
S_G8_Ndist <- S_G8 %>% 
  gather("Gene", "S", 4:11) %>% 
  group_by(vaccineStrain, compareStrain, Gene) %>% 
  summarise(Smean = mean(abs(S)), Smedian = median(abs(S)), Serror = mean(abs(S))-median(abs(S)))

## Export the intermediate results
G8_NS <- cbind(N_G8_Ndist, S_G8_Ndist[,4:6])

## Make year column
G8_NS$Year <- as.numeric(str_sub(gsub("_NA", "", G8_NS$compareStrain),-4,-1))


########### G9 #################

## Data import ##
HA_G9 <- read.table("Data/Genetic Distances/Verification_702/HA/HA_G9_distances_2.txt", header = T)
NA_G9 <- read.table("Data/Genetic Distances/Verification_702/NA/NA_G9_distances_2.txt", header = T)
NP_G9 <- read.table("Data/Genetic Distances/Verification_702/NP/NP_G9_distances_2.txt", header = T)
M_G9 <- read.table("Data/Genetic Distances/Verification_702/M/M_G9_distances_2.txt", header = T)
PB1_G9 <- read.table("Data/Genetic Distances/Verification_702/PB1/PB1_G9_distances_2.txt", header = T)
PB2_G9 <- read.table("Data/Genetic Distances/Verification_702/PB2/PB2_G9_distances_2.txt", header = T)
PA_G9 <- read.table("Data/Genetic Distances/Verification_702/PA/PA_G9_distances_2.txt", header = T)
NS_G9 <- read.table("Data/Genetic Distances/Verification_702/NS/NS_G9_distances_2.txt", header = T)




## Data Sort ##
HA_G9 <- HA_G9[order(HA_G9$tree.nr, HA_G9$compareStrain),]
NA_G9 <- NA_G9[order(NA_G9$tree.nr, NA_G9$compareStrain),]
NP_G9 <- NP_G9[order(NP_G9$tree.nr, NP_G9$compareStrain),]
M_G9 <- M_G9[order(M_G9$tree.nr, M_G9$compareStrain),]
PB1_G9 <- PB1_G9[order(PB1_G9$tree.nr, PB1_G9$compareStrain),]
PB2_G9 <- PB2_G9[order(PB2_G9$tree.nr, PB2_G9$compareStrain),]
NS_G9 <- NS_G9[order(NS_G9$tree.nr, NS_G9$compareStrain),]
PA_G9 <- PA_G9[order(PA_G9$tree.nr, PA_G9$compareStrain),]

## Combine all data with N
N_G9 <- HA_G9[,c(1,4:6)]
colnames(N_G9)[4] <- "HA"
N_G9$"NA" <- NA_G9[,6]
N_G9$"NP" <- NP_G9[,6]
N_G9$"M" <- M_G9[,6]
N_G9$"PB1" <- PB1_G9[,6]
N_G9$"PB2" <- PB2_G9[,6]
N_G9$"NS" <- NS_G9[,6]
N_G9$"PA" <- PA_G9[,6]


## Summarize the 901 tree values
N_G9_Ndist <- N_G9 %>% 
  gather("Gene", "N", 4:11) %>% 
  group_by(vaccineStrain, compareStrain, Gene) %>% 
  summarise(Nmean = mean(abs(N)), Nmedian = median(abs(N)), Nerror = mean(abs(N))-median(abs(N)))

########## Realm of Syn ########

## Combine all data with S
S_G9 <- HA_G9[,c(1,4,5,7)]
colnames(S_G9)[4] <- "HA"
S_G9$"NA" <- NA_G9[,7]
S_G9$"NP" <- NP_G9[,7]
S_G9$"M" <- M_G9[,7]
S_G9$"PB1" <- PB1_G9[,7]
S_G9$"PB2" <- PB2_G9[,7]
S_G9$"NS" <- NS_G9[,7]
S_G9$"PA" <- PA_G9[,7]

## Summarize the 901 tree values
S_G9_Ndist <- S_G9 %>% 
  gather("Gene", "S", 4:11) %>% 
  group_by(vaccineStrain, compareStrain, Gene) %>% 
  summarise(Smean = mean(abs(S)), Smedian = median(abs(S)), Serror = mean(abs(S))-median(abs(S)))

## Export the intermediate results
G9_NS <- cbind(N_G9_Ndist, S_G9_Ndist[,4:6])

## Make year column
G9_NS$Year <- as.numeric(str_sub(gsub("_NA", "", G9_NS$compareStrain),-4,-1))


########### G10 #################

## Data import ##
HA_G10 <- read.table("Data/Genetic Distances/Verification_702/HA/HA_G10_distances_2.txt", header = T)
NA_G10 <- read.table("Data/Genetic Distances/Verification_702/NA/NA_G10_distances_2.txt", header = T)
NP_G10 <- read.table("Data/Genetic Distances/Verification_702/NP/NP_G10_distances_2.txt", header = T)
M_G10 <- read.table("Data/Genetic Distances/Verification_702/M/M_G10_distances_2.txt", header = T)
PB1_G10 <- read.table("Data/Genetic Distances/Verification_702/PB1/PB1_G10_distances_2.txt", header = T)
PB2_G10 <- read.table("Data/Genetic Distances/Verification_702/PB2/PB2_G10_distances_2.txt", header = T)
PA_G10 <- read.table("Data/Genetic Distances/Verification_702/PA/PA_G10_distances_2.txt", header = T)
NS_G10 <- read.table("Data/Genetic Distances/Verification_702/NS/NS_G10_distances_2.txt", header = T)




## Data Sort ##
HA_G10 <- HA_G10[order(HA_G10$tree.nr, HA_G10$compareStrain),]
NA_G10 <- NA_G10[order(NA_G10$tree.nr, NA_G10$compareStrain),]
NP_G10 <- NP_G10[order(NP_G10$tree.nr, NP_G10$compareStrain),]
M_G10 <- M_G10[order(M_G10$tree.nr, M_G10$compareStrain),]
PB1_G10 <- PB1_G10[order(PB1_G10$tree.nr, PB1_G10$compareStrain),]
PB2_G10 <- PB2_G10[order(PB2_G10$tree.nr, PB2_G10$compareStrain),]
NS_G10 <- NS_G10[order(NS_G10$tree.nr, NS_G10$compareStrain),]
PA_G10 <- PA_G10[order(PA_G10$tree.nr, PA_G10$compareStrain),]

## Combine all data with N
N_G10 <- HA_G10[,c(1,4:6)]
colnames(N_G10)[4] <- "HA"
N_G10$"NA" <- NA_G10[,6]
N_G10$"NP" <- NP_G10[,6]
N_G10$"M" <- M_G10[,6]
N_G10$"PB1" <- PB1_G10[,6]
N_G10$"PB2" <- PB2_G10[,6]
N_G10$"NS" <- NS_G10[,6]
N_G10$"PA" <- PA_G10[,6]


## Summarize the 901 tree values
N_G10_Ndist <- N_G10 %>% 
  gather("Gene", "N", 4:11) %>% 
  group_by(vaccineStrain, compareStrain, Gene) %>% 
  summarise(Nmean = mean(abs(N)), Nmedian = median(abs(N)), Nerror = mean(abs(N))-median(abs(N)))

########## Realm of Syn ########

## Combine all data with S
S_G10 <- HA_G10[,c(1,4,5,7)]
colnames(S_G10)[4] <- "HA"
S_G10$"NA" <- NA_G10[,7]
S_G10$"NP" <- NP_G10[,7]
S_G10$"M" <- M_G10[,7]
S_G10$"PB1" <- PB1_G10[,7]
S_G10$"PB2" <- PB2_G10[,7]
S_G10$"NS" <- NS_G10[,7]
S_G10$"PA" <- PA_G10[,7]

## Summarize the 901 tree values
S_G10_Ndist <- S_G10 %>% 
  gather("Gene", "S", 4:11) %>% 
  group_by(vaccineStrain, compareStrain, Gene) %>% 
  summarise(Smean = mean(abs(S)), Smedian = median(abs(S)), Serror = mean(abs(S))-median(abs(S)))

## Export the intermediate results
G10_NS <- cbind(N_G10_Ndist, S_G10_Ndist[,4:6])

## Make year column
G10_NS$Year <- as.numeric(str_sub(gsub("_NA", "", G10_NS$compareStrain),-4,-1))


########### G11 #################

## Data import ##
HA_G11 <- read.table("Data/Genetic Distances/Verification_702/HA/HA_G11_distances_2.txt", header = T)
NA_G11 <- read.table("Data/Genetic Distances/Verification_702/NA/NA_G11_distances_2.txt", header = T)
NP_G11 <- read.table("Data/Genetic Distances/Verification_702/NP/NP_G11_distances_2.txt", header = T)
M_G11 <- read.table("Data/Genetic Distances/Verification_702/M/M_G11_distances_2.txt", header = T)
PB1_G11 <- read.table("Data/Genetic Distances/Verification_702/PB1/PB1_G11_distances_2.txt", header = T)
PB2_G11 <- read.table("Data/Genetic Distances/Verification_702/PB2/PB2_G11_distances_2.txt", header = T)
PA_G11 <- read.table("Data/Genetic Distances/Verification_702/PA/PA_G11_distances_2.txt", header = T)
NS_G11 <- read.table("Data/Genetic Distances/Verification_702/NS/NS_G11_distances_2.txt", header = T)




## Data Sort ##
HA_G11 <- HA_G11[order(HA_G11$tree.nr, HA_G11$compareStrain),]
NA_G11 <- NA_G11[order(NA_G11$tree.nr, NA_G11$compareStrain),]
NP_G11 <- NP_G11[order(NP_G11$tree.nr, NP_G11$compareStrain),]
M_G11 <- M_G11[order(M_G11$tree.nr, M_G11$compareStrain),]
PB1_G11 <- PB1_G11[order(PB1_G11$tree.nr, PB1_G11$compareStrain),]
PB2_G11 <- PB2_G11[order(PB2_G11$tree.nr, PB2_G11$compareStrain),]
NS_G11 <- NS_G11[order(NS_G11$tree.nr, NS_G11$compareStrain),]
PA_G11 <- PA_G11[order(PA_G11$tree.nr, PA_G11$compareStrain),]

## Combine all data with N
N_G11 <- HA_G11[,c(1,4:6)]
colnames(N_G11)[4] <- "HA"
N_G11$"NA" <- NA_G11[,6]
N_G11$"NP" <- NP_G11[,6]
N_G11$"M" <- M_G11[,6]
N_G11$"PB1" <- PB1_G11[,6]
N_G11$"PB2" <- PB2_G11[,6]
N_G11$"NS" <- NS_G11[,6]
N_G11$"PA" <- PA_G11[,6]


## Summarize the 901 tree values
N_G11_Ndist <- N_G11 %>% 
  gather("Gene", "N", 4:11) %>% 
  group_by(vaccineStrain, compareStrain, Gene) %>% 
  summarise(Nmean = mean(abs(N)), Nmedian = median(abs(N)), Nerror = mean(abs(N))-median(abs(N)))

########## Realm of Syn ########

## Combine all data with S
S_G11 <- HA_G11[,c(1,4,5,7)]
colnames(S_G11)[4] <- "HA"
S_G11$"NA" <- NA_G11[,7]
S_G11$"NP" <- NP_G11[,7]
S_G11$"M" <- M_G11[,7]
S_G11$"PB1" <- PB1_G11[,7]
S_G11$"PB2" <- PB2_G11[,7]
S_G11$"NS" <- NS_G11[,7]
S_G11$"PA" <- PA_G11[,7]

## Summarize the 901 tree values
S_G11_Ndist <- S_G11 %>% 
  gather("Gene", "S", 4:11) %>% 
  group_by(vaccineStrain, compareStrain, Gene) %>% 
  summarise(Smean = mean(abs(S)), Smedian = median(abs(S)), Serror = mean(abs(S))-median(abs(S)))

## Export the intermediate results
G11_NS <- cbind(N_G11_Ndist, S_G11_Ndist[,4:6])

## Make year column
G11_NS$Year <- as.numeric(str_sub(gsub("_NA", "", G11_NS$compareStrain),-4,-1))


########### G12 #################

## Data import ##
HA_G12 <- read.table("Data/Genetic Distances/Verification_702/HA/HA_G12_distances_2.txt", header = T)
NA_G12 <- read.table("Data/Genetic Distances/Verification_702/NA/NA_G12_distances_2.txt", header = T)
NP_G12 <- read.table("Data/Genetic Distances/Verification_702/NP/NP_G12_distances_2.txt", header = T)
M_G12 <- read.table("Data/Genetic Distances/Verification_702/M/M_G12_distances_2.txt", header = T)
PB1_G12 <- read.table("Data/Genetic Distances/Verification_702/PB1/PB1_G12_distances_2.txt", header = T)
PB2_G12 <- read.table("Data/Genetic Distances/Verification_702/PB2/PB2_G12_distances_2.txt", header = T)
PA_G12 <- read.table("Data/Genetic Distances/Verification_702/PA/PA_G12_distances_2.txt", header = T)
NS_G12 <- read.table("Data/Genetic Distances/Verification_702/NS/NS_G12_distances_2.txt", header = T)




## Data Sort ##
HA_G12 <- HA_G12[order(HA_G12$tree.nr, HA_G12$compareStrain),]
NA_G12 <- NA_G12[order(NA_G12$tree.nr, NA_G12$compareStrain),]
NP_G12 <- NP_G12[order(NP_G12$tree.nr, NP_G12$compareStrain),]
M_G12 <- M_G12[order(M_G12$tree.nr, M_G12$compareStrain),]
PB1_G12 <- PB1_G12[order(PB1_G12$tree.nr, PB1_G12$compareStrain),]
PB2_G12 <- PB2_G12[order(PB2_G12$tree.nr, PB2_G12$compareStrain),]
NS_G12 <- NS_G12[order(NS_G12$tree.nr, NS_G12$compareStrain),]
PA_G12 <- PA_G12[order(PA_G12$tree.nr, PA_G12$compareStrain),]

## Combine all data with N
N_G12 <- HA_G12[,c(1,4:6)]
colnames(N_G12)[4] <- "HA"
N_G12$"NA" <- NA_G12[,6]
N_G12$"NP" <- NP_G12[,6]
N_G12$"M" <- M_G12[,6]
N_G12$"PB1" <- PB1_G12[,6]
N_G12$"PB2" <- PB2_G12[,6]
N_G12$"NS" <- NS_G12[,6]
N_G12$"PA" <- PA_G12[,6]


## Summarize the 901 tree values
N_G12_Ndist <- N_G12 %>% 
  gather("Gene", "N", 4:11) %>% 
  group_by(vaccineStrain, compareStrain, Gene) %>% 
  summarise(Nmean = mean(abs(N)), Nmedian = median(abs(N)), Nerror = mean(abs(N))-median(abs(N)))

########## Realm of Syn ########

## Combine all data with S
S_G12 <- HA_G12[,c(1,4,5,7)]
colnames(S_G12)[4] <- "HA"
S_G12$"NA" <- NA_G12[,7]
S_G12$"NP" <- NP_G12[,7]
S_G12$"M" <- M_G12[,7]
S_G12$"PB1" <- PB1_G12[,7]
S_G12$"PB2" <- PB2_G12[,7]
S_G12$"NS" <- NS_G12[,7]
S_G12$"PA" <- PA_G12[,7]

## Summarize the 901 tree values
S_G12_Ndist <- S_G12 %>% 
  gather("Gene", "S", 4:11) %>% 
  group_by(vaccineStrain, compareStrain, Gene) %>% 
  summarise(Smean = mean(abs(S)), Smedian = median(abs(S)), Serror = mean(abs(S))-median(abs(S)))

## Export the intermediate results
G12_NS <- cbind(N_G12_Ndist, S_G12_Ndist[,4:6])

## Make year column
G12_NS$Year <- as.numeric(str_sub(gsub("_NA", "", G12_NS$compareStrain),-4,-1))




########### G13 #################

## Data import ##
HA_G13 <- read.table("Data/Genetic Distances/Verification_702/HA/HA_G13_distances_2.txt", header = T)
NA_G13 <- read.table("Data/Genetic Distances/Verification_702/NA/NA_G13_distances_2.txt", header = T)
NP_G13 <- read.table("Data/Genetic Distances/Verification_702/NP/NP_G13_distances_2.txt", header = T)
M_G13 <- read.table("Data/Genetic Distances/Verification_702/M/M_G13_distances_2.txt", header = T)
PB1_G13 <- read.table("Data/Genetic Distances/Verification_702/PB1/PB1_G13_distances_2.txt", header = T)
PB2_G13 <- read.table("Data/Genetic Distances/Verification_702/PB2/PB2_G13_distances_2.txt", header = T)
PA_G13 <- read.table("Data/Genetic Distances/Verification_702/PA/PA_G13_distances_2.txt", header = T)
NS_G13 <- read.table("Data/Genetic Distances/Verification_702/NS/NS_G13_distances_2.txt", header = T)




## Data Sort ##
HA_G13 <- HA_G13[order(HA_G13$tree.nr, HA_G13$compareStrain),]
NA_G13 <- NA_G13[order(NA_G13$tree.nr, NA_G13$compareStrain),]
NP_G13 <- NP_G13[order(NP_G13$tree.nr, NP_G13$compareStrain),]
M_G13 <- M_G13[order(M_G13$tree.nr, M_G13$compareStrain),]
PB1_G13 <- PB1_G13[order(PB1_G13$tree.nr, PB1_G13$compareStrain),]
PB2_G13 <- PB2_G13[order(PB2_G13$tree.nr, PB2_G13$compareStrain),]
NS_G13 <- NS_G13[order(NS_G13$tree.nr, NS_G13$compareStrain),]
PA_G13 <- PA_G13[order(PA_G13$tree.nr, PA_G13$compareStrain),]

## Combine all data with N
N_G13 <- HA_G13[,c(1,4:6)]
colnames(N_G13)[4] <- "HA"
N_G13$"NA" <- NA_G13[,6]
N_G13$"NP" <- NP_G13[,6]
N_G13$"M" <- M_G13[,6]
N_G13$"PB1" <- PB1_G13[,6]
N_G13$"PB2" <- PB2_G13[,6]
N_G13$"NS" <- NS_G13[,6]
N_G13$"PA" <- PA_G13[,6]


## Summarize the 901 tree values
N_G13_Ndist <- N_G13 %>% 
  gather("Gene", "N", 4:11) %>% 
  group_by(vaccineStrain, compareStrain, Gene) %>% 
  summarise(Nmean = mean(abs(N)), Nmedian = median(abs(N)), Nerror = mean(abs(N))-median(abs(N)))

########## Realm of Syn ########

## Combine all data with S
S_G13 <- HA_G13[,c(1,4,5,7)]
colnames(S_G13)[4] <- "HA"
S_G13$"NA" <- NA_G13[,7]
S_G13$"NP" <- NP_G13[,7]
S_G13$"M" <- M_G13[,7]
S_G13$"PB1" <- PB1_G13[,7]
S_G13$"PB2" <- PB2_G13[,7]
S_G13$"NS" <- NS_G13[,7]
S_G13$"PA" <- PA_G13[,7]

## Summarize the 901 tree values
S_G13_Ndist <- S_G13 %>% 
  gather("Gene", "S", 4:11) %>% 
  group_by(vaccineStrain, compareStrain, Gene) %>% 
  summarise(Smean = mean(abs(S)), Smedian = median(abs(S)), Serror = mean(abs(S))-median(abs(S)))

## Export the intermediate results
G13_NS <- cbind(N_G13_Ndist, S_G13_Ndist[,4:6])

## Make year column
G13_NS$Year <- as.numeric(str_sub(gsub("_NA", "", G13_NS$compareStrain),-4,-1))

## Make Group column
G13_NS$Group <- "G13"




######### Summarize ##########

## Make Group column
G1_NS$Group <- "G1"
G2_NS$Group <- "G2"
G3_NS$Group <- "G3"
G4_NS$Group <- "G4"
G5_NS$Group <- "G5"
G6_NS$Group <- "G6"
G7_NS$Group <- "G7"
G8_NS$Group <- "G8"
G9_NS$Group <- "G9"
G10_NS$Group <- "G10"
G11_NS$Group <- "G11"
G12_NS$Group <- "G12"
G13_NS$Group <- "G13"

#rbind
AllG_NS<- rbind(G1_NS, G2_NS, G3_NS, G4_NS, G5_NS, G6_NS, G7_NS,
                G8_NS, G9_NS, G10_NS, G11_NS, G12_NS, G13_NS)

write.csv(AllG_NS,file="Data/Genetic Distances/Verification_702/Summary/AllG_NS.csv")


