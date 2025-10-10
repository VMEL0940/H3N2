### Influ H3N2 GenDist Proj ###
## Package ##
library(tidyverse)
library(factoextra)
library(FactoMineR)
library(ggpubr)

## Set working directory ##
setwd("/01_GenDisFlu/")

## Data import ##
HA_G1 <- read.table("Data/Genetic Distances/HA distances/HA_G1_distances_2.txt", header = T)

## Pretest - Check Distribution of N value by Tree
HA_G1 %>% 
  filter(vaccineStrain == "EPI103320_A_Moscow_10_1999_NA_NA_NA_NA" & 
           compareStrain == "EPI103320_A_Moscow_10_1999_NA_NA_NA_NA") %>% 
  ggplot(aes(x=as.factor(N))) +
    geom_histogram(stat = "count") +
    theme_bw()
  
## Extract the "Noise AA sub values ##



## N subs ##
Nnse <- as.data.frame(HA_G1$tree.nr)
Nnse$value <- 0
HA_G1$N3 <- 0

for (i in 1:nrow(HA_G1))
  {Nnse$value[i] <- HA_G1 %>% 
     filter(vaccineStrain == compareStrain) %>% 
     filter(tree.nr == HA_G1$tree.nr[i]) %>% 
    .$N
  HA_G1$N3[i] <- ifelse(Nnse$value[i] == 0, HA_G1$N[i], round(HA_G1$N[i] / Nnse$value[i], 0))
  print(HA_G1$N3[i])
}

HA_G1$N4 <- 0

for (i in 1:nrow(HA_G1))
{Nnse$value[i] <- HA_G1 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == HA_G1$tree.nr[i]) %>% 
  .$N
HA_G1$N4[i] <- abs(HA_G1$N[i] - Nnse$value[i])
print(HA_G1$N4[i])
}


HA_G1$N5 <- 0

for (i in 1:nrow(HA_G1))
{Nnse$value[i] <- HA_G1 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == HA_G1$tree.nr[i]) %>% 
  .$N
HA_G1$N5[i] <- ifelse(HA_G1$N[i] - Nnse$value[i] < 0, HA_G1$N[i] + Nnse$value[i], HA_G1$N[i] - Nnse$value[i])
print(HA_G1$N5[i])
}



HA_G1$N6 <- 0

for (i in 1:nrow(HA_G1))
{Nnse$value[i] <- HA_G1 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == HA_G1$tree.nr[i]) %>% 
  .$N
HA_G1$N6[i] <- ifelse(HA_G1$N[i] - Nnse$value[i] < 0, HA_G1$N[i], HA_G1$N[i] - Nnse$value[i])
print(HA_G1$N6[i])
}



Fig_HAG1_N6 <- HA_G1 %>%  
  ggplot(aes(x=as.factor(N6))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_HAG1_N5 <- HA_G1 %>%  
  ggplot(aes(x=as.factor(N5))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_HAG1_N4 <- HA_G1 %>%  
  ggplot(aes(x=as.factor(N4))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_HAG1_N3 <- HA_G1 %>%  
  ggplot(aes(x=as.factor(N3))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_HAG1_N2 <- HA_G1 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_HAG1_N <- HA_G1 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_HAG1_N, Fig_HAG1_N2, Fig_HAG1_N3, Fig_HAG1_N4, Fig_HAG1_N5, Fig_HAG1_N6,
          labels = c("Original N", "Substract", "Divide",
                     "Absolute value", "Relative distance", "Partial substract"),
          nrow = 6)

## Make Absolute values 

## Data import ##
NA_G1 <- read.table("Data/Genetic Distances/NA distances/NA_G1_distances_2.txt", header = T)
NA_G2 <- read.table("Data/Genetic Distances/NA distances/NA_G2_distances_2.txt", header = T)
NA_G3 <- read.table("Data/Genetic Distances/NA distances/NA_G3_distances_2.txt", header = T)
NA_G4 <- read.table("Data/Genetic Distances/NA distances/NA_G4_distances_2.txt", header = T)
NA_G5 <- read.table("Data/Genetic Distances/NA distances/NA_G5_distances_2.txt", header = T)
NA_G6 <- read.table("Data/Genetic Distances/NA distances/NA_G6_distances_2.txt", header = T)
NA_G7 <- read.table("Data/Genetic Distances/NA distances/NA_G7_distances_2.txt", header = T)
NA_G8 <- read.table("Data/Genetic Distances/NA distances/NA_G8_distances_2.txt", header = T)
NA_G9 <- read.table("Data/Genetic Distances/NA distances/NA_G9_distances_2.txt", header = T)
NA_G10 <- read.table("Data/Genetic Distances/NA distances/NA_G10_distances_2.txt", header = T)

## Make Absolute column
NA_G1$N3 <- abs(NA_G1$N2)
NA_G1$S3 <- abs(NA_G1$S2)

## Make Absolute column
NA_G2$N3 <- abs(NA_G2$N2)
NA_G2$S3 <- abs(NA_G2$S2)

## Make Absolute column
NA_G3$N3 <- abs(NA_G3$N2)
NA_G3$S3 <- abs(NA_G3$S2)

## Make Absolute column
NA_G4$N3 <- abs(NA_G4$N2)
NA_G4$S3 <- abs(NA_G4$S2)

## Make Absolute column
NA_G5$N3 <- abs(NA_G5$N2)
NA_G5$S3 <- abs(NA_G5$S2)

## Make Absolute column
NA_G6$N3 <- abs(NA_G6$N2)
NA_G6$S3 <- abs(NA_G6$S2)

## Make Absolute column
NA_G7$N3 <- abs(NA_G7$N2)
NA_G7$S3 <- abs(NA_G7$S2)

## Make Absolute column
NA_G8$N3 <- abs(NA_G8$N2)
NA_G8$S3 <- abs(NA_G8$S2)

## Make Absolute column
NA_G9$N3 <- abs(NA_G9$N2)
NA_G9$S3 <- abs(NA_G9$S2)

## Make Absolute column
NA_G10$N3 <- abs(NA_G10$N2)
NA_G10$S3 <- abs(NA_G10$S2)

## Write Table import ##
write.table(NA_G1,file="Data/Genetic Distances/NA distances/NA_G1_distances_3.txt")
write.table(NA_G2,file="Data/Genetic Distances/NA distances/NA_G2_distances_3.txt")
write.table(NA_G3,file="Data/Genetic Distances/NA distances/NA_G3_distances_3.txt")
write.table(NA_G4,file="Data/Genetic Distances/NA distances/NA_G4_distances_3.txt")
write.table(NA_G5,file="Data/Genetic Distances/NA distances/NA_G5_distances_3.txt")
write.table(NA_G6,file="Data/Genetic Distances/NA distances/NA_G6_distances_3.txt")
write.table(NA_G7,file="Data/Genetic Distances/NA distances/NA_G7_distances_3.txt")
write.table(NA_G8,file="Data/Genetic Distances/NA distances/NA_G8_distances_3.txt")
write.table(NA_G9,file="Data/Genetic Distances/NA distances/NA_G9_distances_3.txt")
write.table(NA_G10,file="Data/Genetic Distances/NA distances/NA_G10_distances_3.txt")




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
N_G10 <- HA_G10[,c(1,4,5,8)]
colnames(N_G10)[4] <- "HA"
N_G10$"NA" <- NA_G10[,8]
N_G10$"NP" <- NP_G10[,8]
N_G10$"M" <- M_G10[,8]
N_G10$"PB1" <- PB1_G10[,8]
N_G10$"PB2" <- PB2_G10[,8]
N_G10$"NS" <- NS_G10[,8]
N_G10$"PA" <- PA_G10[,8]

## Summarize the 901 tree values
N_G10_Ndist <- N_G10 %>% 
  gather("Gene", "N", 4:11) %>% 
  group_by(vaccineStrain, compareStrain, Gene) %>% 
  summarise(Nmean = mean(N), Nmedian = median(N), Nerror = mean(N)-median(N))

########## Realm of Syn ########

## Combine all data with S
S_G10 <- HA_G10[,c(1,4,5,9)]
colnames(S_G10)[4] <- "HA"
S_G10$"NA" <- NA_G10[,9]
S_G10$"NP" <- NP_G10[,9]
S_G10$"M" <- M_G10[,9]
S_G10$"PB1" <- PB1_G10[,9]
S_G10$"PB2" <- PB2_G10[,9]
S_G10$"NS" <- NS_G10[,9]
S_G10$"PA" <- PA_G10[,9]

## Summarize the 901 tree values
S_G10_Ndist <- S_G10 %>% 
  gather("Gene", "S", 4:11) %>% 
  group_by(vaccineStrain, compareStrain, Gene) %>% 
  summarise(Smean = mean(S), Smedian = median(S), Serror = mean(S)-median(S))

## Export the intermediate results
G10_NS <- cbind(N_G10_Ndist, S_G10_Ndist[,4:6])

## Make year column
G10_NS$Year <- as.numeric(str_sub(gsub("_NA", "", G10_NS$compareStrain),-4,-1))

write.table(G10_NS,file="Data/Genetic Distances/GeneCompete/G10_NS_2.txt")



