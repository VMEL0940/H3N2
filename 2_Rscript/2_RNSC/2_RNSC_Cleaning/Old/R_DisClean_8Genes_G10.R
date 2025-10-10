### Influ H3N2 GenDist Proj ###
## Package ##
library(tidyverse)
library(factoextra)
library(FactoMineR)
library(ggpubr)

## Set working directory ##
setwd("/1_GenDisFlu/")

########### HA !!!!!! #################

## Data import ##
HA_G10 <- read.table("Data/Genetic Distances/HA distances/HA_G10_distances.txt", header = T)
## Build redudant column ##
HA_G10 <- HA_G10[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(HA_G10$tree.nr)
Nnse$value <- 0
HA_G10$N2 <- 0

for (i in 1:nrow(HA_G10))
  {Nnse$value[i] <- HA_G10 %>% 
     filter(vaccineStrain == compareStrain) %>% 
     filter(tree.nr == HA_G10$tree.nr[i]) %>% 
    .$N
  HA_G10$N2[i] <- HA_G10$N[i] - Nnse$value[i]
  print(HA_G10$N2[i])
}

Fig_HAG10_N2 <- HA_G10 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_HAG10_N <- HA_G10 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_HAG10_N2, Fig_HAG10_N, 
          labels = c("N2", "N"),
          nrow = 2)


## S subs ##

HA_G10$S2 <- 0

for (i in 1:nrow(HA_G10))
{Nnse$value[i] <- HA_G10 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == HA_G10$tree.nr[i]) %>% 
  .$S
HA_G10$S2[i] <- HA_G10$S[i] - Nnse$value[i]
print(HA_G10$S2[i])
}

Fig_HAG10_S2 <- HA_G10 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_HAG10_S <- HA_G10 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_HAG10_S2, Fig_HAG10_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(HA_G10,file="Data/Genetic Distances/HA distances/HA_G10_distances_2.txt")


########### NA !!!!!! #################

## Data import ##
NA_G10 <- read.table("Data/Genetic Distances/NA distances/NA_G10_distances.txt", header = T)
## Build redudant column ##
NA_G10 <- NA_G10[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NA_G10$tree.nr)
Nnse$value <- 0
NA_G10$N2 <- 0

for (i in 1:nrow(NA_G10))
{Nnse$value[i] <- NA_G10 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NA_G10$tree.nr[i]) %>% 
  .$N
NA_G10$N2[i] <- NA_G10$N[i] - Nnse$value[i]
print(NA_G10$N2[i])
}

Fig_NAG10_N2 <- NA_G10 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_NAG10_N <- NA_G10 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_NAG10_N2, Fig_NAG10_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

NA_G10$S2 <- 0

for (i in 1:nrow(NA_G10))
{Nnse$value[i] <- NA_G10 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NA_G10$tree.nr[i]) %>% 
  .$S
NA_G10$S2[i] <- NA_G10$S[i] - Nnse$value[i]
print(NA_G10$S2[i])
}

Fig_NAG10_S2 <- NA_G10 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_NAG10_S <- NA_G10 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_NAG10_S2, Fig_NAG10_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(NA_G10,file="Data/Genetic Distances/NA distances/NA_G10_distances_2.txt")


########### M!!!!!! #################

## Data import ##
M_G10 <- read.table("Data/Genetic Distances/M distances/M_G10_distances.txt", header = T)
## Build redudant column ##
M_G10 <- M_G10[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(M_G10$tree.nr)
Nnse$value <- 0
M_G10$N2 <- 0

for (i in 1:nrow(M_G10))
{Nnse$value[i] <- M_G10 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == M_G10$tree.nr[i]) %>% 
  .$N
M_G10$N2[i] <- M_G10$N[i] - Nnse$value[i]
print(M_G10$N2[i])
}

Fig_MG10_N2 <- M_G10 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_MG10_N <- M_G10 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_MG10_N2, Fig_MG10_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

M_G10$S2 <- 0

for (i in 1:nrow(M_G10))
{Nnse$value[i] <- M_G10 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == M_G10$tree.nr[i]) %>% 
  .$S
M_G10$S2[i] <- M_G10$S[i] - Nnse$value[i]
print(M_G10$S2[i])
}

Fig_MG10_S2 <- M_G10 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_MG10_S <- M_G10 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_MG10_S2, Fig_MG10_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(M_G10,file="Data/Genetic Distances/M distances/M_G10_distances_2.txt")


########### NP!!!!!! #################

## Data import ##
NP_G10 <- read.table("Data/Genetic Distances/NP distances/NP_G10_distances.txt", header = T)
## Build redudant column ##
NP_G10 <- NP_G10[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NP_G10$tree.nr)
Nnse$value <- 0
NP_G10$N2 <- 0

for (i in 1:nrow(NP_G10))
{Nnse$value[i] <- NP_G10 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NP_G10$tree.nr[i]) %>% 
  .$N
NP_G10$N2[i] <- NP_G10$N[i] - Nnse$value[i]
print(NP_G10$N2[i])
}

Fig_NPG10_N2 <- NP_G10 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_NPG10_N <- NP_G10 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_NPG10_N2, Fig_NPG10_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

NP_G10$S2 <- 0

for (i in 1:nrow(NP_G10))
{Nnse$value[i] <- NP_G10 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NP_G10$tree.nr[i]) %>% 
  .$S
NP_G10$S2[i] <- NP_G10$S[i] - Nnse$value[i]
print(NP_G10$S2[i])
}

Fig_NPG10_S2 <- NP_G10 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_NPG10_S <- NP_G10 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_NPG10_S2, Fig_NPG10_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(NP_G10,file="Data/Genetic Distances/NP distances/NP_G10_distances_2.txt", row.names = F)


########### NS!!!!!! #################


## Data import ##
NS_G10 <- read.table("Data/Genetic Distances/NS distances/NS_G10_distances.txt", header = T)
## Build redudant column ##
NS_G10 <- NS_G10[,1:5]


## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NS_G10$tree.nr)
Nnse$value <- 0
NS_G10$N2 <- NS_G10$N

for (i in 1:nrow(NS_G10))
  {Nnse$value[i] <- NS_G10 %>% 
     filter(vaccineStrain == compareStrain) %>% 
     filter(tree.nr == NS_G10$tree.nr[i]) %>% 
    .$N
  NS_G10$N2[i] <- NS_G10$N[i] - Nnse$value[i]
  print(NS_G10$N2[i])
}

Fig_NSG10_N2 <- NS_G10 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_NSG10_N <- NS_G10 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_NSG10_N2, Fig_NSG10_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

NS_G10$S2 <- NS_G10$S

for (i in 1:nrow(NS_G10))
{Nnse$value[i] <- NS_G10 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NS_G10$tree.nr[i]) %>% 
  .$S
NS_G10$S2[i] <- NS_G10$S[i] - Nnse$value[i]
print(NS_G10$S2[i])
}

Fig_NSG10_S2 <- NS_G10 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_NSG10_S <- NS_G10 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_NSG10_S2, Fig_NSG10_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(NS_G10,file="Data/Genetic Distances/NS distances/NS_G10_distances_2.txt")


########### PA!!!!!! #################

## Data import ##
PA_G10 <- read.table("Data/Genetic Distances/PA distances/PA_G10_distances.txt", header = T)
## Build redudant column ##
PA_G10 <- PA_G10[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(PA_G10$tree.nr)
Nnse$value <- 0
PA_G10$N2 <- 0

for (i in 1:nrow(PA_G10))
{Nnse$value[i] <- PA_G10 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PA_G10$tree.nr[i]) %>% 
  .$N
PA_G10$N2[i] <- PA_G10$N[i] - Nnse$value[i]
print(PA_G10$N2[i])
}

Fig_PAG10_N2 <- PA_G10 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_PAG10_N <- PA_G10 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_PAG10_N2, Fig_PAG10_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

PA_G10$S2 <- 0

for (i in 1:nrow(PA_G10))
{Nnse$value[i] <- PA_G10 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PA_G10$tree.nr[i]) %>% 
  .$S
PA_G10$S2[i] <- PA_G10$S[i] - Nnse$value[i]
print(PA_G10$S2[i])
}

Fig_PAG10_S2 <- PA_G10 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_PAG10_S <- PA_G10 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_PAG10_S2, Fig_PAG10_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(PA_G10,file="Data/Genetic Distances/PA distances/PA_G10_distances_2.txt", row.names = F)



########### PB1!!!!!! #################


## Data import ##
PB1_G10 <- read.table("Data/Genetic Distances/PB1 distances/PB1_G10_distances.txt", header = T)
## Build redudant column ##
PB1_G10 <- PB1_G10[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(PB1_G10$tree.nr)
Nnse$value <- 0
PB1_G10$N2 <- 0

for (i in 1:nrow(PB1_G10))
{Nnse$value[i] <- PB1_G10 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB1_G10$tree.nr[i]) %>% 
  .$N
PB1_G10$N2[i] <- PB1_G10$N[i] - Nnse$value[i]
print(PB1_G10$N2[i])
}

Fig_PB1G10_N2 <- PB1_G10 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_PB1G10_N <- PB1_G10 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_PB1G10_N2, Fig_PB1G10_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

PB1_G10$S2 <- 0

for (i in 1:nrow(PB1_G10))
{Nnse$value[i] <- PB1_G10 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB1_G10$tree.nr[i]) %>% 
  .$S
PB1_G10$S2[i] <- PB1_G10$S[i] - Nnse$value[i]
print(PB1_G10$S2[i])
}

Fig_PB1G10_S2 <- PB1_G10 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_PB1G10_S <- PB1_G10 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_PB1G10_S2, Fig_PB1G10_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(PB1_G10,file="Data/Genetic Distances/PB1 distances/PB1_G10_distances_2.txt", row.names = F)


########### PB2!!!!!! #################


## Data import ##
PB2_G10 <- read.table("Data/Genetic Distances/PB2 distances/PB2_G10_distances.txt", header = T)
## Build redudant column ##
PB2_G10 <- PB2_G10[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(PB2_G10$tree.nr)
Nnse$value <- 0
PB2_G10$N2 <- 0

for (i in 1:nrow(PB2_G10))
{Nnse$value[i] <- PB2_G10 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB2_G10$tree.nr[i]) %>% 
  .$N
PB2_G10$N2[i] <- PB2_G10$N[i] - Nnse$value[i]
print(PB2_G10$N2[i])
}

Fig_PB2G10_N2 <- PB2_G10 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_PB2G10_N <- PB2_G10 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_PB2G10_N2, Fig_PB2G10_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

PB2_G10$S2 <- 0

for (i in 1:nrow(PB2_G10))
{Nnse$value[i] <- PB2_G10 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB2_G10$tree.nr[i]) %>% 
  .$S
PB2_G10$S2[i] <- PB2_G10$S[i] - Nnse$value[i]
print(PB2_G10$S2[i])
}

Fig_PB2G10_S2 <- PB2_G10 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_PB2G10_S <- PB2_G10 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_PB2G10_S2, Fig_PB2G10_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(PB2_G10,file="Data/Genetic Distances/PB2 distances/PB2_G10_distances_2.txt", row.names = F)


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
  summarise(Nmean = mean(N), Nmedian = median(N), Nerror = mean(N)-median(N))

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
  summarise(Smean = mean(S), Smedian = median(S), Serror = mean(S)-median(S))

## Export the intermediate results
G10_NS <- cbind(N_G10_Ndist, S_G10_Ndist[,4:6])

## Make year column
G10_NS$Year <- as.numeric(str_sub(gsub("_NA", "", G10_NS$compareStrain),-4,-1))

write.table(G10_NS,file="Data/Genetic Distances/GeneCompete/G10_NS.txt")


