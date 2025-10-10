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
HA_G4 <- read.table("Data/Genetic Distances/HA distances/HA_G4_distances.txt", header = T)
## Build redudant column ##
HA_G4 <- HA_G4[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(HA_G4$tree.nr)
Nnse$value <- 0
HA_G4$N2 <- 0

for (i in 1:nrow(HA_G4))
  {Nnse$value[i] <- HA_G4 %>% 
     filter(vaccineStrain == compareStrain) %>% 
     filter(tree.nr == HA_G4$tree.nr[i]) %>% 
    .$N
  HA_G4$N2[i] <- HA_G4$N[i] - Nnse$value[i]
  print(HA_G4$N2[i])
}

Fig_HAG4_N2 <- HA_G4 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_HAG4_N <- HA_G4 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_HAG4_N2, Fig_HAG4_N, 
          labels = c("N2", "N"),
          nrow = 2)


## S subs ##

HA_G4$S2 <- 0

for (i in 1:nrow(HA_G4))
{Nnse$value[i] <- HA_G4 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == HA_G4$tree.nr[i]) %>% 
  .$S
HA_G4$S2[i] <- HA_G4$S[i] - Nnse$value[i]
print(HA_G4$S2[i])
}

Fig_HAG4_S2 <- HA_G4 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_HAG4_S <- HA_G4 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_HAG4_S2, Fig_HAG4_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(HA_G4,file="Data/Genetic Distances/HA distances/HA_G4_distances_2.txt")


########### NA !!!!!! #################

## Data import ##
NA_G4 <- read.table("Data/Genetic Distances/NA distances/NA_G4_distances.txt", header = T)
## Build redudant column ##
NA_G4 <- NA_G4[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NA_G4$tree.nr)
Nnse$value <- 0
NA_G4$N2 <- 0

for (i in 1:nrow(NA_G4))
{Nnse$value[i] <- NA_G4 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NA_G4$tree.nr[i]) %>% 
  .$N
NA_G4$N2[i] <- NA_G4$N[i] - Nnse$value[i]
print(NA_G4$N2[i])
}

Fig_NAG4_N2 <- NA_G4 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_NAG4_N <- NA_G4 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_NAG4_N2, Fig_NAG4_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

NA_G4$S2 <- 0

for (i in 1:nrow(NA_G4))
{Nnse$value[i] <- NA_G4 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NA_G4$tree.nr[i]) %>% 
  .$S
NA_G4$S2[i] <- NA_G4$S[i] - Nnse$value[i]
print(NA_G4$S2[i])
}

Fig_NAG4_S2 <- NA_G4 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_NAG4_S <- NA_G4 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_NAG4_S2, Fig_NAG4_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(NA_G4,file="Data/Genetic Distances/NA distances/NA_G4_distances_2.txt")


########### M!!!!!! #################

## Data import ##
M_G4 <- read.table("Data/Genetic Distances/M distances/M_G4_distances.txt", header = T)
## Build redudant column ##
M_G4 <- M_G4[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(M_G4$tree.nr)
Nnse$value <- 0
M_G4$N2 <- 0

for (i in 1:nrow(M_G4))
{Nnse$value[i] <- M_G4 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == M_G4$tree.nr[i]) %>% 
  .$N
M_G4$N2[i] <- M_G4$N[i] - Nnse$value[i]
print(M_G4$N2[i])
}

Fig_MG4_N2 <- M_G4 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_MG4_N <- M_G4 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_MG4_N2, Fig_MG4_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

M_G4$S2 <- 0

for (i in 1:nrow(M_G4))
{Nnse$value[i] <- M_G4 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == M_G4$tree.nr[i]) %>% 
  .$S
M_G4$S2[i] <- M_G4$S[i] - Nnse$value[i]
print(M_G4$S2[i])
}

Fig_MG4_S2 <- M_G4 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_MG4_S <- M_G4 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_MG4_S2, Fig_MG4_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(M_G4,file="Data/Genetic Distances/M distances/M_G4_distances_2.txt")


########### NP!!!!!! #################

## Data import ##
NP_G4 <- read.table("Data/Genetic Distances/NP distances/NP_G4_distances.txt", header = T)
## Build redudant column ##
NP_G4 <- NP_G4[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NP_G4$tree.nr)
Nnse$value <- 0
NP_G4$N2 <- 0

for (i in 1:nrow(NP_G4))
{Nnse$value[i] <- NP_G4 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NP_G4$tree.nr[i]) %>% 
  .$N
NP_G4$N2[i] <- NP_G4$N[i] - Nnse$value[i]
print(NP_G4$N2[i])
}

Fig_NPG4_N2 <- NP_G4 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_NPG4_N <- NP_G4 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_NPG4_N2, Fig_NPG4_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

NP_G4$S2 <- 0

for (i in 1:nrow(NP_G4))
{Nnse$value[i] <- NP_G4 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NP_G4$tree.nr[i]) %>% 
  .$S
NP_G4$S2[i] <- NP_G4$S[i] - Nnse$value[i]
print(NP_G4$S2[i])
}

Fig_NPG4_S2 <- NP_G4 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_NPG4_S <- NP_G4 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_NPG4_S2, Fig_NPG4_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(NP_G4,file="Data/Genetic Distances/NP distances/NP_G4_distances_2.txt", row.names = F)


########### NS!!!!!! #################


## Data import ##
NS_G4 <- read.table("Data/Genetic Distances/NS distances/NS_G4_distances.txt", header = T)
## Build redudant column ##
NS_G4 <- NS_G4[,1:5]


## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NS_G4$tree.nr)
Nnse$value <- 0
NS_G4$N2 <- NS_G4$N

for (i in 1:nrow(NS_G4))
  {Nnse$value[i] <- NS_G4 %>% 
     filter(vaccineStrain == compareStrain) %>% 
     filter(tree.nr == NS_G4$tree.nr[i]) %>% 
    .$N
  NS_G4$N2[i] <- NS_G4$N[i] - Nnse$value[i]
  print(NS_G4$N2[i])
}

Fig_NSG4_N2 <- NS_G4 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_NSG4_N <- NS_G4 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_NSG4_N2, Fig_NSG4_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

NS_G4$S2 <- NS_G4$S

for (i in 1:nrow(NS_G4))
{Nnse$value[i] <- NS_G4 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NS_G4$tree.nr[i]) %>% 
  .$S
NS_G4$S2[i] <- NS_G4$S[i] - Nnse$value[i]
print(NS_G4$S2[i])
}

Fig_NSG4_S2 <- NS_G4 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_NSG4_S <- NS_G4 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_NSG4_S2, Fig_NSG4_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(NS_G4,file="Data/Genetic Distances/NS distances/NS_G4_distances_2.txt")


########### PA!!!!!! #################

## Data import ##
PA_G4 <- read.table("Data/Genetic Distances/PA distances/PA_G4_distances.txt", header = T)
## Build redudant column ##
PA_G4 <- PA_G4[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(PA_G4$tree.nr)
Nnse$value <- 0
PA_G4$N2 <- 0

for (i in 1:nrow(PA_G4))
{Nnse$value[i] <- PA_G4 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PA_G4$tree.nr[i]) %>% 
  .$N
PA_G4$N2[i] <- PA_G4$N[i] - Nnse$value[i]
print(PA_G4$N2[i])
}

Fig_PAG4_N2 <- PA_G4 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_PAG4_N <- PA_G4 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_PAG4_N2, Fig_PAG4_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

PA_G4$S2 <- 0

for (i in 1:nrow(PA_G4))
{Nnse$value[i] <- PA_G4 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PA_G4$tree.nr[i]) %>% 
  .$S
PA_G4$S2[i] <- PA_G4$S[i] - Nnse$value[i]
print(PA_G4$S2[i])
}

Fig_PAG4_S2 <- PA_G4 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_PAG4_S <- PA_G4 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_PAG4_S2, Fig_PAG4_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(PA_G4,file="Data/Genetic Distances/PA distances/PA_G4_distances_2.txt", row.names = F)



########### PB1!!!!!! #################


## Data import ##
PB1_G4 <- read.table("Data/Genetic Distances/PB1 distances/PB1_G4_distances.txt", header = T)
## Build redudant column ##
PB1_G4 <- PB1_G4[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(PB1_G4$tree.nr)
Nnse$value <- 0
PB1_G4$N2 <- 0

for (i in 1:nrow(PB1_G4))
{Nnse$value[i] <- PB1_G4 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB1_G4$tree.nr[i]) %>% 
  .$N
PB1_G4$N2[i] <- PB1_G4$N[i] - Nnse$value[i]
print(PB1_G4$N2[i])
}

Fig_PB1G4_N2 <- PB1_G4 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_PB1G4_N <- PB1_G4 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_PB1G4_N2, Fig_PB1G4_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

PB1_G4$S2 <- 0

for (i in 1:nrow(PB1_G4))
{Nnse$value[i] <- PB1_G4 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB1_G4$tree.nr[i]) %>% 
  .$S
PB1_G4$S2[i] <- PB1_G4$S[i] - Nnse$value[i]
print(PB1_G4$S2[i])
}

Fig_PB1G4_S2 <- PB1_G4 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_PB1G4_S <- PB1_G4 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_PB1G4_S2, Fig_PB1G4_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(PB1_G4,file="Data/Genetic Distances/PB1 distances/PB1_G4_distances_2.txt", row.names = F)


########### PB2!!!!!! #################


## Data import ##
PB2_G4 <- read.table("Data/Genetic Distances/PB2 distances/PB2_G4_distances.txt", header = T)
## Build redudant column ##
PB2_G4 <- PB2_G4[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(PB2_G4$tree.nr)
Nnse$value <- 0
PB2_G4$N2 <- 0

for (i in 1:nrow(PB2_G4))
{Nnse$value[i] <- PB2_G4 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB2_G4$tree.nr[i]) %>% 
  .$N
PB2_G4$N2[i] <- PB2_G4$N[i] - Nnse$value[i]
print(PB2_G4$N2[i])
}

Fig_PB2G4_N2 <- PB2_G4 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_PB2G4_N <- PB2_G4 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_PB2G4_N2, Fig_PB2G4_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

PB2_G4$S2 <- 0

for (i in 1:nrow(PB2_G4))
{Nnse$value[i] <- PB2_G4 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB2_G4$tree.nr[i]) %>% 
  .$S
PB2_G4$S2[i] <- PB2_G4$S[i] - Nnse$value[i]
print(PB2_G4$S2[i])
}

Fig_PB2G4_S2 <- PB2_G4 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_PB2G4_S <- PB2_G4 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_PB2G4_S2, Fig_PB2G4_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(PB2_G4,file="Data/Genetic Distances/PB2 distances/PB2_G4_distances_2.txt", row.names = F)

