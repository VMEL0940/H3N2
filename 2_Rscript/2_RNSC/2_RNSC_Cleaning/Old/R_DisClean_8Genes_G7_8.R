### Influ H3N2 GenDist Proj ###
## Package ##
library(tidyverse)
library(factoextra)
library(FactoMineR)
library(ggpubr)

## Set working directory ##
setwd("/01_GenDisFlu/")

########### HA !!!!!! #################

## Data import ##
HA_G7 <- read.table("Data/Genetic Distances/HA distances/HA_G7_distances.txt", header = T)
## Build redudant column ##
HA_G7 <- HA_G7[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(HA_G7$tree.nr)
Nnse$value <- 0
HA_G7$N2 <- 0

for (i in 1:nrow(HA_G7))
  {Nnse$value[i] <- HA_G7 %>% 
     filter(vaccineStrain == compareStrain) %>% 
     filter(tree.nr == HA_G7$tree.nr[i]) %>% 
    .$N
  HA_G7$N2[i] <- HA_G7$N[i] - Nnse$value[i]
  print(HA_G7$N2[i])
}

Fig_HAG7_N2 <- HA_G7 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_HAG7_N <- HA_G7 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_HAG7_N2, Fig_HAG7_N, 
          labels = c("N2", "N"),
          nrow = 2)


## S subs ##

HA_G7$S2 <- 0

for (i in 1:nrow(HA_G7))
{Nnse$value[i] <- HA_G7 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == HA_G7$tree.nr[i]) %>% 
  .$S
HA_G7$S2[i] <- HA_G7$S[i] - Nnse$value[i]
print(HA_G7$S2[i])
}

Fig_HAG7_S2 <- HA_G7 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_HAG7_S <- HA_G7 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_HAG7_S2, Fig_HAG7_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(HA_G7,file="Data/Genetic Distances/HA distances/HA_G7_distances_2.txt")


########### NA !!!!!! #################

## Data import ##
NA_G7 <- read.table("Data/Genetic Distances/NA distances/NA_G7_distances.txt", header = T)
## Build redudant column ##
NA_G7 <- NA_G7[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NA_G7$tree.nr)
Nnse$value <- 0
NA_G7$N2 <- 0

for (i in 1:nrow(NA_G7))
{Nnse$value[i] <- NA_G7 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NA_G7$tree.nr[i]) %>% 
  .$N
NA_G7$N2[i] <- NA_G7$N[i] - Nnse$value[i]
print(NA_G7$N2[i])
}

Fig_NAG7_N2 <- NA_G7 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_NAG7_N <- NA_G7 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_NAG7_N2, Fig_NAG7_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

NA_G7$S2 <- 0

for (i in 1:nrow(NA_G7))
{Nnse$value[i] <- NA_G7 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NA_G7$tree.nr[i]) %>% 
  .$S
NA_G7$S2[i] <- NA_G7$S[i] - Nnse$value[i]
print(NA_G7$S2[i])
}

Fig_NAG7_S2 <- NA_G7 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_NAG7_S <- NA_G7 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_NAG7_S2, Fig_NAG7_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(NA_G7,file="Data/Genetic Distances/NA distances/NA_G7_distances_2.txt")


########### M!!!!!! #################

## Data import ##
M_G7 <- read.table("Data/Genetic Distances/M distances/M_G7_distances.txt", header = T)
## Build redudant column ##
M_G7 <- M_G7[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(M_G7$tree.nr)
Nnse$value <- 0
M_G7$N2 <- 0

for (i in 1:nrow(M_G7))
{Nnse$value[i] <- M_G7 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == M_G7$tree.nr[i]) %>% 
  .$N
M_G7$N2[i] <- M_G7$N[i] - Nnse$value[i]
print(M_G7$N2[i])
}

Fig_MG7_N2 <- M_G7 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_MG7_N <- M_G7 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_MG7_N2, Fig_MG7_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

M_G7$S2 <- 0

for (i in 1:nrow(M_G7))
{Nnse$value[i] <- M_G7 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == M_G7$tree.nr[i]) %>% 
  .$S
M_G7$S2[i] <- M_G7$S[i] - Nnse$value[i]
print(M_G7$S2[i])
}

Fig_MG7_S2 <- M_G7 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_MG7_S <- M_G7 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_MG7_S2, Fig_MG7_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(M_G7,file="Data/Genetic Distances/M distances/M_G7_distances_2.txt")


########### NP!!!!!! #################

## Data import ##
NP_G7 <- read.table("Data/Genetic Distances/NP distances/NP_G7_distances.txt", header = T)
## Build redudant column ##
NP_G7 <- NP_G7[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NP_G7$tree.nr)
Nnse$value <- 0
NP_G7$N2 <- 0

for (i in 1:nrow(NP_G7))
{Nnse$value[i] <- NP_G7 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NP_G7$tree.nr[i]) %>% 
  .$N
NP_G7$N2[i] <- NP_G7$N[i] - Nnse$value[i]
print(NP_G7$N2[i])
}

Fig_NPG7_N2 <- NP_G7 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_NPG7_N <- NP_G7 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_NPG7_N2, Fig_NPG7_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

NP_G7$S2 <- 0

for (i in 1:nrow(NP_G7))
{Nnse$value[i] <- NP_G7 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NP_G7$tree.nr[i]) %>% 
  .$S
NP_G7$S2[i] <- NP_G7$S[i] - Nnse$value[i]
print(NP_G7$S2[i])
}

Fig_NPG7_S2 <- NP_G7 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_NPG7_S <- NP_G7 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_NPG7_S2, Fig_NPG7_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(NP_G7,file="Data/Genetic Distances/NP distances/NP_G7_distances_2.txt", row.names = F)


########### NS!!!!!! #################


## Data import ##
NS_G7 <- read.table("Data/Genetic Distances/NS distances/NS_G7_distances.txt", header = T)
## Build redudant column ##
NS_G7 <- NS_G7[,1:5]


## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NS_G7$tree.nr)
Nnse$value <- 0
NS_G7$N2 <- NS_G7$N

for (i in 1:nrow(NS_G7))
  {Nnse$value[i] <- NS_G7 %>% 
     filter(vaccineStrain == compareStrain) %>% 
     filter(tree.nr == NS_G7$tree.nr[i]) %>% 
    .$N
  NS_G7$N2[i] <- NS_G7$N[i] - Nnse$value[i]
  print(NS_G7$N2[i])
}

Fig_NSG7_N2 <- NS_G7 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_NSG7_N <- NS_G7 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_NSG7_N2, Fig_NSG7_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

NS_G7$S2 <- NS_G7$S

for (i in 1:nrow(NS_G7))
{Nnse$value[i] <- NS_G7 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NS_G7$tree.nr[i]) %>% 
  .$S
NS_G7$S2[i] <- NS_G7$S[i] - Nnse$value[i]
print(NS_G7$S2[i])
}

Fig_NSG7_S2 <- NS_G7 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_NSG7_S <- NS_G7 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_NSG7_S2, Fig_NSG7_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(NS_G7,file="Data/Genetic Distances/NS distances/NS_G7_distances_2.txt")


########### PA!!!!!! #################

## Data import ##
PA_G7 <- read.table("Data/Genetic Distances/PA distances/PA_G7_distances.txt", header = T)
## Build redudant column ##
PA_G7 <- PA_G7[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(PA_G7$tree.nr)
Nnse$value <- 0
PA_G7$N2 <- 0

for (i in 1:nrow(PA_G7))
{Nnse$value[i] <- PA_G7 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PA_G7$tree.nr[i]) %>% 
  .$N
PA_G7$N2[i] <- PA_G7$N[i] - Nnse$value[i]
print(PA_G7$N2[i])
}

Fig_PAG7_N2 <- PA_G7 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_PAG7_N <- PA_G7 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_PAG7_N2, Fig_PAG7_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

PA_G7$S2 <- 0

for (i in 1:nrow(PA_G7))
{Nnse$value[i] <- PA_G7 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PA_G7$tree.nr[i]) %>% 
  .$S
PA_G7$S2[i] <- PA_G7$S[i] - Nnse$value[i]
print(PA_G7$S2[i])
}

Fig_PAG7_S2 <- PA_G7 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_PAG7_S <- PA_G7 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_PAG7_S2, Fig_PAG7_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(PA_G7,file="Data/Genetic Distances/PA distances/PA_G7_distances_2.txt", row.names = F)



########### PB1!!!!!! #################


## Data import ##
PB1_G7 <- read.table("Data/Genetic Distances/PB1 distances/PB1_G7_distances.txt", header = T)
## Build redudant column ##
PB1_G7 <- PB1_G7[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(PB1_G7$tree.nr)
Nnse$value <- 0
PB1_G7$N2 <- 0

for (i in 1:nrow(PB1_G7))
{Nnse$value[i] <- PB1_G7 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB1_G7$tree.nr[i]) %>% 
  .$N
PB1_G7$N2[i] <- PB1_G7$N[i] - Nnse$value[i]
print(PB1_G7$N2[i])
}

Fig_PB1G7_N2 <- PB1_G7 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_PB1G7_N <- PB1_G7 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_PB1G7_N2, Fig_PB1G7_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

PB1_G7$S2 <- 0

for (i in 1:nrow(PB1_G7))
{Nnse$value[i] <- PB1_G7 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB1_G7$tree.nr[i]) %>% 
  .$S
PB1_G7$S2[i] <- PB1_G7$S[i] - Nnse$value[i]
print(PB1_G7$S2[i])
}

Fig_PB1G7_S2 <- PB1_G7 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_PB1G7_S <- PB1_G7 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_PB1G7_S2, Fig_PB1G7_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(PB1_G7,file="Data/Genetic Distances/PB1 distances/PB1_G7_distances_2.txt", row.names = F)


########### PB2!!!!!! #################


## Data import ##
PB2_G7 <- read.table("Data/Genetic Distances/PB2 distances/PB2_G7_distances.txt", header = T)
## Build redudant column ##
PB2_G7 <- PB2_G7[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(PB2_G7$tree.nr)
Nnse$value <- 0
PB2_G7$N2 <- 0

for (i in 1:nrow(PB2_G7))
{Nnse$value[i] <- PB2_G7 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB2_G7$tree.nr[i]) %>% 
  .$N
PB2_G7$N2[i] <- PB2_G7$N[i] - Nnse$value[i]
print(PB2_G7$N2[i])
}

Fig_PB2G7_N2 <- PB2_G7 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_PB2G7_N <- PB2_G7 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_PB2G7_N2, Fig_PB2G7_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

PB2_G7$S2 <- 0

for (i in 1:nrow(PB2_G7))
{Nnse$value[i] <- PB2_G7 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB2_G7$tree.nr[i]) %>% 
  .$S
PB2_G7$S2[i] <- PB2_G7$S[i] - Nnse$value[i]
print(PB2_G7$S2[i])
}

Fig_PB2G7_S2 <- PB2_G7 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_PB2G7_S <- PB2_G7 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_PB2G7_S2, Fig_PB2G7_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(PB2_G7,file="Data/Genetic Distances/PB2 distances/PB2_G7_distances_2.txt", row.names = F)


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
  summarise(Nmean = mean(N), Nmedian = median(N), Nerror = mean(N)-median(N))

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
  summarise(Smean = mean(S), Smedian = median(S), Serror = mean(S)-median(S))

## Export the intermediate results
G7_NS <- cbind(N_G7_Ndist, S_G7_Ndist[,4:6])

## Make year column
G7_NS$Year <- as.numeric(str_sub(gsub("_NA", "", G7_NS$compareStrain),-4,-1))

write.table(G7_NS,file="Data/Genetic Distances/GeneCompete/G7_NS.txt")



## Data import ##
HA_G8 <- read.table("Data/Genetic Distances/HA distances/HA_G8_distances.txt", header = T)
## Build redudant column ##
HA_G8 <- HA_G8[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(HA_G8$tree.nr)
Nnse$value <- 0
HA_G8$N2 <- 0

for (i in 1:nrow(HA_G8))
{Nnse$value[i] <- HA_G8 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == HA_G8$tree.nr[i]) %>% 
  .$N
HA_G8$N2[i] <- HA_G8$N[i] - Nnse$value[i]
print(HA_G8$N2[i])
}

Fig_HAG8_N2 <- HA_G8 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_HAG8_N <- HA_G8 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_HAG8_N2, Fig_HAG8_N, 
          labels = c("N2", "N"),
          nrow = 2)


## S subs ##

HA_G8$S2 <- 0

for (i in 1:nrow(HA_G8))
{Nnse$value[i] <- HA_G8 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == HA_G8$tree.nr[i]) %>% 
  .$S
HA_G8$S2[i] <- HA_G8$S[i] - Nnse$value[i]
print(HA_G8$S2[i])
}

Fig_HAG8_S2 <- HA_G8 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_HAG8_S <- HA_G8 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_HAG8_S2, Fig_HAG8_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(HA_G8,file="Data/Genetic Distances/HA distances/HA_G8_distances_2.txt")


########### NA !!!!!! #################

## Data import ##
NA_G8 <- read.table("Data/Genetic Distances/NA distances/NA_G8_distances.txt", header = T)
## Build redudant column ##
NA_G8 <- NA_G8[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NA_G8$tree.nr)
Nnse$value <- 0
NA_G8$N2 <- 0

for (i in 1:nrow(NA_G8))
{Nnse$value[i] <- NA_G8 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NA_G8$tree.nr[i]) %>% 
  .$N
NA_G8$N2[i] <- NA_G8$N[i] - Nnse$value[i]
print(NA_G8$N2[i])
}

Fig_NAG8_N2 <- NA_G8 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_NAG8_N <- NA_G8 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_NAG8_N2, Fig_NAG8_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

NA_G8$S2 <- 0

for (i in 1:nrow(NA_G8))
{Nnse$value[i] <- NA_G8 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NA_G8$tree.nr[i]) %>% 
  .$S
NA_G8$S2[i] <- NA_G8$S[i] - Nnse$value[i]
print(NA_G8$S2[i])
}

Fig_NAG8_S2 <- NA_G8 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_NAG8_S <- NA_G8 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_NAG8_S2, Fig_NAG8_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(NA_G8,file="Data/Genetic Distances/NA distances/NA_G8_distances_2.txt")


########### M!!!!!! #################

## Data import ##
M_G8 <- read.table("Data/Genetic Distances/M distances/M_G8_distances.txt", header = T)
## Build redudant column ##
M_G8 <- M_G8[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(M_G8$tree.nr)
Nnse$value <- 0
M_G8$N2 <- 0

for (i in 1:nrow(M_G8))
{Nnse$value[i] <- M_G8 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == M_G8$tree.nr[i]) %>% 
  .$N
M_G8$N2[i] <- M_G8$N[i] - Nnse$value[i]
print(M_G8$N2[i])
}

Fig_MG8_N2 <- M_G8 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_MG8_N <- M_G8 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_MG8_N2, Fig_MG8_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

M_G8$S2 <- 0

for (i in 1:nrow(M_G8))
{Nnse$value[i] <- M_G8 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == M_G8$tree.nr[i]) %>% 
  .$S
M_G8$S2[i] <- M_G8$S[i] - Nnse$value[i]
print(M_G8$S2[i])
}

Fig_MG8_S2 <- M_G8 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_MG8_S <- M_G8 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_MG8_S2, Fig_MG8_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(M_G8,file="Data/Genetic Distances/M distances/M_G8_distances_2.txt")


########### NP!!!!!! #################

## Data import ##
NP_G8 <- read.table("Data/Genetic Distances/NP distances/NP_G8_distances.txt", header = T)
## Build redudant column ##
NP_G8 <- NP_G8[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NP_G8$tree.nr)
Nnse$value <- 0
NP_G8$N2 <- 0

for (i in 1:nrow(NP_G8))
{Nnse$value[i] <- NP_G8 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NP_G8$tree.nr[i]) %>% 
  .$N
NP_G8$N2[i] <- NP_G8$N[i] - Nnse$value[i]
print(NP_G8$N2[i])
}

Fig_NPG8_N2 <- NP_G8 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_NPG8_N <- NP_G8 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_NPG8_N2, Fig_NPG8_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

NP_G8$S2 <- 0

for (i in 1:nrow(NP_G8))
{Nnse$value[i] <- NP_G8 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NP_G8$tree.nr[i]) %>% 
  .$S
NP_G8$S2[i] <- NP_G8$S[i] - Nnse$value[i]
print(NP_G8$S2[i])
}

Fig_NPG8_S2 <- NP_G8 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_NPG8_S <- NP_G8 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_NPG8_S2, Fig_NPG8_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(NP_G8,file="Data/Genetic Distances/NP distances/NP_G8_distances_2.txt", row.names = F)


########### NS!!!!!! #################


## Data import ##
NS_G8 <- read.table("Data/Genetic Distances/NS distances/NS_G8_distances.txt", header = T)
## Build redudant column ##
NS_G8 <- NS_G8[,1:5]


## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NS_G8$tree.nr)
Nnse$value <- 0
NS_G8$N2 <- NS_G8$N

for (i in 1:nrow(NS_G8))
{Nnse$value[i] <- NS_G8 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NS_G8$tree.nr[i]) %>% 
  .$N
NS_G8$N2[i] <- NS_G8$N[i] - Nnse$value[i]
print(NS_G8$N2[i])
}

Fig_NSG8_N2 <- NS_G8 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_NSG8_N <- NS_G8 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_NSG8_N2, Fig_NSG8_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

NS_G8$S2 <- NS_G8$S

for (i in 1:nrow(NS_G8))
{Nnse$value[i] <- NS_G8 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NS_G8$tree.nr[i]) %>% 
  .$S
NS_G8$S2[i] <- NS_G8$S[i] - Nnse$value[i]
print(NS_G8$S2[i])
}

Fig_NSG8_S2 <- NS_G8 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_NSG8_S <- NS_G8 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_NSG8_S2, Fig_NSG8_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(NS_G8,file="Data/Genetic Distances/NS distances/NS_G8_distances_2.txt")


########### PA!!!!!! #################

## Data import ##
PA_G8 <- read.table("Data/Genetic Distances/PA distances/PA_G8_distances.txt", header = T)
## Build redudant column ##
PA_G8 <- PA_G8[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(PA_G8$tree.nr)
Nnse$value <- 0
PA_G8$N2 <- 0

for (i in 1:nrow(PA_G8))
{Nnse$value[i] <- PA_G8 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PA_G8$tree.nr[i]) %>% 
  .$N
PA_G8$N2[i] <- PA_G8$N[i] - Nnse$value[i]
print(PA_G8$N2[i])
}

Fig_PAG8_N2 <- PA_G8 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_PAG8_N <- PA_G8 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_PAG8_N2, Fig_PAG8_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

PA_G8$S2 <- 0

for (i in 1:nrow(PA_G8))
{Nnse$value[i] <- PA_G8 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PA_G8$tree.nr[i]) %>% 
  .$S
PA_G8$S2[i] <- PA_G8$S[i] - Nnse$value[i]
print(PA_G8$S2[i])
}

Fig_PAG8_S2 <- PA_G8 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_PAG8_S <- PA_G8 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_PAG8_S2, Fig_PAG8_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(PA_G8,file="Data/Genetic Distances/PA distances/PA_G8_distances_2.txt", row.names = F)



########### PB1!!!!!! #################


## Data import ##
PB1_G8 <- read.table("Data/Genetic Distances/PB1 distances/PB1_G8_distances.txt", header = T)
## Build redudant column ##
PB1_G8 <- PB1_G8[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(PB1_G8$tree.nr)
Nnse$value <- 0
PB1_G8$N2 <- 0

for (i in 1:nrow(PB1_G8))
{Nnse$value[i] <- PB1_G8 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB1_G8$tree.nr[i]) %>% 
  .$N
PB1_G8$N2[i] <- PB1_G8$N[i] - Nnse$value[i]
print(PB1_G8$N2[i])
}

Fig_PB1G8_N2 <- PB1_G8 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_PB1G8_N <- PB1_G8 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_PB1G8_N2, Fig_PB1G8_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

PB1_G8$S2 <- 0

for (i in 1:nrow(PB1_G8))
{Nnse$value[i] <- PB1_G8 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB1_G8$tree.nr[i]) %>% 
  .$S
PB1_G8$S2[i] <- PB1_G8$S[i] - Nnse$value[i]
print(PB1_G8$S2[i])
}

Fig_PB1G8_S2 <- PB1_G8 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_PB1G8_S <- PB1_G8 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_PB1G8_S2, Fig_PB1G8_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(PB1_G8,file="Data/Genetic Distances/PB1 distances/PB1_G8_distances_2.txt", row.names = F)


########### PB2!!!!!! #################


## Data import ##
PB2_G8 <- read.table("Data/Genetic Distances/PB2 distances/PB2_G8_distances.txt", header = T)
## Build redudant column ##
PB2_G8 <- PB2_G8[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(PB2_G8$tree.nr)
Nnse$value <- 0
PB2_G8$N2 <- 0

for (i in 1:nrow(PB2_G8))
{Nnse$value[i] <- PB2_G8 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB2_G8$tree.nr[i]) %>% 
  .$N
PB2_G8$N2[i] <- PB2_G8$N[i] - Nnse$value[i]
print(PB2_G8$N2[i])
}

Fig_PB2G8_N2 <- PB2_G8 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_PB2G8_N <- PB2_G8 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_PB2G8_N2, Fig_PB2G8_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

PB2_G8$S2 <- 0

for (i in 1:nrow(PB2_G8))
{Nnse$value[i] <- PB2_G8 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB2_G8$tree.nr[i]) %>% 
  .$S
PB2_G8$S2[i] <- PB2_G8$S[i] - Nnse$value[i]
print(PB2_G8$S2[i])
}

Fig_PB2G8_S2 <- PB2_G8 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_PB2G8_S <- PB2_G8 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_PB2G8_S2, Fig_PB2G8_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(PB2_G8,file="Data/Genetic Distances/PB2 distances/PB2_G8_distances_2.txt", row.names = F)


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
  summarise(Nmean = mean(N), Nmedian = median(N), Nerror = mean(N)-median(N))

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
  summarise(Smean = mean(S), Smedian = median(S), Serror = mean(S)-median(S))

## Export the intermediate results
G8_NS <- cbind(N_G8_Ndist, S_G8_Ndist[,4:6])

## Make year column
G8_NS$Year <- as.numeric(str_sub(gsub("_NA", "", G8_NS$compareStrain),-4,-1))

write.table(G8_NS,file="Data/Genetic Distances/GeneCompete/G8_NS.txt")

