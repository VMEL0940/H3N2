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
HA_G5 <- read.table("Data/Genetic Distances/HA distances/HA_G5_distances.txt", header = T)
## Build redudant column ##
HA_G5 <- HA_G5[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(HA_G5$tree.nr)
Nnse$value <- 0
HA_G5$N2 <- 0

for (i in 1:nrow(HA_G5))
  {Nnse$value[i] <- HA_G5 %>% 
     filter(vaccineStrain == compareStrain) %>% 
     filter(tree.nr == HA_G5$tree.nr[i]) %>% 
    .$N
  HA_G5$N2[i] <- HA_G5$N[i] - Nnse$value[i]
  print(HA_G5$N2[i])
}

Fig_HAG5_N2 <- HA_G5 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_HAG5_N <- HA_G5 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_HAG5_N2, Fig_HAG5_N, 
          labels = c("N2", "N"),
          nrow = 2)


## S subs ##

HA_G5$S2 <- 0

for (i in 1:nrow(HA_G5))
{Nnse$value[i] <- HA_G5 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == HA_G5$tree.nr[i]) %>% 
  .$S
HA_G5$S2[i] <- HA_G5$S[i] - Nnse$value[i]
print(HA_G5$S2[i])
}

Fig_HAG5_S2 <- HA_G5 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_HAG5_S <- HA_G5 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_HAG5_S2, Fig_HAG5_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(HA_G5,file="Data/Genetic Distances/HA distances/HA_G5_distances_2.txt")


########### NA !!!!!! #################

## Data import ##
NA_G5 <- read.table("Data/Genetic Distances/NA distances/NA_G5_distances.txt", header = T)
## Build redudant column ##
NA_G5 <- NA_G5[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NA_G5$tree.nr)
Nnse$value <- 0
NA_G5$N2 <- 0

for (i in 1:nrow(NA_G5))
{Nnse$value[i] <- NA_G5 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NA_G5$tree.nr[i]) %>% 
  .$N
NA_G5$N2[i] <- NA_G5$N[i] - Nnse$value[i]
print(NA_G5$N2[i])
}

Fig_NAG5_N2 <- NA_G5 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_NAG5_N <- NA_G5 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_NAG5_N2, Fig_NAG5_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

NA_G5$S2 <- 0

for (i in 1:nrow(NA_G5))
{Nnse$value[i] <- NA_G5 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NA_G5$tree.nr[i]) %>% 
  .$S
NA_G5$S2[i] <- NA_G5$S[i] - Nnse$value[i]
print(NA_G5$S2[i])
}

Fig_NAG5_S2 <- NA_G5 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_NAG5_S <- NA_G5 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_NAG5_S2, Fig_NAG5_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(NA_G5,file="Data/Genetic Distances/NA distances/NA_G5_distances_2.txt")


########### M!!!!!! #################

## Data import ##
M_G5 <- read.table("Data/Genetic Distances/M distances/M_G5_distances.txt", header = T)
## Build redudant column ##
M_G5 <- M_G5[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(M_G5$tree.nr)
Nnse$value <- 0
M_G5$N2 <- 0

for (i in 1:nrow(M_G5))
{Nnse$value[i] <- M_G5 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == M_G5$tree.nr[i]) %>% 
  .$N
M_G5$N2[i] <- M_G5$N[i] - Nnse$value[i]
print(M_G5$N2[i])
}

Fig_MG5_N2 <- M_G5 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_MG5_N <- M_G5 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_MG5_N2, Fig_MG5_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

M_G5$S2 <- 0

for (i in 1:nrow(M_G5))
{Nnse$value[i] <- M_G5 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == M_G5$tree.nr[i]) %>% 
  .$S
M_G5$S2[i] <- M_G5$S[i] - Nnse$value[i]
print(M_G5$S2[i])
}

Fig_MG5_S2 <- M_G5 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_MG5_S <- M_G5 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_MG5_S2, Fig_MG5_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(M_G5,file="Data/Genetic Distances/M distances/M_G5_distances_2.txt")


########### NP!!!!!! #################

## Data import ##
NP_G5 <- read.table("Data/Genetic Distances/NP distances/NP_G5_distances.txt", header = T)
## Build redudant column ##
NP_G5 <- NP_G5[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NP_G5$tree.nr)
Nnse$value <- 0
NP_G5$N2 <- 0

for (i in 1:nrow(NP_G5))
{Nnse$value[i] <- NP_G5 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NP_G5$tree.nr[i]) %>% 
  .$N
NP_G5$N2[i] <- NP_G5$N[i] - Nnse$value[i]
print(NP_G5$N2[i])
}

Fig_NPG5_N2 <- NP_G5 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_NPG5_N <- NP_G5 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_NPG5_N2, Fig_NPG5_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

NP_G5$S2 <- 0

for (i in 1:nrow(NP_G5))
{Nnse$value[i] <- NP_G5 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NP_G5$tree.nr[i]) %>% 
  .$S
NP_G5$S2[i] <- NP_G5$S[i] - Nnse$value[i]
print(NP_G5$S2[i])
}

Fig_NPG5_S2 <- NP_G5 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_NPG5_S <- NP_G5 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_NPG5_S2, Fig_NPG5_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(NP_G5,file="Data/Genetic Distances/NP distances/NP_G5_distances_2.txt", row.names = F)


########### NS!!!!!! #################


## Data import ##
NS_G5 <- read.table("Data/Genetic Distances/NS distances/NS_G5_distances.txt", header = T)
## Build redudant column ##
NS_G5 <- NS_G5[,1:5]


## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NS_G5$tree.nr)
Nnse$value <- 0
NS_G5$N2 <- NS_G5$N

for (i in 1:nrow(NS_G5))
  {Nnse$value[i] <- NS_G5 %>% 
     filter(vaccineStrain == compareStrain) %>% 
     filter(tree.nr == NS_G5$tree.nr[i]) %>% 
    .$N
  NS_G5$N2[i] <- NS_G5$N[i] - Nnse$value[i]
  print(NS_G5$N2[i])
}

Fig_NSG5_N2 <- NS_G5 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_NSG5_N <- NS_G5 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_NSG5_N2, Fig_NSG5_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

NS_G5$S2 <- NS_G5$S

for (i in 1:nrow(NS_G5))
{Nnse$value[i] <- NS_G5 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NS_G5$tree.nr[i]) %>% 
  .$S
NS_G5$S2[i] <- NS_G5$S[i] - Nnse$value[i]
print(NS_G5$S2[i])
}

Fig_NSG5_S2 <- NS_G5 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_NSG5_S <- NS_G5 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_NSG5_S2, Fig_NSG5_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(NS_G5,file="Data/Genetic Distances/NS distances/NS_G5_distances_2.txt")


########### PA!!!!!! #################

## Data import ##
PA_G5 <- read.table("Data/Genetic Distances/PA distances/PA_G5_distances.txt", header = T)
## Build redudant column ##
PA_G5 <- PA_G5[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(PA_G5$tree.nr)
Nnse$value <- 0
PA_G5$N2 <- 0

for (i in 1:nrow(PA_G5))
{Nnse$value[i] <- PA_G5 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PA_G5$tree.nr[i]) %>% 
  .$N
PA_G5$N2[i] <- PA_G5$N[i] - Nnse$value[i]
print(PA_G5$N2[i])
}

Fig_PAG5_N2 <- PA_G5 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_PAG5_N <- PA_G5 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_PAG5_N2, Fig_PAG5_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

PA_G5$S2 <- 0

for (i in 1:nrow(PA_G5))
{Nnse$value[i] <- PA_G5 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PA_G5$tree.nr[i]) %>% 
  .$S
PA_G5$S2[i] <- PA_G5$S[i] - Nnse$value[i]
print(PA_G5$S2[i])
}

Fig_PAG5_S2 <- PA_G5 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_PAG5_S <- PA_G5 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_PAG5_S2, Fig_PAG5_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(PA_G5,file="Data/Genetic Distances/PA distances/PA_G5_distances_2.txt", row.names = F)



########### PB1!!!!!! #################


## Data import ##
PB1_G5 <- read.table("Data/Genetic Distances/PB1 distances/PB1_G5_distances.txt", header = T)
## Build redudant column ##
PB1_G5 <- PB1_G5[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(PB1_G5$tree.nr)
Nnse$value <- 0
PB1_G5$N2 <- 0

for (i in 1:nrow(PB1_G5))
{Nnse$value[i] <- PB1_G5 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB1_G5$tree.nr[i]) %>% 
  .$N
PB1_G5$N2[i] <- PB1_G5$N[i] - Nnse$value[i]
print(PB1_G5$N2[i])
}

Fig_PB1G5_N2 <- PB1_G5 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_PB1G5_N <- PB1_G5 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_PB1G5_N2, Fig_PB1G5_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

PB1_G5$S2 <- 0

for (i in 1:nrow(PB1_G5))
{Nnse$value[i] <- PB1_G5 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB1_G5$tree.nr[i]) %>% 
  .$S
PB1_G5$S2[i] <- PB1_G5$S[i] - Nnse$value[i]
print(PB1_G5$S2[i])
}

Fig_PB1G5_S2 <- PB1_G5 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_PB1G5_S <- PB1_G5 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_PB1G5_S2, Fig_PB1G5_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(PB1_G5,file="Data/Genetic Distances/PB1 distances/PB1_G5_distances_2.txt", row.names = F)


########### PB2!!!!!! #################


## Data import ##
PB2_G5 <- read.table("Data/Genetic Distances/PB2 distances/PB2_G5_distances.txt", header = T)
## Build redudant column ##
PB2_G5 <- PB2_G5[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(PB2_G5$tree.nr)
Nnse$value <- 0
PB2_G5$N2 <- 0

for (i in 1:nrow(PB2_G5))
{Nnse$value[i] <- PB2_G5 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB2_G5$tree.nr[i]) %>% 
  .$N
PB2_G5$N2[i] <- PB2_G5$N[i] - Nnse$value[i]
print(PB2_G5$N2[i])
}

Fig_PB2G5_N2 <- PB2_G5 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_PB2G5_N <- PB2_G5 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_PB2G5_N2, Fig_PB2G5_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

PB2_G5$S2 <- 0

for (i in 1:nrow(PB2_G5))
{Nnse$value[i] <- PB2_G5 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB2_G5$tree.nr[i]) %>% 
  .$S
PB2_G5$S2[i] <- PB2_G5$S[i] - Nnse$value[i]
print(PB2_G5$S2[i])
}

Fig_PB2G5_S2 <- PB2_G5 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_PB2G5_S <- PB2_G5 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_PB2G5_S2, Fig_PB2G5_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(PB2_G5,file="Data/Genetic Distances/PB2 distances/PB2_G5_distances_2.txt", row.names = F)


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
  summarise(Nmean = mean(N), Nmedian = median(N), Nerror = mean(N)-median(N))

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
  summarise(Smean = mean(S), Smedian = median(S), Serror = mean(S)-median(S))

## Export the intermediate results
G5_NS <- cbind(N_G5_Ndist, S_G5_Ndist[,4:6])

## Make year column
G5_NS$Year <- as.numeric(str_sub(gsub("_NA", "", G5_NS$compareStrain),-4,-1))

write.table(G5_NS,file="Data/Genetic Distances/GeneCompete/G5_NS.txt")


## Data import ##
HA_G6 <- read.table("Data/Genetic Distances/HA distances/HA_G6_distances.txt", header = T)
## Build redudant column ##
HA_G6 <- HA_G6[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(HA_G6$tree.nr)
Nnse$value <- 0
HA_G6$N2 <- 0

for (i in 1:nrow(HA_G6))
{Nnse$value[i] <- HA_G6 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == HA_G6$tree.nr[i]) %>% 
  .$N
HA_G6$N2[i] <- HA_G6$N[i] - Nnse$value[i]
print(HA_G6$N2[i])
}

Fig_HAG6_N2 <- HA_G6 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_HAG6_N <- HA_G6 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_HAG6_N2, Fig_HAG6_N, 
          labels = c("N2", "N"),
          nrow = 2)


## S subs ##

HA_G6$S2 <- 0

for (i in 1:nrow(HA_G6))
{Nnse$value[i] <- HA_G6 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == HA_G6$tree.nr[i]) %>% 
  .$S
HA_G6$S2[i] <- HA_G6$S[i] - Nnse$value[i]
print(HA_G6$S2[i])
}

Fig_HAG6_S2 <- HA_G6 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_HAG6_S <- HA_G6 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_HAG6_S2, Fig_HAG6_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(HA_G6,file="Data/Genetic Distances/HA distances/HA_G6_distances_2.txt")


########### NA !!!!!! #################

## Data import ##
NA_G6 <- read.table("Data/Genetic Distances/NA distances/NA_G6_distances.txt", header = T)
## Build redudant column ##
NA_G6 <- NA_G6[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NA_G6$tree.nr)
Nnse$value <- 0
NA_G6$N2 <- 0

for (i in 1:nrow(NA_G6))
{Nnse$value[i] <- NA_G6 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NA_G6$tree.nr[i]) %>% 
  .$N
NA_G6$N2[i] <- NA_G6$N[i] - Nnse$value[i]
print(NA_G6$N2[i])
}

Fig_NAG6_N2 <- NA_G6 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_NAG6_N <- NA_G6 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_NAG6_N2, Fig_NAG6_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

NA_G6$S2 <- 0

for (i in 1:nrow(NA_G6))
{Nnse$value[i] <- NA_G6 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NA_G6$tree.nr[i]) %>% 
  .$S
NA_G6$S2[i] <- NA_G6$S[i] - Nnse$value[i]
print(NA_G6$S2[i])
}

Fig_NAG6_S2 <- NA_G6 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_NAG6_S <- NA_G6 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_NAG6_S2, Fig_NAG6_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(NA_G6,file="Data/Genetic Distances/NA distances/NA_G6_distances_2.txt")


########### M!!!!!! #################

## Data import ##
M_G6 <- read.table("Data/Genetic Distances/M distances/M_G6_distances.txt", header = T)
## Build redudant column ##
M_G6 <- M_G6[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(M_G6$tree.nr)
Nnse$value <- 0
M_G6$N2 <- 0

for (i in 1:nrow(M_G6))
{Nnse$value[i] <- M_G6 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == M_G6$tree.nr[i]) %>% 
  .$N
M_G6$N2[i] <- M_G6$N[i] - Nnse$value[i]
print(M_G6$N2[i])
}

Fig_MG6_N2 <- M_G6 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_MG6_N <- M_G6 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_MG6_N2, Fig_MG6_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

M_G6$S2 <- 0

for (i in 1:nrow(M_G6))
{Nnse$value[i] <- M_G6 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == M_G6$tree.nr[i]) %>% 
  .$S
M_G6$S2[i] <- M_G6$S[i] - Nnse$value[i]
print(M_G6$S2[i])
}

Fig_MG6_S2 <- M_G6 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_MG6_S <- M_G6 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_MG6_S2, Fig_MG6_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(M_G6,file="Data/Genetic Distances/M distances/M_G6_distances_2.txt")


########### NP!!!!!! #################

## Data import ##
NP_G6 <- read.table("Data/Genetic Distances/NP distances/NP_G6_distances.txt", header = T)
## Build redudant column ##
NP_G6 <- NP_G6[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NP_G6$tree.nr)
Nnse$value <- 0
NP_G6$N2 <- 0

for (i in 1:nrow(NP_G6))
{Nnse$value[i] <- NP_G6 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NP_G6$tree.nr[i]) %>% 
  .$N
NP_G6$N2[i] <- NP_G6$N[i] - Nnse$value[i]
print(NP_G6$N2[i])
}

Fig_NPG6_N2 <- NP_G6 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_NPG6_N <- NP_G6 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_NPG6_N2, Fig_NPG6_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

NP_G6$S2 <- 0

for (i in 1:nrow(NP_G6))
{Nnse$value[i] <- NP_G6 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NP_G6$tree.nr[i]) %>% 
  .$S
NP_G6$S2[i] <- NP_G6$S[i] - Nnse$value[i]
print(NP_G6$S2[i])
}

Fig_NPG6_S2 <- NP_G6 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_NPG6_S <- NP_G6 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_NPG6_S2, Fig_NPG6_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(NP_G6,file="Data/Genetic Distances/NP distances/NP_G6_distances_2.txt", row.names = F)


########### NS!!!!!! #################


## Data import ##
NS_G6 <- read.table("Data/Genetic Distances/NS distances/NS_G6_distances.txt", header = T)
## Build redudant column ##
NS_G6 <- NS_G6[,1:5]


## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NS_G6$tree.nr)
Nnse$value <- 0
NS_G6$N2 <- NS_G6$N

for (i in 1:nrow(NS_G6))
{Nnse$value[i] <- NS_G6 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NS_G6$tree.nr[i]) %>% 
  .$N
NS_G6$N2[i] <- NS_G6$N[i] - Nnse$value[i]
print(NS_G6$N2[i])
}

Fig_NSG6_N2 <- NS_G6 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_NSG6_N <- NS_G6 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_NSG6_N2, Fig_NSG6_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

NS_G6$S2 <- NS_G6$S

for (i in 1:nrow(NS_G6))
{Nnse$value[i] <- NS_G6 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NS_G6$tree.nr[i]) %>% 
  .$S
NS_G6$S2[i] <- NS_G6$S[i] - Nnse$value[i]
print(NS_G6$S2[i])
}

Fig_NSG6_S2 <- NS_G6 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_NSG6_S <- NS_G6 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_NSG6_S2, Fig_NSG6_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(NS_G6,file="Data/Genetic Distances/NS distances/NS_G6_distances_2.txt")


########### PA!!!!!! #################

## Data import ##
PA_G6 <- read.table("Data/Genetic Distances/PA distances/PA_G6_distances.txt", header = T)
## Build redudant column ##
PA_G6 <- PA_G6[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(PA_G6$tree.nr)
Nnse$value <- 0
PA_G6$N2 <- 0

for (i in 1:nrow(PA_G6))
{Nnse$value[i] <- PA_G6 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PA_G6$tree.nr[i]) %>% 
  .$N
PA_G6$N2[i] <- PA_G6$N[i] - Nnse$value[i]
print(PA_G6$N2[i])
}

Fig_PAG6_N2 <- PA_G6 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_PAG6_N <- PA_G6 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_PAG6_N2, Fig_PAG6_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

PA_G6$S2 <- 0

for (i in 1:nrow(PA_G6))
{Nnse$value[i] <- PA_G6 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PA_G6$tree.nr[i]) %>% 
  .$S
PA_G6$S2[i] <- PA_G6$S[i] - Nnse$value[i]
print(PA_G6$S2[i])
}

Fig_PAG6_S2 <- PA_G6 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_PAG6_S <- PA_G6 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_PAG6_S2, Fig_PAG6_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(PA_G6,file="Data/Genetic Distances/PA distances/PA_G6_distances_2.txt", row.names = F)



########### PB1!!!!!! #################


## Data import ##
PB1_G6 <- read.table("Data/Genetic Distances/PB1 distances/PB1_G6_distances.txt", header = T)
## Build redudant column ##
PB1_G6 <- PB1_G6[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(PB1_G6$tree.nr)
Nnse$value <- 0
PB1_G6$N2 <- 0

for (i in 1:nrow(PB1_G6))
{Nnse$value[i] <- PB1_G6 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB1_G6$tree.nr[i]) %>% 
  .$N
PB1_G6$N2[i] <- PB1_G6$N[i] - Nnse$value[i]
print(PB1_G6$N2[i])
}

Fig_PB1G6_N2 <- PB1_G6 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_PB1G6_N <- PB1_G6 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_PB1G6_N2, Fig_PB1G6_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

PB1_G6$S2 <- 0

for (i in 1:nrow(PB1_G6))
{Nnse$value[i] <- PB1_G6 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB1_G6$tree.nr[i]) %>% 
  .$S
PB1_G6$S2[i] <- PB1_G6$S[i] - Nnse$value[i]
print(PB1_G6$S2[i])
}

Fig_PB1G6_S2 <- PB1_G6 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_PB1G6_S <- PB1_G6 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_PB1G6_S2, Fig_PB1G6_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(PB1_G6,file="Data/Genetic Distances/PB1 distances/PB1_G6_distances_2.txt", row.names = F)


########### PB2!!!!!! #################


## Data import ##
PB2_G6 <- read.table("Data/Genetic Distances/PB2 distances/PB2_G6_distances.txt", header = T)
## Build redudant column ##
PB2_G6 <- PB2_G6[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(PB2_G6$tree.nr)
Nnse$value <- 0
PB2_G6$N2 <- 0

for (i in 1:nrow(PB2_G6))
{Nnse$value[i] <- PB2_G6 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB2_G6$tree.nr[i]) %>% 
  .$N
PB2_G6$N2[i] <- PB2_G6$N[i] - Nnse$value[i]
print(PB2_G6$N2[i])
}

Fig_PB2G6_N2 <- PB2_G6 %>%  
  ggplot(aes(x=as.factor(N2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_PB2G6_N <- PB2_G6 %>%  
  ggplot(aes(x=as.factor(N))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_PB2G6_N2, Fig_PB2G6_N, 
          labels = c("N2", "N"),
          nrow = 2)
## S subs ##

PB2_G6$S2 <- 0

for (i in 1:nrow(PB2_G6))
{Nnse$value[i] <- PB2_G6 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB2_G6$tree.nr[i]) %>% 
  .$S
PB2_G6$S2[i] <- PB2_G6$S[i] - Nnse$value[i]
print(PB2_G6$S2[i])
}

Fig_PB2G6_S2 <- PB2_G6 %>%  
  ggplot(aes(x=as.factor(S2))) +
  geom_histogram(stat = "count") +
  theme_bw()

Fig_PB2G6_S <- PB2_G6 %>%  
  ggplot(aes(x=as.factor(S))) +
  geom_histogram(stat = "count") +
  theme_bw()

ggarrange(Fig_PB2G6_S2, Fig_PB2G6_S, 
          labels = c("S2", "S"),
          nrow = 2)

write.table(PB2_G6,file="Data/Genetic Distances/PB2 distances/PB2_G6_distances_2.txt", row.names = F)


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
  summarise(Nmean = mean(N), Nmedian = median(N), Nerror = mean(N)-median(N))

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
  summarise(Smean = mean(S), Smedian = median(S), Serror = mean(S)-median(S))

## Export the intermediate results
G6_NS <- cbind(N_G6_Ndist, S_G6_Ndist[,4:6])

## Make year column
G6_NS$Year <- as.numeric(str_sub(gsub("_NA", "", G6_NS$compareStrain),-4,-1))

write.table(G6_NS,file="Data/Genetic Distances/GeneCompete/G6_NS.txt")


