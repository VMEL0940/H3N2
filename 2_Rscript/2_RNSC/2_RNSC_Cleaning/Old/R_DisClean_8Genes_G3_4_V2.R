### Influ H3N2 GenDist Proj ###
## Package ##
library(tidyverse)

## Set working directory ##
setwd("/01_GenDisFlu/")

########### HA !!!!!! #################

## Data import ##
HA_G3 <- read.table("Data/Genetic Distances/HA distances/HA_G3_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(HA_G3$tree.nr)
Nnse$value <- 0
HA_G3$N3 <- 0

for (i in 1:nrow(HA_G3))
{Nnse$value[i] <- HA_G3 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == HA_G3$tree.nr[i]) %>% 
  .$N
HA_G3$N3[i] <- abs(HA_G3$N[i] - Nnse$value[i])
print(HA_G3$N3[i])
}


## S subs ##

HA_G3$S3 <- 0

for (i in 1:nrow(HA_G3))
{Nnse$value[i] <- HA_G3 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == HA_G3$tree.nr[i]) %>% 
  .$S
HA_G3$S3[i] <- abs(HA_G3$S[i] - Nnse$value[i])
print(HA_G3$S3[i])
}


write.table(HA_G3,file="Data/Genetic Distances/HA distances/HA_G3_distances_3.txt")


########### NA !!!!!! #################

## Data import ##
NA_G3 <- read.table("Data/Genetic Distances/NA distances/NA_G3_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NA_G3$tree.nr)
Nnse$value <- 0
NA_G3$N3 <- 0

for (i in 1:nrow(NA_G3))
{Nnse$value[i] <- NA_G3 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NA_G3$tree.nr[i]) %>% 
  .$N
NA_G3$N3[i] <- abs(NA_G3$N[i] - Nnse$value[i])
print(NA_G3$N3[i])
}


## S subs ##

NA_G3$S3 <- 0

for (i in 1:nrow(NA_G3))
{Nnse$value[i] <- NA_G3 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NA_G3$tree.nr[i]) %>% 
  .$S
NA_G3$S3[i] <- abs(NA_G3$S[i] - Nnse$value[i])
print(NA_G3$S3[i])
}


write.table(NA_G3,file="Data/Genetic Distances/NA distances/NA_G3_distances_3.txt")


########### M!!!!!! #################

## Data import ##
M_G3 <- read.table("Data/Genetic Distances/M distances/M_G3_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(M_G3$tree.nr)
Nnse$value <- 0
M_G3$N3 <- 0

for (i in 1:nrow(M_G3))
{Nnse$value[i] <- M_G3 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == M_G3$tree.nr[i]) %>% 
  .$N
M_G3$N3[i] <- abs(M_G3$N[i] - Nnse$value[i])
print(M_G3$N3[i])
}


## S subs ##

M_G3$S3 <- 0

for (i in 1:nrow(M_G3))
{Nnse$value[i] <- M_G3 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == M_G3$tree.nr[i]) %>% 
  .$S
M_G3$S3[i] <- abs(M_G3$S[i] - Nnse$value[i])
print(M_G3$S3[i])
}


write.table(M_G3,file="Data/Genetic Distances/M distances/M_G3_distances_3.txt")

########### NP!!!!!! #################

## Data import ##
NP_G3 <- read.table("Data/Genetic Distances/NP distances/NP_G3_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NP_G3$tree.nr)
Nnse$value <- 0
NP_G3$N3 <- 0

for (i in 1:nrow(NP_G3))
{Nnse$value[i] <- NP_G3 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NP_G3$tree.nr[i]) %>% 
  .$N
NP_G3$N3[i] <- abs(NP_G3$N[i] - Nnse$value[i])
print(NP_G3$N3[i])
}


## S subs ##

NP_G3$S3 <- 0

for (i in 1:nrow(NP_G3))
{Nnse$value[i] <- NP_G3 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NP_G3$tree.nr[i]) %>% 
  .$S
NP_G3$S3[i] <- abs(NP_G3$S[i] - Nnse$value[i])
print(NP_G3$S3[i])
}


write.table(NP_G3,file="Data/Genetic Distances/NP distances/NP_G3_distances_3.txt")


########### NS!!!!!! #################


## Data import ##
NS_G3 <- read.table("Data/Genetic Distances/NS distances/NS_G3_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NS_G3$tree.nr)
Nnse$value <- 0
NS_G3$N3 <- 0

for (i in 1:nrow(NS_G3))
{Nnse$value[i] <- NS_G3 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NS_G3$tree.nr[i]) %>% 
  .$N
NS_G3$N3[i] <- abs(NS_G3$N[i] - Nnse$value[i])
print(NS_G3$N3[i])
}


## S subs ##

NS_G3$S3 <- 0

for (i in 1:nrow(NS_G3))
{Nnse$value[i] <- NS_G3 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NS_G3$tree.nr[i]) %>% 
  .$S
NS_G3$S3[i] <- abs(NS_G3$S[i] - Nnse$value[i])
print(NS_G3$S3[i])
}


write.table(NS_G3,file="Data/Genetic Distances/NS distances/NS_G3_distances_3.txt")

########### PA!!!!!! #################

## Data import ##
PA_G3 <- read.table("Data/Genetic Distances/PA distances/PA_G3_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(PA_G3$tree.nr)
Nnse$value <- 0
PA_G3$N3 <- 0

for (i in 1:nrow(PA_G3))
{Nnse$value[i] <- PA_G3 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PA_G3$tree.nr[i]) %>% 
  .$N
PA_G3$N3[i] <- abs(PA_G3$N[i] - Nnse$value[i])
print(PA_G3$N3[i])
}


## S subs ##

PA_G3$S3 <- 0

for (i in 1:nrow(PA_G3))
{Nnse$value[i] <- PA_G3 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PA_G3$tree.nr[i]) %>% 
  .$S
PA_G3$S3[i] <- abs(PA_G3$S[i] - Nnse$value[i])
print(PA_G3$S3[i])
}


write.table(PA_G3,file="Data/Genetic Distances/PA distances/PA_G3_distances_3.txt")


########### PB1!!!!!! #################


## Data import ##
PB1_G3 <- read.table("Data/Genetic Distances/PB1 distances/PB1_G3_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(PB1_G3$tree.nr)
Nnse$value <- 0
PB1_G3$N3 <- 0

for (i in 1:nrow(PB1_G3))
{Nnse$value[i] <- PB1_G3 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB1_G3$tree.nr[i]) %>% 
  .$N
PB1_G3$N3[i] <- abs(PB1_G3$N[i] - Nnse$value[i])
print(PB1_G3$N3[i])
}


## S subs ##

PB1_G3$S3 <- 0

for (i in 1:nrow(PB1_G3))
{Nnse$value[i] <- PB1_G3 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB1_G3$tree.nr[i]) %>% 
  .$S
PB1_G3$S3[i] <- abs(PB1_G3$S[i] - Nnse$value[i])
print(PB1_G3$S3[i])
}


write.table(PB1_G3,file="Data/Genetic Distances/PB1 distances/PB1_G3_distances_3.txt")


########### PB2!!!!!! #################


## Data import ##
PB2_G3 <- read.table("Data/Genetic Distances/PB2 distances/PB2_G3_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(PB2_G3$tree.nr)
Nnse$value <- 0
PB2_G3$N3 <- 0

for (i in 1:nrow(PB2_G3))
{Nnse$value[i] <- PB2_G3 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB2_G3$tree.nr[i]) %>% 
  .$N
PB2_G3$N3[i] <- abs(PB2_G3$N[i] - Nnse$value[i])
print(PB2_G3$N3[i])
}


## S subs ##

PB2_G3$S3 <- 0

for (i in 1:nrow(PB2_G3))
{Nnse$value[i] <- PB2_G3 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB2_G3$tree.nr[i]) %>% 
  .$S
PB2_G3$S3[i] <- abs(PB2_G3$S[i] - Nnse$value[i])
print(PB2_G3$S3[i])
}


write.table(PB2_G3,file="Data/Genetic Distances/PB2 distances/PB2_G3_distances_3.txt")


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
N_G3 <- HA_G3[,c(1,4,5,8)]
colnames(N_G3)[4] <- "HA"
N_G3$"NA" <- NA_G3[,8]
N_G3$"NP" <- NP_G3[,8]
N_G3$"M" <- M_G3[,8]
N_G3$"PB1" <- PB1_G3[,8]
N_G3$"PB2" <- PB2_G3[,8]
N_G3$"NS" <- NS_G3[,8]
N_G3$"PA" <- PA_G3[,8]

## Summarize the 901 tree values
N_G3_Ndist <- N_G3 %>% 
  gather("Gene", "N", 4:11) %>% 
  group_by(vaccineStrain, compareStrain, Gene) %>% 
  summarise(Nmean = mean(N), Nmedian = median(N), Nerror = mean(N)-median(N))

########## Realm of Syn ########

## Combine all data with S
S_G3 <- HA_G3[,c(1,4,5,9)]
colnames(S_G3)[4] <- "HA"
S_G3$"NA" <- NA_G3[,9]
S_G3$"NP" <- NP_G3[,9]
S_G3$"M" <- M_G3[,9]
S_G3$"PB1" <- PB1_G3[,9]
S_G3$"PB2" <- PB2_G3[,9]
S_G3$"NS" <- NS_G3[,9]
S_G3$"PA" <- PA_G3[,9]

## Summarize the 901 tree values
S_G3_Ndist <- S_G3 %>% 
  gather("Gene", "S", 4:11) %>% 
  group_by(vaccineStrain, compareStrain, Gene) %>% 
  summarise(Smean = mean(S), Smedian = median(S), Serror = mean(S)-median(S))

## Export the intermediate results
G3_NS <- cbind(N_G3_Ndist, S_G3_Ndist[,4:6])

## Make year column
G3_NS$Year <- as.numeric(str_sub(gsub("_NA", "", G3_NS$compareStrain),-4,-1))

write.table(G3_NS,file="Data/Genetic Distances/GeneCompete/G3_NS_2.txt")


################## G4 ####################

########### HA !!!!!! #################

## Data import ##
HA_G4 <- read.table("Data/Genetic Distances/HA distances/HA_G4_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(HA_G4$tree.nr)
Nnse$value <- 0
HA_G4$N3 <- 0

for (i in 1:nrow(HA_G4))
{Nnse$value[i] <- HA_G4 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == HA_G4$tree.nr[i]) %>% 
  .$N
HA_G4$N3[i] <- abs(HA_G4$N[i] - Nnse$value[i])
print(HA_G4$N3[i])
}


## S subs ##

HA_G4$S3 <- 0

for (i in 1:nrow(HA_G4))
{Nnse$value[i] <- HA_G4 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == HA_G4$tree.nr[i]) %>% 
  .$S
HA_G4$S3[i] <- abs(HA_G4$S[i] - Nnse$value[i])
print(HA_G4$S3[i])
}


write.table(HA_G4,file="Data/Genetic Distances/HA distances/HA_G4_distances_3.txt")


########### NA !!!!!! #################

## Data import ##
NA_G4 <- read.table("Data/Genetic Distances/NA distances/NA_G4_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NA_G4$tree.nr)
Nnse$value <- 0
NA_G4$N3 <- 0

for (i in 1:nrow(NA_G4))
{Nnse$value[i] <- NA_G4 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NA_G4$tree.nr[i]) %>% 
  .$N
NA_G4$N3[i] <- abs(NA_G4$N[i] - Nnse$value[i])
print(NA_G4$N3[i])
}


## S subs ##

NA_G4$S3 <- 0

for (i in 1:nrow(NA_G4))
{Nnse$value[i] <- NA_G4 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NA_G4$tree.nr[i]) %>% 
  .$S
NA_G4$S3[i] <- abs(NA_G4$S[i] - Nnse$value[i])
print(NA_G4$S3[i])
}


write.table(NA_G4,file="Data/Genetic Distances/NA distances/NA_G4_distances_3.txt")


########### M!!!!!! #################

## Data import ##
M_G4 <- read.table("Data/Genetic Distances/M distances/M_G4_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(M_G4$tree.nr)
Nnse$value <- 0
M_G4$N3 <- 0

for (i in 1:nrow(M_G4))
{Nnse$value[i] <- M_G4 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == M_G4$tree.nr[i]) %>% 
  .$N
M_G4$N3[i] <- abs(M_G4$N[i] - Nnse$value[i])
print(M_G4$N3[i])
}


## S subs ##

M_G4$S3 <- 0

for (i in 1:nrow(M_G4))
{Nnse$value[i] <- M_G4 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == M_G4$tree.nr[i]) %>% 
  .$S
M_G4$S3[i] <- abs(M_G4$S[i] - Nnse$value[i])
print(M_G4$S3[i])
}


write.table(M_G4,file="Data/Genetic Distances/M distances/M_G4_distances_3.txt")

########### NP!!!!!! #################

## Data import ##
NP_G4 <- read.table("Data/Genetic Distances/NP distances/NP_G4_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NP_G4$tree.nr)
Nnse$value <- 0
NP_G4$N3 <- 0

for (i in 1:nrow(NP_G4))
{Nnse$value[i] <- NP_G4 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NP_G4$tree.nr[i]) %>% 
  .$N
NP_G4$N3[i] <- abs(NP_G4$N[i] - Nnse$value[i])
print(NP_G4$N3[i])
}


## S subs ##

NP_G4$S3 <- 0

for (i in 1:nrow(NP_G4))
{Nnse$value[i] <- NP_G4 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NP_G4$tree.nr[i]) %>% 
  .$S
NP_G4$S3[i] <- abs(NP_G4$S[i] - Nnse$value[i])
print(NP_G4$S3[i])
}


write.table(NP_G4,file="Data/Genetic Distances/NP distances/NP_G4_distances_3.txt")


########### NS!!!!!! #################


## Data import ##
NS_G4 <- read.table("Data/Genetic Distances/NS distances/NS_G4_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NS_G4$tree.nr)
Nnse$value <- 0
NS_G4$N3 <- 0

for (i in 1:nrow(NS_G4))
{Nnse$value[i] <- NS_G4 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NS_G4$tree.nr[i]) %>% 
  .$N
NS_G4$N3[i] <- abs(NS_G4$N[i] - Nnse$value[i])
print(NS_G4$N3[i])
}


## S subs ##

NS_G4$S3 <- 0

for (i in 1:nrow(NS_G4))
{Nnse$value[i] <- NS_G4 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NS_G4$tree.nr[i]) %>% 
  .$S
NS_G4$S3[i] <- abs(NS_G4$S[i] - Nnse$value[i])
print(NS_G4$S3[i])
}


write.table(NS_G4,file="Data/Genetic Distances/NS distances/NS_G4_distances_3.txt")

########### PA!!!!!! #################

## Data import ##
PA_G4 <- read.table("Data/Genetic Distances/PA distances/PA_G4_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(PA_G4$tree.nr)
Nnse$value <- 0
PA_G4$N3 <- 0

for (i in 1:nrow(PA_G4))
{Nnse$value[i] <- PA_G4 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PA_G4$tree.nr[i]) %>% 
  .$N
PA_G4$N3[i] <- abs(PA_G4$N[i] - Nnse$value[i])
print(PA_G4$N3[i])
}


## S subs ##

PA_G4$S3 <- 0

for (i in 1:nrow(PA_G4))
{Nnse$value[i] <- PA_G4 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PA_G4$tree.nr[i]) %>% 
  .$S
PA_G4$S3[i] <- abs(PA_G4$S[i] - Nnse$value[i])
print(PA_G4$S3[i])
}


write.table(PA_G4,file="Data/Genetic Distances/PA distances/PA_G4_distances_3.txt")


########### PB1!!!!!! #################


## Data import ##
PB1_G4 <- read.table("Data/Genetic Distances/PB1 distances/PB1_G4_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(PB1_G4$tree.nr)
Nnse$value <- 0
PB1_G4$N3 <- 0

for (i in 1:nrow(PB1_G4))
{Nnse$value[i] <- PB1_G4 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB1_G4$tree.nr[i]) %>% 
  .$N
PB1_G4$N3[i] <- abs(PB1_G4$N[i] - Nnse$value[i])
print(PB1_G4$N3[i])
}


## S subs ##

PB1_G4$S3 <- 0

for (i in 1:nrow(PB1_G4))
{Nnse$value[i] <- PB1_G4 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB1_G4$tree.nr[i]) %>% 
  .$S
PB1_G4$S3[i] <- abs(PB1_G4$S[i] - Nnse$value[i])
print(PB1_G4$S3[i])
}


write.table(PB1_G4,file="Data/Genetic Distances/PB1 distances/PB1_G4_distances_3.txt")


########### PB2!!!!!! #################


## Data import ##
PB2_G4 <- read.table("Data/Genetic Distances/PB2 distances/PB2_G4_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(PB2_G4$tree.nr)
Nnse$value <- 0
PB2_G4$N3 <- 0

for (i in 1:nrow(PB2_G4))
{Nnse$value[i] <- PB2_G4 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB2_G4$tree.nr[i]) %>% 
  .$N
PB2_G4$N3[i] <- abs(PB2_G4$N[i] - Nnse$value[i])
print(PB2_G4$N3[i])
}


## S subs ##

PB2_G4$S3 <- 0

for (i in 1:nrow(PB2_G4))
{Nnse$value[i] <- PB2_G4 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB2_G4$tree.nr[i]) %>% 
  .$S
PB2_G4$S3[i] <- abs(PB2_G4$S[i] - Nnse$value[i])
print(PB2_G4$S3[i])
}


write.table(PB2_G4,file="Data/Genetic Distances/PB2 distances/PB2_G4_distances_3.txt")


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
N_G4 <- HA_G4[,c(1,4,5,8)]
colnames(N_G4)[4] <- "HA"
N_G4$"NA" <- NA_G4[,8]
N_G4$"NP" <- NP_G4[,8]
N_G4$"M" <- M_G4[,8]
N_G4$"PB1" <- PB1_G4[,8]
N_G4$"PB2" <- PB2_G4[,8]
N_G4$"NS" <- NS_G4[,8]
N_G4$"PA" <- PA_G4[,8]

## Summarize the 901 tree values
N_G4_Ndist <- N_G4 %>% 
  gather("Gene", "N", 4:11) %>% 
  group_by(vaccineStrain, compareStrain, Gene) %>% 
  summarise(Nmean = mean(N), Nmedian = median(N), Nerror = mean(N)-median(N))

########## Realm of Syn ########

## Combine all data with S
S_G4 <- HA_G4[,c(1,4,5,9)]
colnames(S_G4)[4] <- "HA"
S_G4$"NA" <- NA_G4[,9]
S_G4$"NP" <- NP_G4[,9]
S_G4$"M" <- M_G4[,9]
S_G4$"PB1" <- PB1_G4[,9]
S_G4$"PB2" <- PB2_G4[,9]
S_G4$"NS" <- NS_G4[,9]
S_G4$"PA" <- PA_G4[,9]

## Summarize the 901 tree values
S_G4_Ndist <- S_G4 %>% 
  gather("Gene", "S", 4:11) %>% 
  group_by(vaccineStrain, compareStrain, Gene) %>% 
  summarise(Smean = mean(S), Smedian = median(S), Serror = mean(S)-median(S))

## Export the intermediate results
G4_NS <- cbind(N_G4_Ndist, S_G4_Ndist[,4:6])

## Make year column
G4_NS$Year <- as.numeric(str_sub(gsub("_NA", "", G4_NS$compareStrain),-4,-1))

write.table(G4_NS,file="Data/Genetic Distances/GeneCompete/G4_NS_2.txt")






