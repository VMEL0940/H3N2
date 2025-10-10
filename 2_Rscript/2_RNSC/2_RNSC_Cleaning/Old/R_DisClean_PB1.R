### Influ H3N2 GenDist Proj ###
## Package ##
library(tidyverse)

## Set working directory ##
setwd("01_GenDisFlu/Data/Genetic Distances/Verification_702/")

######## G1 ###########

## Data import ##
PB1_G1 <- read.table("PB1/PB1_G1_distances.txt", header = T)
## Build redudant column ##
PB1_G1 <- PB1_G1[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(PB1_G1$tree.nr)
Nnse$value <- 0
PB1_G1$N2 <- 0

for (i in 1:nrow(PB1_G1))
{Nnse$value[i] <- PB1_G1 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB1_G1$tree.nr[i]) %>% 
  .$N
PB1_G1$N2[i] <- abs(PB1_G1$N[i] - Nnse$value[i])
print(PB1_G1$N2[i])
}

## S subs ##

PB1_G1$S2 <- 0

for (i in 1:nrow(PB1_G1))
{Nnse$value[i] <- PB1_G1 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB1_G1$tree.nr[i]) %>% 
  .$S
PB1_G1$S2[i] <- abs(PB1_G1$S[i] - Nnse$value[i])
print(PB1_G1$S2[i])
}


write.table(PB1_G1,file="PB1/PB1_G1_distances_2.txt")


######## G2 ###########


## Data import ##
PB1_G2 <- read.table("PB1/PB1_G2_distances.txt", header = T)
## Build redudant column ##
PB1_G2 <- PB1_G2[,1:5]


## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(PB1_G2$tree.nr)
Nnse$value <- 0
PB1_G2$N2 <- 0

for (i in 1:nrow(PB1_G2))
{Nnse$value[i] <- PB1_G2 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB1_G2$tree.nr[i]) %>% 
  .$N
PB1_G2$N2[i] <- abs(PB1_G2$N[i] - Nnse$value[i])
print(PB1_G2$N2[i])
}


## S subs ##

PB1_G2$S2 <- 0

for (i in 1:nrow(PB1_G2))
{Nnse$value[i] <- PB1_G2 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB1_G2$tree.nr[i]) %>% 
  .$S
PB1_G2$S2[i] <- abs(PB1_G2$S[i] - Nnse$value[i])
print(PB1_G2$S2[i])
}


write.table(PB1_G2,file="PB1/PB1_G2_distances_2.txt")





######## G3 ###########


## Data import ##
PB1_G3 <- read.table("PB1/PB1_G3_distances.txt", header = T)
## Build redudant column ##
PB1_G3 <- PB1_G3[,1:5]


## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(PB1_G3$tree.nr)
Nnse$value <- 0
PB1_G3$N2 <- 0

for (i in 1:nrow(PB1_G3))
{Nnse$value[i] <- PB1_G3 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB1_G3$tree.nr[i]) %>% 
  .$N
PB1_G3$N2[i] <- abs(PB1_G3$N[i] - Nnse$value[i])
print(PB1_G3$N2[i])
}


## S subs ##

PB1_G3$S2 <- 0

for (i in 1:nrow(PB1_G3))
{Nnse$value[i] <- PB1_G3 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB1_G3$tree.nr[i]) %>% 
  .$S
PB1_G3$S2[i] <- abs(PB1_G3$S[i] - Nnse$value[i])
print(PB1_G3$S2[i])
}

write.table(PB1_G3,file="PB1/PB1_G3_distances_2.txt")





######## G4 ###########


## Data import ##
PB1_G4 <- read.table("PB1/PB1_G4_distances.txt", header = T)
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
PB1_G4$N2[i] <- abs(PB1_G4$N[i] - Nnse$value[i])
print(PB1_G4$N2[i])
}

## S subs ##

PB1_G4$S2 <- 0

for (i in 1:nrow(PB1_G4))
{Nnse$value[i] <- PB1_G4 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB1_G4$tree.nr[i]) %>% 
  .$S
PB1_G4$S2[i] <- abs(PB1_G4$S[i] - Nnse$value[i])
print(PB1_G4$S2[i])
}


write.table(PB1_G4,file="PB1/PB1_G4_distances_2.txt")



######## G5 ###########


## Data import ##
PB1_G5 <- read.table("PB1/PB1_G5_distances.txt", header = T)
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
PB1_G5$N2[i] <- abs(PB1_G5$N[i] - Nnse$value[i])
print(PB1_G5$N2[i])
}

## S subs ##

PB1_G5$S2 <- 0

for (i in 1:nrow(PB1_G5))
{Nnse$value[i] <- PB1_G5 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB1_G5$tree.nr[i]) %>% 
  .$S
PB1_G5$S2[i] <- abs(PB1_G5$S[i] - Nnse$value[i])
print(PB1_G5$S2[i])
}


write.table(PB1_G5,file="PB1/PB1_G5_distances_2.txt")



######## G6 ###########


## Data import ##
PB1_G6 <- read.table("PB1/PB1_G6_distances.txt", header = T)
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
PB1_G6$N2[i] <- abs(PB1_G6$N[i] - Nnse$value[i])
print(PB1_G6$N2[i])
}

## S subs ##

PB1_G6$S2 <- 0

for (i in 1:nrow(PB1_G6))
{Nnse$value[i] <- PB1_G6 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB1_G6$tree.nr[i]) %>% 
  .$S
PB1_G6$S2[i] <- abs(PB1_G6$S[i] - Nnse$value[i])
print(PB1_G6$S2[i])
}


write.table(PB1_G6,file="PB1/PB1_G6_distances_2.txt")



######## G7 ###########


## Data import ##
PB1_G7 <- read.table("PB1/PB1_G7_distances.txt", header = T)
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
PB1_G7$N2[i] <- abs(PB1_G7$N[i] - Nnse$value[i])
print(PB1_G7$N2[i])
}

## S subs ##

PB1_G7$S2 <- 0

for (i in 1:nrow(PB1_G7))
{Nnse$value[i] <- PB1_G7 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB1_G7$tree.nr[i]) %>% 
  .$S
PB1_G7$S2[i] <- abs(PB1_G7$S[i] - Nnse$value[i])
print(PB1_G7$S2[i])
}


write.table(PB1_G7,file="PB1/PB1_G7_distances_2.txt")




######## G8 ###########


## Data import ##
PB1_G8 <- read.table("PB1/PB1_G8_distances.txt", header = T)
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
PB1_G8$N2[i] <- abs(PB1_G8$N[i] - Nnse$value[i])
print(PB1_G8$N2[i])
}

## S subs ##

PB1_G8$S2 <- 0

for (i in 1:nrow(PB1_G8))
{Nnse$value[i] <- PB1_G8 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB1_G8$tree.nr[i]) %>% 
  .$S
PB1_G8$S2[i] <- abs(PB1_G8$S[i] - Nnse$value[i])
print(PB1_G8$S2[i])
}


write.table(PB1_G8,file="PB1/PB1_G8_distances_2.txt")




######## G9 ###########


## Data import ##
PB1_G9 <- read.table("PB1/PB1_G9_distances.txt", header = T)
## Build redudant column ##
PB1_G9 <- PB1_G9[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(PB1_G9$tree.nr)
Nnse$value <- 0
PB1_G9$N2 <- 0

for (i in 1:nrow(PB1_G9))
{Nnse$value[i] <- PB1_G9 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB1_G9$tree.nr[i]) %>% 
  .$N
PB1_G9$N2[i] <- abs(PB1_G9$N[i] - Nnse$value[i])
print(PB1_G9$N2[i])
}

## S subs ##

PB1_G9$S2 <- 0

for (i in 1:nrow(PB1_G9))
{Nnse$value[i] <- PB1_G9 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB1_G9$tree.nr[i]) %>% 
  .$S
PB1_G9$S2[i] <- abs(PB1_G9$S[i] - Nnse$value[i])
print(PB1_G9$S2[i])
}


write.table(PB1_G9,file="PB1/PB1_G9_distances_2.txt")




######## G10 ###########


## Data import ##
PB1_G10 <- read.table("PB1/PB1_G10_distances.txt", header = T)
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
PB1_G10$N2[i] <- abs(PB1_G10$N[i] - Nnse$value[i])
print(PB1_G10$N2[i])
}

## S subs ##

PB1_G10$S2 <- 0

for (i in 1:nrow(PB1_G10))
{Nnse$value[i] <- PB1_G10 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB1_G10$tree.nr[i]) %>% 
  .$S
PB1_G10$S2[i] <- abs(PB1_G10$S[i] - Nnse$value[i])
print(PB1_G10$S2[i])
}


write.table(PB1_G10,file="PB1/PB1_G10_distances_2.txt")




######## G11 ###########


## Data import ##
PB1_G11 <- read.table("PB1/PB1_G11_distances.txt", header = T)
## Build redudant column ##
PB1_G11 <- PB1_G11[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(PB1_G11$tree.nr)
Nnse$value <- 0
PB1_G11$N2 <- 0

for (i in 1:nrow(PB1_G11))
{Nnse$value[i] <- PB1_G11 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB1_G11$tree.nr[i]) %>% 
  .$N
PB1_G11$N2[i] <- abs(PB1_G11$N[i] - Nnse$value[i])
print(PB1_G11$N2[i])
}

## S subs ##

PB1_G11$S2 <- 0

for (i in 1:nrow(PB1_G11))
{Nnse$value[i] <- PB1_G11 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB1_G11$tree.nr[i]) %>% 
  .$S
PB1_G11$S2[i] <- abs(PB1_G11$S[i] - Nnse$value[i])
print(PB1_G11$S2[i])
}


write.table(PB1_G11,file="PB1/PB1_G11_distances_2.txt")




######## G12 ###########


## Data import ##
PB1_G12 <- read.table("PB1/PB1_G12_distances.txt", header = T)
## Build redudant column ##
PB1_G12 <- PB1_G12[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(PB1_G12$tree.nr)
Nnse$value <- 0
PB1_G12$N2 <- 0

for (i in 1:nrow(PB1_G12))
{Nnse$value[i] <- PB1_G12 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB1_G12$tree.nr[i]) %>% 
  .$N
PB1_G12$N2[i] <- abs(PB1_G12$N[i] - Nnse$value[i])
print(PB1_G12$N2[i])
}

## S subs ##

PB1_G12$S2 <- 0

for (i in 1:nrow(PB1_G12))
{Nnse$value[i] <- PB1_G12 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB1_G12$tree.nr[i]) %>% 
  .$S
PB1_G12$S2[i] <- abs(PB1_G12$S[i] - Nnse$value[i])
print(PB1_G12$S2[i])
}


write.table(PB1_G12,file="PB1/PB1_G12_distances_2.txt")




######## G13 ###########


## Data import ##
PB1_G13 <- read.table("PB1/PB1_G13_distances.txt", header = T)
## Build redudant column ##
PB1_G13 <- PB1_G13[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(PB1_G13$tree.nr)
Nnse$value <- 0
PB1_G13$N2 <- 0

for (i in 1:nrow(PB1_G13))
{Nnse$value[i] <- PB1_G13 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB1_G13$tree.nr[i]) %>% 
  .$N
PB1_G13$N2[i] <- abs(PB1_G13$N[i] - Nnse$value[i])
print(PB1_G13$N2[i])
}

## S subs ##

PB1_G13$S2 <- 0

for (i in 1:nrow(PB1_G13))
{Nnse$value[i] <- PB1_G13 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PB1_G13$tree.nr[i]) %>% 
  .$S
PB1_G13$S2[i] <- abs(PB1_G13$S[i] - Nnse$value[i])
print(PB1_G13$S2[i])
}


write.table(PB1_G13,file="PB1/PB1_G13_distances_2.txt")
