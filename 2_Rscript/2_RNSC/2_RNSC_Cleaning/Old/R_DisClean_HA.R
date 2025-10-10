### Influ H3N2 GenDist Proj ###
## Package ##
library(tidyverse)

## Set working directory ##
setwd("01_GenDisFlu/Data/Genetic Distances/Verification_702/")


######## G1 ###########

## Data import ##
HA_G1 <- read.table("HA/HA_G1_distances.txt", header = T)
## Build redudant column ##
HA_G1 <- HA_G1[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(HA_G1$tree.nr)
Nnse$value <- 0
HA_G1$N2 <- 0

for (i in 1:nrow(HA_G1))
  {Nnse$value[i] <- HA_G1 %>% 
     filter(vaccineStrain == compareStrain) %>% 
     filter(tree.nr == HA_G1$tree.nr[i]) %>% 
    .$N
  HA_G1$N2[i] <- abs(HA_G1$N[i] - Nnse$value[i])
  print(HA_G1$N2[i])
}

## S subs ##

HA_G1$S2 <- 0

for (i in 1:nrow(HA_G1))
{Nnse$value[i] <- HA_G1 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == HA_G1$tree.nr[i]) %>% 
  .$S
HA_G1$S2[i] <- abs(HA_G1$S[i] - Nnse$value[i])
print(HA_G1$S2[i])
}


write.table(HA_G1,file="HA/HA_G1_distances_2.txt")


######## G2 ###########


## Data import ##
HA_G2 <- read.table("HA/HA_G2_distances.txt", header = T)
## Build redudant column ##
HA_G2 <- HA_G2[,1:5]


## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(HA_G2$tree.nr)
Nnse$value <- 0
HA_G2$N2 <- 0

for (i in 1:nrow(HA_G2))
{Nnse$value[i] <- HA_G2 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == HA_G2$tree.nr[i]) %>% 
  .$N
HA_G2$N2[i] <- abs(HA_G2$N[i] - Nnse$value[i])
print(HA_G2$N2[i])
}


## S subs ##

HA_G2$S2 <- 0

for (i in 1:nrow(HA_G2))
{Nnse$value[i] <- HA_G2 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == HA_G2$tree.nr[i]) %>% 
  .$S
HA_G2$S2[i] <- abs(HA_G2$S[i] - Nnse$value[i])
print(HA_G2$S2[i])
}


write.table(HA_G2,file="HA/HA_G2_distances_2.txt")





######## G3 ###########


## Data import ##
HA_G3 <- read.table("HA/HA_G3_distances.txt", header = T)
## Build redudant column ##
HA_G3 <- HA_G3[,1:5]


## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(HA_G3$tree.nr)
Nnse$value <- 0
HA_G3$N2 <- 0

for (i in 1:nrow(HA_G3))
{Nnse$value[i] <- HA_G3 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == HA_G3$tree.nr[i]) %>% 
  .$N
HA_G3$N2[i] <- abs(HA_G3$N[i] - Nnse$value[i])
print(HA_G3$N2[i])
}


## S subs ##

HA_G3$S2 <- 0

for (i in 1:nrow(HA_G3))
{Nnse$value[i] <- HA_G3 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == HA_G3$tree.nr[i]) %>% 
  .$S
HA_G3$S2[i] <- abs(HA_G3$S[i] - Nnse$value[i])
print(HA_G3$S2[i])
}

write.table(HA_G3,file="HA/HA_G3_distances_2.txt")





######## G4 ###########


## Data import ##
HA_G4 <- read.table("HA/HA_G4_distances.txt", header = T)
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
HA_G4$N2[i] <- abs(HA_G4$N[i] - Nnse$value[i])
print(HA_G4$N2[i])
}

## S subs ##

HA_G4$S2 <- 0

for (i in 1:nrow(HA_G4))
{Nnse$value[i] <- HA_G4 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == HA_G4$tree.nr[i]) %>% 
  .$S
HA_G4$S2[i] <- abs(HA_G4$S[i] - Nnse$value[i])
print(HA_G4$S2[i])
}


write.table(HA_G4,file="HA/HA_G4_distances_2.txt")



######## G5 ###########


## Data import ##
HA_G5 <- read.table("HA/HA_G5_distances.txt", header = T)
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
HA_G5$N2[i] <- abs(HA_G5$N[i] - Nnse$value[i])
print(HA_G5$N2[i])
}

## S subs ##

HA_G5$S2 <- 0

for (i in 1:nrow(HA_G5))
{Nnse$value[i] <- HA_G5 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == HA_G5$tree.nr[i]) %>% 
  .$S
HA_G5$S2[i] <- abs(HA_G5$S[i] - Nnse$value[i])
print(HA_G5$S2[i])
}


write.table(HA_G5,file="HA/HA_G5_distances_2.txt")



######## G6 ###########


## Data import ##
HA_G6 <- read.table("HA/HA_G6_distances.txt", header = T)
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
HA_G6$N2[i] <- abs(HA_G6$N[i] - Nnse$value[i])
print(HA_G6$N2[i])
}

## S subs ##

HA_G6$S2 <- 0

for (i in 1:nrow(HA_G6))
{Nnse$value[i] <- HA_G6 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == HA_G6$tree.nr[i]) %>% 
  .$S
HA_G6$S2[i] <- abs(HA_G6$S[i] - Nnse$value[i])
print(HA_G6$S2[i])
}


write.table(HA_G6,file="HA/HA_G6_distances_2.txt")



######## G7 ###########


## Data import ##
HA_G7 <- read.table("HA/HA_G7_distances.txt", header = T)
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
HA_G7$N2[i] <- abs(HA_G7$N[i] - Nnse$value[i])
print(HA_G7$N2[i])
}

## S subs ##

HA_G7$S2 <- 0

for (i in 1:nrow(HA_G7))
{Nnse$value[i] <- HA_G7 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == HA_G7$tree.nr[i]) %>% 
  .$S
HA_G7$S2[i] <- abs(HA_G7$S[i] - Nnse$value[i])
print(HA_G7$S2[i])
}


write.table(HA_G7,file="HA/HA_G7_distances_2.txt")




######## G8 ###########


## Data import ##
HA_G8 <- read.table("HA/HA_G8_distances.txt", header = T)
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
HA_G8$N2[i] <- abs(HA_G8$N[i] - Nnse$value[i])
print(HA_G8$N2[i])
}

## S subs ##

HA_G8$S2 <- 0

for (i in 1:nrow(HA_G8))
{Nnse$value[i] <- HA_G8 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == HA_G8$tree.nr[i]) %>% 
  .$S
HA_G8$S2[i] <- abs(HA_G8$S[i] - Nnse$value[i])
print(HA_G8$S2[i])
}


write.table(HA_G8,file="HA/HA_G8_distances_2.txt")




######## G9 ###########


## Data import ##
HA_G9 <- read.table("HA/HA_G9_distances.txt", header = T)
## Build redudant column ##
HA_G9 <- HA_G9[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(HA_G9$tree.nr)
Nnse$value <- 0
HA_G9$N2 <- 0

for (i in 1:nrow(HA_G9))
{Nnse$value[i] <- HA_G9 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == HA_G9$tree.nr[i]) %>% 
  .$N
HA_G9$N2[i] <- abs(HA_G9$N[i] - Nnse$value[i])
print(HA_G9$N2[i])
}

## S subs ##

HA_G9$S2 <- 0

for (i in 1:nrow(HA_G9))
{Nnse$value[i] <- HA_G9 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == HA_G9$tree.nr[i]) %>% 
  .$S
HA_G9$S2[i] <- abs(HA_G9$S[i] - Nnse$value[i])
print(HA_G9$S2[i])
}


write.table(HA_G9,file="HA/HA_G9_distances_2.txt")




######## G10 ###########


## Data import ##
HA_G10 <- read.table("HA/HA_G10_distances.txt", header = T)
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
HA_G10$N2[i] <- abs(HA_G10$N[i] - Nnse$value[i])
print(HA_G10$N2[i])
}

## S subs ##

HA_G10$S2 <- 0

for (i in 1:nrow(HA_G10))
{Nnse$value[i] <- HA_G10 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == HA_G10$tree.nr[i]) %>% 
  .$S
HA_G10$S2[i] <- abs(HA_G10$S[i] - Nnse$value[i])
print(HA_G10$S2[i])
}


write.table(HA_G10,file="HA/HA_G10_distances_2.txt")




######## G11 ###########


## Data import ##
HA_G11 <- read.table("HA/HA_G11_distances.txt", header = T)
## Build redudant column ##
HA_G11 <- HA_G11[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(HA_G11$tree.nr)
Nnse$value <- 0
HA_G11$N2 <- 0

for (i in 1:nrow(HA_G11))
{Nnse$value[i] <- HA_G11 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == HA_G11$tree.nr[i]) %>% 
  .$N
HA_G11$N2[i] <- abs(HA_G11$N[i] - Nnse$value[i])
print(HA_G11$N2[i])
}

## S subs ##

HA_G11$S2 <- 0

for (i in 1:nrow(HA_G11))
{Nnse$value[i] <- HA_G11 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == HA_G11$tree.nr[i]) %>% 
  .$S
HA_G11$S2[i] <- abs(HA_G11$S[i] - Nnse$value[i])
print(HA_G11$S2[i])
}


write.table(HA_G11,file="HA/HA_G11_distances_2.txt")




######## G12 ###########


## Data import ##
HA_G12 <- read.table("HA/HA_G12_distances.txt", header = T)
## Build redudant column ##
HA_G12 <- HA_G12[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(HA_G12$tree.nr)
Nnse$value <- 0
HA_G12$N2 <- 0

for (i in 1:nrow(HA_G12))
{Nnse$value[i] <- HA_G12 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == HA_G12$tree.nr[i]) %>% 
  .$N
HA_G12$N2[i] <- abs(HA_G12$N[i] - Nnse$value[i])
print(HA_G12$N2[i])
}

## S subs ##

HA_G12$S2 <- 0

for (i in 1:nrow(HA_G12))
{Nnse$value[i] <- HA_G12 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == HA_G12$tree.nr[i]) %>% 
  .$S
HA_G12$S2[i] <- abs(HA_G12$S[i] - Nnse$value[i])
print(HA_G12$S2[i])
}


write.table(HA_G12,file="HA/HA_G12_distances_2.txt")




######## G13 ###########


## Data import ##
HA_G13 <- read.table("HA/HA_G13_distances.txt", header = T)
## Build redudant column ##
HA_G13 <- HA_G13[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(HA_G13$tree.nr)
Nnse$value <- 0
HA_G13$N2 <- 0

for (i in 1:nrow(HA_G13))
{Nnse$value[i] <- HA_G13 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == HA_G13$tree.nr[i]) %>% 
  .$N
HA_G13$N2[i] <- abs(HA_G13$N[i] - Nnse$value[i])
print(HA_G13$N2[i])
}

## S subs ##

HA_G13$S2 <- 0

for (i in 1:nrow(HA_G13))
{Nnse$value[i] <- HA_G13 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == HA_G13$tree.nr[i]) %>% 
  .$S
HA_G13$S2[i] <- abs(HA_G13$S[i] - Nnse$value[i])
print(HA_G13$S2[i])
}


write.table(HA_G13,file="HA/HA_G13_distances_2.txt")
