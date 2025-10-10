### Influ H3N2 GenDist Proj ###
## Package ##
library(tidyverse)

## Set working directory ##
setwd("01_GenDisFlu/Data/Genetic Distances/Verification_702/")


######## G1 ###########

## Data import ##
NP_G1 <- read.table("NP/NP_G1_distances.txt", header = T)
## Build redudant column ##
NP_G1 <- NP_G1[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NP_G1$tree.nr)
Nnse$value <- 0
NP_G1$N2 <- 0

for (i in 1:nrow(NP_G1))
{Nnse$value[i] <- NP_G1 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NP_G1$tree.nr[i]) %>% 
  .$N
NP_G1$N2[i] <- abs(NP_G1$N[i] - Nnse$value[i])
print(NP_G1$N2[i])
}

## S subs ##

NP_G1$S2 <- 0

for (i in 1:nrow(NP_G1))
{Nnse$value[i] <- NP_G1 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NP_G1$tree.nr[i]) %>% 
  .$S
NP_G1$S2[i] <- abs(NP_G1$S[i] - Nnse$value[i])
print(NP_G1$S2[i])
}


write.table(NP_G1,file="NP/NP_G1_distances_2.txt")


######## G2 ###########


## Data import ##
NP_G2 <- read.table("NP/NP_G2_distances.txt", header = T)
## Build redudant column ##
NP_G2 <- NP_G2[,1:5]


## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NP_G2$tree.nr)
Nnse$value <- 0
NP_G2$N2 <- 0

for (i in 1:nrow(NP_G2))
{Nnse$value[i] <- NP_G2 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NP_G2$tree.nr[i]) %>% 
  .$N
NP_G2$N2[i] <- abs(NP_G2$N[i] - Nnse$value[i])
print(NP_G2$N2[i])
}


## S subs ##

NP_G2$S2 <- 0

for (i in 1:nrow(NP_G2))
{Nnse$value[i] <- NP_G2 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NP_G2$tree.nr[i]) %>% 
  .$S
NP_G2$S2[i] <- abs(NP_G2$S[i] - Nnse$value[i])
print(NP_G2$S2[i])
}


write.table(NP_G2,file="NP/NP_G2_distances_2.txt")





######## G3 ###########


## Data import ##
NP_G3 <- read.table("NP/NP_G3_distances.txt", header = T)
## Build redudant column ##
NP_G3 <- NP_G3[,1:5]


## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NP_G3$tree.nr)
Nnse$value <- 0
NP_G3$N2 <- 0

for (i in 1:nrow(NP_G3))
{Nnse$value[i] <- NP_G3 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NP_G3$tree.nr[i]) %>% 
  .$N
NP_G3$N2[i] <- abs(NP_G3$N[i] - Nnse$value[i])
print(NP_G3$N2[i])
}


## S subs ##

NP_G3$S2 <- 0

for (i in 1:nrow(NP_G3))
{Nnse$value[i] <- NP_G3 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NP_G3$tree.nr[i]) %>% 
  .$S
NP_G3$S2[i] <- abs(NP_G3$S[i] - Nnse$value[i])
print(NP_G3$S2[i])
}

write.table(NP_G3,file="NP/NP_G3_distances_2.txt")





######## G4 ###########


## Data import ##
NP_G4 <- read.table("NP/NP_G4_distances.txt", header = T)
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
NP_G4$N2[i] <- abs(NP_G4$N[i] - Nnse$value[i])
print(NP_G4$N2[i])
}

## S subs ##

NP_G4$S2 <- 0

for (i in 1:nrow(NP_G4))
{Nnse$value[i] <- NP_G4 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NP_G4$tree.nr[i]) %>% 
  .$S
NP_G4$S2[i] <- abs(NP_G4$S[i] - Nnse$value[i])
print(NP_G4$S2[i])
}


write.table(NP_G4,file="NP/NP_G4_distances_2.txt")



######## G5 ###########


## Data import ##
NP_G5 <- read.table("NP/NP_G5_distances.txt", header = T)
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
NP_G5$N2[i] <- abs(NP_G5$N[i] - Nnse$value[i])
print(NP_G5$N2[i])
}

## S subs ##

NP_G5$S2 <- 0

for (i in 1:nrow(NP_G5))
{Nnse$value[i] <- NP_G5 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NP_G5$tree.nr[i]) %>% 
  .$S
NP_G5$S2[i] <- abs(NP_G5$S[i] - Nnse$value[i])
print(NP_G5$S2[i])
}


write.table(NP_G5,file="NP/NP_G5_distances_2.txt")



######## G6 ###########


## Data import ##
NP_G6 <- read.table("NP/NP_G6_distances.txt", header = T)
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
NP_G6$N2[i] <- abs(NP_G6$N[i] - Nnse$value[i])
print(NP_G6$N2[i])
}

## S subs ##

NP_G6$S2 <- 0

for (i in 1:nrow(NP_G6))
{Nnse$value[i] <- NP_G6 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NP_G6$tree.nr[i]) %>% 
  .$S
NP_G6$S2[i] <- abs(NP_G6$S[i] - Nnse$value[i])
print(NP_G6$S2[i])
}


write.table(NP_G6,file="NP/NP_G6_distances_2.txt")



######## G7 ###########


## Data import ##
NP_G7 <- read.table("NP/NP_G7_distances.txt", header = T)
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
NP_G7$N2[i] <- abs(NP_G7$N[i] - Nnse$value[i])
print(NP_G7$N2[i])
}

## S subs ##

NP_G7$S2 <- 0

for (i in 1:nrow(NP_G7))
{Nnse$value[i] <- NP_G7 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NP_G7$tree.nr[i]) %>% 
  .$S
NP_G7$S2[i] <- abs(NP_G7$S[i] - Nnse$value[i])
print(NP_G7$S2[i])
}


write.table(NP_G7,file="NP/NP_G7_distances_2.txt")




######## G8 ###########


## Data import ##
NP_G8 <- read.table("NP/NP_G8_distances.txt", header = T)
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
NP_G8$N2[i] <- abs(NP_G8$N[i] - Nnse$value[i])
print(NP_G8$N2[i])
}

## S subs ##

NP_G8$S2 <- 0

for (i in 1:nrow(NP_G8))
{Nnse$value[i] <- NP_G8 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NP_G8$tree.nr[i]) %>% 
  .$S
NP_G8$S2[i] <- abs(NP_G8$S[i] - Nnse$value[i])
print(NP_G8$S2[i])
}


write.table(NP_G8,file="NP/NP_G8_distances_2.txt")




######## G9 ###########


## Data import ##
NP_G9 <- read.table("NP/NP_G9_distances.txt", header = T)
## Build redudant column ##
NP_G9 <- NP_G9[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NP_G9$tree.nr)
Nnse$value <- 0
NP_G9$N2 <- 0

for (i in 1:nrow(NP_G9))
{Nnse$value[i] <- NP_G9 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NP_G9$tree.nr[i]) %>% 
  .$N
NP_G9$N2[i] <- abs(NP_G9$N[i] - Nnse$value[i])
print(NP_G9$N2[i])
}

## S subs ##

NP_G9$S2 <- 0

for (i in 1:nrow(NP_G9))
{Nnse$value[i] <- NP_G9 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NP_G9$tree.nr[i]) %>% 
  .$S
NP_G9$S2[i] <- abs(NP_G9$S[i] - Nnse$value[i])
print(NP_G9$S2[i])
}


write.table(NP_G9,file="NP/NP_G9_distances_2.txt")




######## G10 ###########


## Data import ##
NP_G10 <- read.table("NP/NP_G10_distances.txt", header = T)
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
NP_G10$N2[i] <- abs(NP_G10$N[i] - Nnse$value[i])
print(NP_G10$N2[i])
}

## S subs ##

NP_G10$S2 <- 0

for (i in 1:nrow(NP_G10))
{Nnse$value[i] <- NP_G10 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NP_G10$tree.nr[i]) %>% 
  .$S
NP_G10$S2[i] <- abs(NP_G10$S[i] - Nnse$value[i])
print(NP_G10$S2[i])
}


write.table(NP_G10,file="NP/NP_G10_distances_2.txt")




######## G11 ###########


## Data import ##
NP_G11 <- read.table("NP/NP_G11_distances.txt", header = T)
## Build redudant column ##
NP_G11 <- NP_G11[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NP_G11$tree.nr)
Nnse$value <- 0
NP_G11$N2 <- 0

for (i in 1:nrow(NP_G11))
{Nnse$value[i] <- NP_G11 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NP_G11$tree.nr[i]) %>% 
  .$N
NP_G11$N2[i] <- abs(NP_G11$N[i] - Nnse$value[i])
print(NP_G11$N2[i])
}

## S subs ##

NP_G11$S2 <- 0

for (i in 1:nrow(NP_G11))
{Nnse$value[i] <- NP_G11 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NP_G11$tree.nr[i]) %>% 
  .$S
NP_G11$S2[i] <- abs(NP_G11$S[i] - Nnse$value[i])
print(NP_G11$S2[i])
}


write.table(NP_G11,file="NP/NP_G11_distances_2.txt")




######## G12 ###########


## Data import ##
NP_G12 <- read.table("NP/NP_G12_distances.txt", header = T)
## Build redudant column ##
NP_G12 <- NP_G12[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NP_G12$tree.nr)
Nnse$value <- 0
NP_G12$N2 <- 0

for (i in 1:nrow(NP_G12))
{Nnse$value[i] <- NP_G12 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NP_G12$tree.nr[i]) %>% 
  .$N
NP_G12$N2[i] <- abs(NP_G12$N[i] - Nnse$value[i])
print(NP_G12$N2[i])
}

## S subs ##

NP_G12$S2 <- 0

for (i in 1:nrow(NP_G12))
{Nnse$value[i] <- NP_G12 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NP_G12$tree.nr[i]) %>% 
  .$S
NP_G12$S2[i] <- abs(NP_G12$S[i] - Nnse$value[i])
print(NP_G12$S2[i])
}


write.table(NP_G12,file="NP/NP_G12_distances_2.txt")




######## G13 ###########


## Data import ##
NP_G13 <- read.table("NP/NP_G13_distances.txt", header = T)
## Build redudant column ##
NP_G13 <- NP_G13[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NP_G13$tree.nr)
Nnse$value <- 0
NP_G13$N2 <- 0

for (i in 1:nrow(NP_G13))
{Nnse$value[i] <- NP_G13 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NP_G13$tree.nr[i]) %>% 
  .$N
NP_G13$N2[i] <- abs(NP_G13$N[i] - Nnse$value[i])
print(NP_G13$N2[i])
}

## S subs ##

NP_G13$S2 <- 0

for (i in 1:nrow(NP_G13))
{Nnse$value[i] <- NP_G13 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NP_G13$tree.nr[i]) %>% 
  .$S
NP_G13$S2[i] <- abs(NP_G13$S[i] - Nnse$value[i])
print(NP_G13$S2[i])
}


write.table(NP_G13,file="NP/NP_G13_distances_2.txt")
