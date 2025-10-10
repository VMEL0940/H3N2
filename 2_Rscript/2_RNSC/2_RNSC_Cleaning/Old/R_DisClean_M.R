### Influ H3N2 GenDist Proj ###
## Package ##
library(tidyverse)

## Set working directory ##
setwd("01_GenDisFlu/Data/Genetic Distances/Verification_702/")

######## G1 ###########

## Data import ##
M_G1 <- read.table("M/M_G1_distances.txt", header = T)
## Build redudant column ##
M_G1 <- M_G1[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(M_G1$tree.nr)
Nnse$value <- 0
M_G1$N2 <- 0

for (i in 1:nrow(M_G1))
{Nnse$value[i] <- M_G1 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == M_G1$tree.nr[i]) %>% 
  .$N
M_G1$N2[i] <- abs(M_G1$N[i] - Nnse$value[i])
print(M_G1$N2[i])
}

## S subs ##

M_G1$S2 <- 0

for (i in 1:nrow(M_G1))
{Nnse$value[i] <- M_G1 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == M_G1$tree.nr[i]) %>% 
  .$S
M_G1$S2[i] <- abs(M_G1$S[i] - Nnse$value[i])
print(M_G1$S2[i])
}


write.table(M_G1,file="M/M_G1_distances_2.txt")


######## G2 ###########


## Data import ##
M_G2 <- read.table("M/M_G2_distances.txt", header = T)
## Build redudant column ##
M_G2 <- M_G2[,1:5]


## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(M_G2$tree.nr)
Nnse$value <- 0
M_G2$N2 <- 0

for (i in 1:nrow(M_G2))
{Nnse$value[i] <- M_G2 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == M_G2$tree.nr[i]) %>% 
  .$N
M_G2$N2[i] <- abs(M_G2$N[i] - Nnse$value[i])
print(M_G2$N2[i])
}


## S subs ##

M_G2$S2 <- 0

for (i in 1:nrow(M_G2))
{Nnse$value[i] <- M_G2 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == M_G2$tree.nr[i]) %>% 
  .$S
M_G2$S2[i] <- abs(M_G2$S[i] - Nnse$value[i])
print(M_G2$S2[i])
}


write.table(M_G2,file="M/M_G2_distances_2.txt")





######## G3 ###########


## Data import ##
M_G3 <- read.table("M/M_G3_distances.txt", header = T)
## Build redudant column ##
M_G3 <- M_G3[,1:5]


## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(M_G3$tree.nr)
Nnse$value <- 0
M_G3$N2 <- 0

for (i in 1:nrow(M_G3))
{Nnse$value[i] <- M_G3 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == M_G3$tree.nr[i]) %>% 
  .$N
M_G3$N2[i] <- abs(M_G3$N[i] - Nnse$value[i])
print(M_G3$N2[i])
}


## S subs ##

M_G3$S2 <- 0

for (i in 1:nrow(M_G3))
{Nnse$value[i] <- M_G3 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == M_G3$tree.nr[i]) %>% 
  .$S
M_G3$S2[i] <- abs(M_G3$S[i] - Nnse$value[i])
print(M_G3$S2[i])
}

write.table(M_G3,file="M/M_G3_distances_2.txt")





######## G4 ###########


## Data import ##
M_G4 <- read.table("M/M_G4_distances.txt", header = T)
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
M_G4$N2[i] <- abs(M_G4$N[i] - Nnse$value[i])
print(M_G4$N2[i])
}

## S subs ##

M_G4$S2 <- 0

for (i in 1:nrow(M_G4))
{Nnse$value[i] <- M_G4 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == M_G4$tree.nr[i]) %>% 
  .$S
M_G4$S2[i] <- abs(M_G4$S[i] - Nnse$value[i])
print(M_G4$S2[i])
}


write.table(M_G4,file="M/M_G4_distances_2.txt")



######## G5 ###########


## Data import ##
M_G5 <- read.table("M/M_G5_distances.txt", header = T)
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
M_G5$N2[i] <- abs(M_G5$N[i] - Nnse$value[i])
print(M_G5$N2[i])
}

## S subs ##

M_G5$S2 <- 0

for (i in 1:nrow(M_G5))
{Nnse$value[i] <- M_G5 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == M_G5$tree.nr[i]) %>% 
  .$S
M_G5$S2[i] <- abs(M_G5$S[i] - Nnse$value[i])
print(M_G5$S2[i])
}


write.table(M_G5,file="M/M_G5_distances_2.txt")



######## G6 ###########


## Data import ##
M_G6 <- read.table("M/M_G6_distances.txt", header = T)
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
M_G6$N2[i] <- abs(M_G6$N[i] - Nnse$value[i])
print(M_G6$N2[i])
}

## S subs ##

M_G6$S2 <- 0

for (i in 1:nrow(M_G6))
{Nnse$value[i] <- M_G6 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == M_G6$tree.nr[i]) %>% 
  .$S
M_G6$S2[i] <- abs(M_G6$S[i] - Nnse$value[i])
print(M_G6$S2[i])
}


write.table(M_G6,file="M/M_G6_distances_2.txt")



######## G7 ###########


## Data import ##
M_G7 <- read.table("M/M_G7_distances.txt", header = T)
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
M_G7$N2[i] <- abs(M_G7$N[i] - Nnse$value[i])
print(M_G7$N2[i])
}

## S subs ##

M_G7$S2 <- 0

for (i in 1:nrow(M_G7))
{Nnse$value[i] <- M_G7 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == M_G7$tree.nr[i]) %>% 
  .$S
M_G7$S2[i] <- abs(M_G7$S[i] - Nnse$value[i])
print(M_G7$S2[i])
}


write.table(M_G7,file="M/M_G7_distances_2.txt")




######## G8 ###########


## Data import ##
M_G8 <- read.table("M/M_G8_distances.txt", header = T)
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
M_G8$N2[i] <- abs(M_G8$N[i] - Nnse$value[i])
print(M_G8$N2[i])
}

## S subs ##

M_G8$S2 <- 0

for (i in 1:nrow(M_G8))
{Nnse$value[i] <- M_G8 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == M_G8$tree.nr[i]) %>% 
  .$S
M_G8$S2[i] <- abs(M_G8$S[i] - Nnse$value[i])
print(M_G8$S2[i])
}


write.table(M_G8,file="M/M_G8_distances_2.txt")




######## G9 ###########


## Data import ##
M_G9 <- read.table("M/M_G9_distances.txt", header = T)
## Build redudant column ##
M_G9 <- M_G9[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(M_G9$tree.nr)
Nnse$value <- 0
M_G9$N2 <- 0

for (i in 1:nrow(M_G9))
{Nnse$value[i] <- M_G9 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == M_G9$tree.nr[i]) %>% 
  .$N
M_G9$N2[i] <- abs(M_G9$N[i] - Nnse$value[i])
print(M_G9$N2[i])
}

## S subs ##

M_G9$S2 <- 0

for (i in 1:nrow(M_G9))
{Nnse$value[i] <- M_G9 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == M_G9$tree.nr[i]) %>% 
  .$S
M_G9$S2[i] <- abs(M_G9$S[i] - Nnse$value[i])
print(M_G9$S2[i])
}


write.table(M_G9,file="M/M_G9_distances_2.txt")




######## G10 ###########


## Data import ##
M_G10 <- read.table("M/M_G10_distances.txt", header = T)
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
M_G10$N2[i] <- abs(M_G10$N[i] - Nnse$value[i])
print(M_G10$N2[i])
}

## S subs ##

M_G10$S2 <- 0

for (i in 1:nrow(M_G10))
{Nnse$value[i] <- M_G10 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == M_G10$tree.nr[i]) %>% 
  .$S
M_G10$S2[i] <- abs(M_G10$S[i] - Nnse$value[i])
print(M_G10$S2[i])
}


write.table(M_G10,file="M/M_G10_distances_2.txt")




######## G11 ###########


## Data import ##
M_G11 <- read.table("M/M_G11_distances.txt", header = T)
## Build redudant column ##
M_G11 <- M_G11[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(M_G11$tree.nr)
Nnse$value <- 0
M_G11$N2 <- 0

for (i in 1:nrow(M_G11))
{Nnse$value[i] <- M_G11 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == M_G11$tree.nr[i]) %>% 
  .$N
M_G11$N2[i] <- abs(M_G11$N[i] - Nnse$value[i])
print(M_G11$N2[i])
}

## S subs ##

M_G11$S2 <- 0

for (i in 1:nrow(M_G11))
{Nnse$value[i] <- M_G11 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == M_G11$tree.nr[i]) %>% 
  .$S
M_G11$S2[i] <- abs(M_G11$S[i] - Nnse$value[i])
print(M_G11$S2[i])
}


write.table(M_G11,file="M/M_G11_distances_2.txt")




######## G12 ###########


## Data import ##
M_G12 <- read.table("M/M_G12_distances.txt", header = T)
## Build redudant column ##
M_G12 <- M_G12[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(M_G12$tree.nr)
Nnse$value <- 0
M_G12$N2 <- 0

for (i in 1:nrow(M_G12))
{Nnse$value[i] <- M_G12 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == M_G12$tree.nr[i]) %>% 
  .$N
M_G12$N2[i] <- abs(M_G12$N[i] - Nnse$value[i])
print(M_G12$N2[i])
}

## S subs ##

M_G12$S2 <- 0

for (i in 1:nrow(M_G12))
{Nnse$value[i] <- M_G12 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == M_G12$tree.nr[i]) %>% 
  .$S
M_G12$S2[i] <- abs(M_G12$S[i] - Nnse$value[i])
print(M_G12$S2[i])
}


write.table(M_G12,file="M/M_G12_distances_2.txt")




######## G13 ###########


## Data import ##
M_G13 <- read.table("M/M_G13_distances.txt", header = T)
## Build redudant column ##
M_G13 <- M_G13[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(M_G13$tree.nr)
Nnse$value <- 0
M_G13$N2 <- 0

for (i in 1:nrow(M_G13))
{Nnse$value[i] <- M_G13 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == M_G13$tree.nr[i]) %>% 
  .$N
M_G13$N2[i] <- abs(M_G13$N[i] - Nnse$value[i])
print(M_G13$N2[i])
}

## S subs ##

M_G13$S2 <- 0

for (i in 1:nrow(M_G13))
{Nnse$value[i] <- M_G13 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == M_G13$tree.nr[i]) %>% 
  .$S
M_G13$S2[i] <- abs(M_G13$S[i] - Nnse$value[i])
print(M_G13$S2[i])
}


write.table(M_G13,file="M/M_G13_distances_2.txt")
