### Influ H3N2 GenDist Proj ###
## Package ##
library(tidyverse)

## Set working directory ##
setwd("01_GenDisFlu/Data/Genetic Distances/Verification_702/")

######## G1 ###########

## Data import ##
NA_G1 <- read.table("NA/NA_G1_distances.txt", header = T)
## Build redudant column ##
NA_G1 <- NA_G1[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NA_G1$tree.nr)
Nnse$value <- 0
NA_G1$N2 <- 0

for (i in 1:nrow(NA_G1))
{Nnse$value[i] <- NA_G1 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NA_G1$tree.nr[i]) %>% 
  .$N
NA_G1$N2[i] <- abs(NA_G1$N[i] - Nnse$value[i])
print(NA_G1$N2[i])
}

## S subs ##

NA_G1$S2 <- 0

for (i in 1:nrow(NA_G1))
{Nnse$value[i] <- NA_G1 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NA_G1$tree.nr[i]) %>% 
  .$S
NA_G1$S2[i] <- abs(NA_G1$S[i] - Nnse$value[i])
print(NA_G1$S2[i])
}


write.table(NA_G1,file="NA/NA_G1_distances_2.txt")


######## G2 ###########


## Data import ##
NA_G2 <- read.table("NA/NA_G2_distances.txt", header = T)
## Build redudant column ##
NA_G2 <- NA_G2[,1:5]


## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NA_G2$tree.nr)
Nnse$value <- 0
NA_G2$N2 <- 0

for (i in 1:nrow(NA_G2))
{Nnse$value[i] <- NA_G2 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NA_G2$tree.nr[i]) %>% 
  .$N
NA_G2$N2[i] <- abs(NA_G2$N[i] - Nnse$value[i])
print(NA_G2$N2[i])
}


## S subs ##

NA_G2$S2 <- 0

for (i in 1:nrow(NA_G2))
{Nnse$value[i] <- NA_G2 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NA_G2$tree.nr[i]) %>% 
  .$S
NA_G2$S2[i] <- abs(NA_G2$S[i] - Nnse$value[i])
print(NA_G2$S2[i])
}


write.table(NA_G2,file="NA/NA_G2_distances_2.txt")





######## G3 ###########


## Data import ##
NA_G3 <- read.table("NA/NA_G3_distances.txt", header = T)
## Build redudant column ##
NA_G3 <- NA_G3[,1:5]


## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NA_G3$tree.nr)
Nnse$value <- 0
NA_G3$N2 <- 0

for (i in 1:nrow(NA_G3))
{Nnse$value[i] <- NA_G3 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NA_G3$tree.nr[i]) %>% 
  .$N
NA_G3$N2[i] <- abs(NA_G3$N[i] - Nnse$value[i])
print(NA_G3$N2[i])
}


## S subs ##

NA_G3$S2 <- 0

for (i in 1:nrow(NA_G3))
{Nnse$value[i] <- NA_G3 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NA_G3$tree.nr[i]) %>% 
  .$S
NA_G3$S2[i] <- abs(NA_G3$S[i] - Nnse$value[i])
print(NA_G3$S2[i])
}

write.table(NA_G3,file="NA/NA_G3_distances_2.txt")





######## G4 ###########


## Data import ##
NA_G4 <- read.table("NA/NA_G4_distances.txt", header = T)
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
NA_G4$N2[i] <- abs(NA_G4$N[i] - Nnse$value[i])
print(NA_G4$N2[i])
}

## S subs ##

NA_G4$S2 <- 0

for (i in 1:nrow(NA_G4))
{Nnse$value[i] <- NA_G4 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NA_G4$tree.nr[i]) %>% 
  .$S
NA_G4$S2[i] <- abs(NA_G4$S[i] - Nnse$value[i])
print(NA_G4$S2[i])
}


write.table(NA_G4,file="NA/NA_G4_distances_2.txt")



######## G5 ###########


## Data import ##
NA_G5 <- read.table("NA/NA_G5_distances.txt", header = T)
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
NA_G5$N2[i] <- abs(NA_G5$N[i] - Nnse$value[i])
print(NA_G5$N2[i])
}

## S subs ##

NA_G5$S2 <- 0

for (i in 1:nrow(NA_G5))
{Nnse$value[i] <- NA_G5 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NA_G5$tree.nr[i]) %>% 
  .$S
NA_G5$S2[i] <- abs(NA_G5$S[i] - Nnse$value[i])
print(NA_G5$S2[i])
}


write.table(NA_G5,file="NA/NA_G5_distances_2.txt")



######## G6 ###########


## Data import ##
NA_G6 <- read.table("NA/NA_G6_distances.txt", header = T)
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
NA_G6$N2[i] <- abs(NA_G6$N[i] - Nnse$value[i])
print(NA_G6$N2[i])
}

## S subs ##

NA_G6$S2 <- 0

for (i in 1:nrow(NA_G6))
{Nnse$value[i] <- NA_G6 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NA_G6$tree.nr[i]) %>% 
  .$S
NA_G6$S2[i] <- abs(NA_G6$S[i] - Nnse$value[i])
print(NA_G6$S2[i])
}


write.table(NA_G6,file="NA/NA_G6_distances_2.txt")



######## G7 ###########


## Data import ##
NA_G7 <- read.table("NA/NA_G7_distances.txt", header = T)
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
NA_G7$N2[i] <- abs(NA_G7$N[i] - Nnse$value[i])
print(NA_G7$N2[i])
}

## S subs ##

NA_G7$S2 <- 0

for (i in 1:nrow(NA_G7))
{Nnse$value[i] <- NA_G7 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NA_G7$tree.nr[i]) %>% 
  .$S
NA_G7$S2[i] <- abs(NA_G7$S[i] - Nnse$value[i])
print(NA_G7$S2[i])
}


write.table(NA_G7,file="NA/NA_G7_distances_2.txt")




######## G8 ###########


## Data import ##
NA_G8 <- read.table("NA/NA_G8_distances.txt", header = T)
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
NA_G8$N2[i] <- abs(NA_G8$N[i] - Nnse$value[i])
print(NA_G8$N2[i])
}

## S subs ##

NA_G8$S2 <- 0

for (i in 1:nrow(NA_G8))
{Nnse$value[i] <- NA_G8 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NA_G8$tree.nr[i]) %>% 
  .$S
NA_G8$S2[i] <- abs(NA_G8$S[i] - Nnse$value[i])
print(NA_G8$S2[i])
}


write.table(NA_G8,file="NA/NA_G8_distances_2.txt")




######## G9 ###########


## Data import ##
NA_G9 <- read.table("NA/NA_G9_distances.txt", header = T)
## Build redudant column ##
NA_G9 <- NA_G9[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NA_G9$tree.nr)
Nnse$value <- 0
NA_G9$N2 <- 0

for (i in 1:nrow(NA_G9))
{Nnse$value[i] <- NA_G9 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NA_G9$tree.nr[i]) %>% 
  .$N
NA_G9$N2[i] <- abs(NA_G9$N[i] - Nnse$value[i])
print(NA_G9$N2[i])
}

## S subs ##

NA_G9$S2 <- 0

for (i in 1:nrow(NA_G9))
{Nnse$value[i] <- NA_G9 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NA_G9$tree.nr[i]) %>% 
  .$S
NA_G9$S2[i] <- abs(NA_G9$S[i] - Nnse$value[i])
print(NA_G9$S2[i])
}


write.table(NA_G9,file="NA/NA_G9_distances_2.txt")




######## G10 ###########


## Data import ##
NA_G10 <- read.table("NA/NA_G10_distances.txt", header = T)
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
NA_G10$N2[i] <- abs(NA_G10$N[i] - Nnse$value[i])
print(NA_G10$N2[i])
}

## S subs ##

NA_G10$S2 <- 0

for (i in 1:nrow(NA_G10))
{Nnse$value[i] <- NA_G10 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NA_G10$tree.nr[i]) %>% 
  .$S
NA_G10$S2[i] <- abs(NA_G10$S[i] - Nnse$value[i])
print(NA_G10$S2[i])
}


write.table(NA_G10,file="NA/NA_G10_distances_2.txt")




######## G11 ###########


## Data import ##
NA_G11 <- read.table("NA/NA_G11_distances.txt", header = T)
## Build redudant column ##
NA_G11 <- NA_G11[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NA_G11$tree.nr)
Nnse$value <- 0
NA_G11$N2 <- 0

for (i in 1:nrow(NA_G11))
{Nnse$value[i] <- NA_G11 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NA_G11$tree.nr[i]) %>% 
  .$N
NA_G11$N2[i] <- abs(NA_G11$N[i] - Nnse$value[i])
print(NA_G11$N2[i])
}

## S subs ##

NA_G11$S2 <- 0

for (i in 1:nrow(NA_G11))
{Nnse$value[i] <- NA_G11 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NA_G11$tree.nr[i]) %>% 
  .$S
NA_G11$S2[i] <- abs(NA_G11$S[i] - Nnse$value[i])
print(NA_G11$S2[i])
}


write.table(NA_G11,file="NA/NA_G11_distances_2.txt")




######## G12 ###########


## Data import ##
NA_G12 <- read.table("NA/NA_G12_distances.txt", header = T)
## Build redudant column ##
NA_G12 <- NA_G12[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NA_G12$tree.nr)
Nnse$value <- 0
NA_G12$N2 <- 0

for (i in 1:nrow(NA_G12))
{Nnse$value[i] <- NA_G12 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NA_G12$tree.nr[i]) %>% 
  .$N
NA_G12$N2[i] <- abs(NA_G12$N[i] - Nnse$value[i])
print(NA_G12$N2[i])
}

## S subs ##

NA_G12$S2 <- 0

for (i in 1:nrow(NA_G12))
{Nnse$value[i] <- NA_G12 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NA_G12$tree.nr[i]) %>% 
  .$S
NA_G12$S2[i] <- abs(NA_G12$S[i] - Nnse$value[i])
print(NA_G12$S2[i])
}


write.table(NA_G12,file="NA/NA_G12_distances_2.txt")




######## G13 ###########


## Data import ##
NA_G13 <- read.table("NA/NA_G13_distances.txt", header = T)
## Build redudant column ##
NA_G13 <- NA_G13[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NA_G13$tree.nr)
Nnse$value <- 0
NA_G13$N2 <- 0

for (i in 1:nrow(NA_G13))
{Nnse$value[i] <- NA_G13 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NA_G13$tree.nr[i]) %>% 
  .$N
NA_G13$N2[i] <- abs(NA_G13$N[i] - Nnse$value[i])
print(NA_G13$N2[i])
}

## S subs ##

NA_G13$S2 <- 0

for (i in 1:nrow(NA_G13))
{Nnse$value[i] <- NA_G13 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NA_G13$tree.nr[i]) %>% 
  .$S
NA_G13$S2[i] <- abs(NA_G13$S[i] - Nnse$value[i])
print(NA_G13$S2[i])
}


write.table(NA_G13,file="NA/NA_G13_distances_2.txt")
