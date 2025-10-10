### Influ H3N2 GenDist Proj ###
## Package ##
library(tidyverse)

## Set working directory ##
setwd("01_GenDisFlu/Data/Genetic Distances/Verification_702/")

######## G1 ###########

## Data import ##
PA_G1 <- read.table("PA/PA_G1_distances.txt", header = T)
## Build redudant column ##
PA_G1 <- PA_G1[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(PA_G1$tree.nr)
Nnse$value <- 0
PA_G1$N2 <- 0

for (i in 1:nrow(PA_G1))
{Nnse$value[i] <- PA_G1 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PA_G1$tree.nr[i]) %>% 
  .$N
PA_G1$N2[i] <- abs(PA_G1$N[i] - Nnse$value[i])
print(PA_G1$N2[i])
}

## S subs ##

PA_G1$S2 <- 0

for (i in 1:nrow(PA_G1))
{Nnse$value[i] <- PA_G1 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PA_G1$tree.nr[i]) %>% 
  .$S
PA_G1$S2[i] <- abs(PA_G1$S[i] - Nnse$value[i])
print(PA_G1$S2[i])
}


write.table(PA_G1,file="PA/PA_G1_distances_2.txt")


######## G2 ###########


## Data import ##
PA_G2 <- read.table("PA/PA_G2_distances.txt", header = T)
## Build redudant column ##
PA_G2 <- PA_G2[,1:5]


## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(PA_G2$tree.nr)
Nnse$value <- 0
PA_G2$N2 <- 0

for (i in 1:nrow(PA_G2))
{Nnse$value[i] <- PA_G2 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PA_G2$tree.nr[i]) %>% 
  .$N
PA_G2$N2[i] <- abs(PA_G2$N[i] - Nnse$value[i])
print(PA_G2$N2[i])
}


## S subs ##

PA_G2$S2 <- 0

for (i in 1:nrow(PA_G2))
{Nnse$value[i] <- PA_G2 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PA_G2$tree.nr[i]) %>% 
  .$S
PA_G2$S2[i] <- abs(PA_G2$S[i] - Nnse$value[i])
print(PA_G2$S2[i])
}


write.table(PA_G2,file="PA/PA_G2_distances_2.txt")





######## G3 ###########


## Data import ##
PA_G3 <- read.table("PA/PA_G3_distances.txt", header = T)
## Build redudant column ##
PA_G3 <- PA_G3[,1:5]


## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(PA_G3$tree.nr)
Nnse$value <- 0
PA_G3$N2 <- 0

for (i in 1:nrow(PA_G3))
{Nnse$value[i] <- PA_G3 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PA_G3$tree.nr[i]) %>% 
  .$N
PA_G3$N2[i] <- abs(PA_G3$N[i] - Nnse$value[i])
print(PA_G3$N2[i])
}


## S subs ##

PA_G3$S2 <- 0

for (i in 1:nrow(PA_G3))
{Nnse$value[i] <- PA_G3 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PA_G3$tree.nr[i]) %>% 
  .$S
PA_G3$S2[i] <- abs(PA_G3$S[i] - Nnse$value[i])
print(PA_G3$S2[i])
}

write.table(PA_G3,file="PA/PA_G3_distances_2.txt")





######## G4 ###########


## Data import ##
PA_G4 <- read.table("PA/PA_G4_distances.txt", header = T)
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
PA_G4$N2[i] <- abs(PA_G4$N[i] - Nnse$value[i])
print(PA_G4$N2[i])
}

## S subs ##

PA_G4$S2 <- 0

for (i in 1:nrow(PA_G4))
{Nnse$value[i] <- PA_G4 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PA_G4$tree.nr[i]) %>% 
  .$S
PA_G4$S2[i] <- abs(PA_G4$S[i] - Nnse$value[i])
print(PA_G4$S2[i])
}


write.table(PA_G4,file="PA/PA_G4_distances_2.txt")



######## G5 ###########


## Data import ##
PA_G5 <- read.table("PA/PA_G5_distances.txt", header = T)
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
PA_G5$N2[i] <- abs(PA_G5$N[i] - Nnse$value[i])
print(PA_G5$N2[i])
}

## S subs ##

PA_G5$S2 <- 0

for (i in 1:nrow(PA_G5))
{Nnse$value[i] <- PA_G5 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PA_G5$tree.nr[i]) %>% 
  .$S
PA_G5$S2[i] <- abs(PA_G5$S[i] - Nnse$value[i])
print(PA_G5$S2[i])
}


write.table(PA_G5,file="PA/PA_G5_distances_2.txt")



######## G6 ###########


## Data import ##
PA_G6 <- read.table("PA/PA_G6_distances.txt", header = T)
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
PA_G6$N2[i] <- abs(PA_G6$N[i] - Nnse$value[i])
print(PA_G6$N2[i])
}

## S subs ##

PA_G6$S2 <- 0

for (i in 1:nrow(PA_G6))
{Nnse$value[i] <- PA_G6 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PA_G6$tree.nr[i]) %>% 
  .$S
PA_G6$S2[i] <- abs(PA_G6$S[i] - Nnse$value[i])
print(PA_G6$S2[i])
}


write.table(PA_G6,file="PA/PA_G6_distances_2.txt")



######## G7 ###########


## Data import ##
PA_G7 <- read.table("PA/PA_G7_distances.txt", header = T)
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
PA_G7$N2[i] <- abs(PA_G7$N[i] - Nnse$value[i])
print(PA_G7$N2[i])
}

## S subs ##

PA_G7$S2 <- 0

for (i in 1:nrow(PA_G7))
{Nnse$value[i] <- PA_G7 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PA_G7$tree.nr[i]) %>% 
  .$S
PA_G7$S2[i] <- abs(PA_G7$S[i] - Nnse$value[i])
print(PA_G7$S2[i])
}


write.table(PA_G7,file="PA/PA_G7_distances_2.txt")




######## G8 ###########


## Data import ##
PA_G8 <- read.table("PA/PA_G8_distances.txt", header = T)
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
PA_G8$N2[i] <- abs(PA_G8$N[i] - Nnse$value[i])
print(PA_G8$N2[i])
}

## S subs ##

PA_G8$S2 <- 0

for (i in 1:nrow(PA_G8))
{Nnse$value[i] <- PA_G8 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PA_G8$tree.nr[i]) %>% 
  .$S
PA_G8$S2[i] <- abs(PA_G8$S[i] - Nnse$value[i])
print(PA_G8$S2[i])
}


write.table(PA_G8,file="PA/PA_G8_distances_2.txt")




######## G9 ###########


## Data import ##
PA_G9 <- read.table("PA/PA_G9_distances.txt", header = T)
## Build redudant column ##
PA_G9 <- PA_G9[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(PA_G9$tree.nr)
Nnse$value <- 0
PA_G9$N2 <- 0

for (i in 1:nrow(PA_G9))
{Nnse$value[i] <- PA_G9 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PA_G9$tree.nr[i]) %>% 
  .$N
PA_G9$N2[i] <- abs(PA_G9$N[i] - Nnse$value[i])
print(PA_G9$N2[i])
}

## S subs ##

PA_G9$S2 <- 0

for (i in 1:nrow(PA_G9))
{Nnse$value[i] <- PA_G9 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PA_G9$tree.nr[i]) %>% 
  .$S
PA_G9$S2[i] <- abs(PA_G9$S[i] - Nnse$value[i])
print(PA_G9$S2[i])
}


write.table(PA_G9,file="PA/PA_G9_distances_2.txt")




######## G10 ###########


## Data import ##
PA_G10 <- read.table("PA/PA_G10_distances.txt", header = T)
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
PA_G10$N2[i] <- abs(PA_G10$N[i] - Nnse$value[i])
print(PA_G10$N2[i])
}

## S subs ##

PA_G10$S2 <- 0

for (i in 1:nrow(PA_G10))
{Nnse$value[i] <- PA_G10 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PA_G10$tree.nr[i]) %>% 
  .$S
PA_G10$S2[i] <- abs(PA_G10$S[i] - Nnse$value[i])
print(PA_G10$S2[i])
}


write.table(PA_G10,file="PA/PA_G10_distances_2.txt")




######## G11 ###########


## Data import ##
PA_G11 <- read.table("PA/PA_G11_distances.txt", header = T)
## Build redudant column ##
PA_G11 <- PA_G11[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(PA_G11$tree.nr)
Nnse$value <- 0
PA_G11$N2 <- 0

for (i in 1:nrow(PA_G11))
{Nnse$value[i] <- PA_G11 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PA_G11$tree.nr[i]) %>% 
  .$N
PA_G11$N2[i] <- abs(PA_G11$N[i] - Nnse$value[i])
print(PA_G11$N2[i])
}

## S subs ##

PA_G11$S2 <- 0

for (i in 1:nrow(PA_G11))
{Nnse$value[i] <- PA_G11 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PA_G11$tree.nr[i]) %>% 
  .$S
PA_G11$S2[i] <- abs(PA_G11$S[i] - Nnse$value[i])
print(PA_G11$S2[i])
}


write.table(PA_G11,file="PA/PA_G11_distances_2.txt")




######## G12 ###########


## Data import ##
PA_G12 <- read.table("PA/PA_G12_distances.txt", header = T)
## Build redudant column ##
PA_G12 <- PA_G12[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(PA_G12$tree.nr)
Nnse$value <- 0
PA_G12$N2 <- 0

for (i in 1:nrow(PA_G12))
{Nnse$value[i] <- PA_G12 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PA_G12$tree.nr[i]) %>% 
  .$N
PA_G12$N2[i] <- abs(PA_G12$N[i] - Nnse$value[i])
print(PA_G12$N2[i])
}

## S subs ##

PA_G12$S2 <- 0

for (i in 1:nrow(PA_G12))
{Nnse$value[i] <- PA_G12 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PA_G12$tree.nr[i]) %>% 
  .$S
PA_G12$S2[i] <- abs(PA_G12$S[i] - Nnse$value[i])
print(PA_G12$S2[i])
}


write.table(PA_G12,file="PA/PA_G12_distances_2.txt")




######## G13 ###########


## Data import ##
PA_G13 <- read.table("PA/PA_G13_distances.txt", header = T)
## Build redudant column ##
PA_G13 <- PA_G13[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(PA_G13$tree.nr)
Nnse$value <- 0
PA_G13$N2 <- 0

for (i in 1:nrow(PA_G13))
{Nnse$value[i] <- PA_G13 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PA_G13$tree.nr[i]) %>% 
  .$N
PA_G13$N2[i] <- abs(PA_G13$N[i] - Nnse$value[i])
print(PA_G13$N2[i])
}

## S subs ##

PA_G13$S2 <- 0

for (i in 1:nrow(PA_G13))
{Nnse$value[i] <- PA_G13 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == PA_G13$tree.nr[i]) %>% 
  .$S
PA_G13$S2[i] <- abs(PA_G13$S[i] - Nnse$value[i])
print(PA_G13$S2[i])
}


write.table(PA_G13,file="PA/PA_G13_distances_2.txt")
