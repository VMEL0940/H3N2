### Influ H3N2 GenDist Proj ###
## Package ##
library(tidyverse)

## Set working directory ##
setwd("01_GenDisFlu/Data/Genetic Distances/Verification_702/")

######## G1 ###########

## Data import ##
NS_G1 <- read.table("NS/NS_G1_distances.txt", header = T)
## Build redudant column ##
NS_G1 <- NS_G1[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NS_G1$tree.nr)
Nnse$value <- 0
NS_G1$N2 <- 0

for (i in 1:nrow(NS_G1))
{Nnse$value[i] <- NS_G1 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NS_G1$tree.nr[i]) %>% 
  .$N
NS_G1$N2[i] <- abs(NS_G1$N[i] - Nnse$value[i])
print(NS_G1$N2[i])
}

## S subs ##

NS_G1$S2 <- 0

for (i in 1:nrow(NS_G1))
{Nnse$value[i] <- NS_G1 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NS_G1$tree.nr[i]) %>% 
  .$S
NS_G1$S2[i] <- abs(NS_G1$S[i] - Nnse$value[i])
print(NS_G1$S2[i])
}


write.table(NS_G1,file="NS/NS_G1_distances_2.txt")


######## G2 ###########


## Data import ##
NS_G2 <- read.table("NS/NS_G2_distances.txt", header = T)
## Build redudant column ##
NS_G2 <- NS_G2[,1:5]


## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NS_G2$tree.nr)
Nnse$value <- 0
NS_G2$N2 <- 0

for (i in 1:nrow(NS_G2))
{Nnse$value[i] <- NS_G2 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NS_G2$tree.nr[i]) %>% 
  .$N
NS_G2$N2[i] <- abs(NS_G2$N[i] - Nnse$value[i])
print(NS_G2$N2[i])
}


## S subs ##

NS_G2$S2 <- 0

for (i in 1:nrow(NS_G2))
{Nnse$value[i] <- NS_G2 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NS_G2$tree.nr[i]) %>% 
  .$S
NS_G2$S2[i] <- abs(NS_G2$S[i] - Nnse$value[i])
print(NS_G2$S2[i])
}


write.table(NS_G2,file="NS/NS_G2_distances_2.txt")





######## G3 ###########


## Data import ##
NS_G3 <- read.table("NS/NS_G3_distances.txt", header = T)
## Build redudant column ##
NS_G3 <- NS_G3[,1:5]


## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NS_G3$tree.nr)
Nnse$value <- 0
NS_G3$N2 <- 0

for (i in 1:nrow(NS_G3))
{Nnse$value[i] <- NS_G3 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NS_G3$tree.nr[i]) %>% 
  .$N
NS_G3$N2[i] <- abs(NS_G3$N[i] - Nnse$value[i])
print(NS_G3$N2[i])
}


## S subs ##

NS_G3$S2 <- 0

for (i in 1:nrow(NS_G3))
{Nnse$value[i] <- NS_G3 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NS_G3$tree.nr[i]) %>% 
  .$S
NS_G3$S2[i] <- abs(NS_G3$S[i] - Nnse$value[i])
print(NS_G3$S2[i])
}

write.table(NS_G3,file="NS/NS_G3_distances_2.txt")





######## G4 ###########


## Data import ##
NS_G4 <- read.table("NS/NS_G4_distances.txt", header = T)
## Build redudant column ##
NS_G4 <- NS_G4[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NS_G4$tree.nr)
Nnse$value <- 0
NS_G4$N2 <- 0

for (i in 1:nrow(NS_G4))
{Nnse$value[i] <- NS_G4 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NS_G4$tree.nr[i]) %>% 
  .$N
NS_G4$N2[i] <- abs(NS_G4$N[i] - Nnse$value[i])
print(NS_G4$N2[i])
}

## S subs ##

NS_G4$S2 <- 0

for (i in 1:nrow(NS_G4))
{Nnse$value[i] <- NS_G4 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NS_G4$tree.nr[i]) %>% 
  .$S
NS_G4$S2[i] <- abs(NS_G4$S[i] - Nnse$value[i])
print(NS_G4$S2[i])
}


write.table(NS_G4,file="NS/NS_G4_distances_2.txt")



######## G5 ###########


## Data import ##
NS_G5 <- read.table("NS/NS_G5_distances.txt", header = T)
## Build redudant column ##
NS_G5 <- NS_G5[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NS_G5$tree.nr)
Nnse$value <- 0
NS_G5$N2 <- 0

for (i in 1:nrow(NS_G5))
{Nnse$value[i] <- NS_G5 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NS_G5$tree.nr[i]) %>% 
  .$N
NS_G5$N2[i] <- abs(NS_G5$N[i] - Nnse$value[i])
print(NS_G5$N2[i])
}

## S subs ##

NS_G5$S2 <- 0

for (i in 1:nrow(NS_G5))
{Nnse$value[i] <- NS_G5 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NS_G5$tree.nr[i]) %>% 
  .$S
NS_G5$S2[i] <- abs(NS_G5$S[i] - Nnse$value[i])
print(NS_G5$S2[i])
}


write.table(NS_G5,file="NS/NS_G5_distances_2.txt")



######## G6 ###########


## Data import ##
NS_G6 <- read.table("NS/NS_G6_distances.txt", header = T)
## Build redudant column ##
NS_G6 <- NS_G6[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NS_G6$tree.nr)
Nnse$value <- 0
NS_G6$N2 <- 0

for (i in 1:nrow(NS_G6))
{Nnse$value[i] <- NS_G6 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NS_G6$tree.nr[i]) %>% 
  .$N
NS_G6$N2[i] <- abs(NS_G6$N[i] - Nnse$value[i])
print(NS_G6$N2[i])
}

## S subs ##

NS_G6$S2 <- 0

for (i in 1:nrow(NS_G6))
{Nnse$value[i] <- NS_G6 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NS_G6$tree.nr[i]) %>% 
  .$S
NS_G6$S2[i] <- abs(NS_G6$S[i] - Nnse$value[i])
print(NS_G6$S2[i])
}


write.table(NS_G6,file="NS/NS_G6_distances_2.txt")



######## G7 ###########


## Data import ##
NS_G7 <- read.table("NS/NS_G7_distances.txt", header = T)
## Build redudant column ##
NS_G7 <- NS_G7[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NS_G7$tree.nr)
Nnse$value <- 0
NS_G7$N2 <- 0

for (i in 1:nrow(NS_G7))
{Nnse$value[i] <- NS_G7 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NS_G7$tree.nr[i]) %>% 
  .$N
NS_G7$N2[i] <- abs(NS_G7$N[i] - Nnse$value[i])
print(NS_G7$N2[i])
}

## S subs ##

NS_G7$S2 <- 0

for (i in 1:nrow(NS_G7))
{Nnse$value[i] <- NS_G7 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NS_G7$tree.nr[i]) %>% 
  .$S
NS_G7$S2[i] <- abs(NS_G7$S[i] - Nnse$value[i])
print(NS_G7$S2[i])
}


write.table(NS_G7,file="NS/NS_G7_distances_2.txt")




######## G8 ###########


## Data import ##
NS_G8 <- read.table("NS/NS_G8_distances.txt", header = T)
## Build redudant column ##
NS_G8 <- NS_G8[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NS_G8$tree.nr)
Nnse$value <- 0
NS_G8$N2 <- 0

for (i in 1:nrow(NS_G8))
{Nnse$value[i] <- NS_G8 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NS_G8$tree.nr[i]) %>% 
  .$N
NS_G8$N2[i] <- abs(NS_G8$N[i] - Nnse$value[i])
print(NS_G8$N2[i])
}

## S subs ##

NS_G8$S2 <- 0

for (i in 1:nrow(NS_G8))
{Nnse$value[i] <- NS_G8 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NS_G8$tree.nr[i]) %>% 
  .$S
NS_G8$S2[i] <- abs(NS_G8$S[i] - Nnse$value[i])
print(NS_G8$S2[i])
}


write.table(NS_G8,file="NS/NS_G8_distances_2.txt")




######## G9 ###########


## Data import ##
NS_G9 <- read.table("NS/NS_G9_distances.txt", header = T)
## Build redudant column ##
NS_G9 <- NS_G9[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NS_G9$tree.nr)
Nnse$value <- 0
NS_G9$N2 <- 0

for (i in 1:nrow(NS_G9))
{Nnse$value[i] <- NS_G9 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NS_G9$tree.nr[i]) %>% 
  .$N
NS_G9$N2[i] <- abs(NS_G9$N[i] - Nnse$value[i])
print(NS_G9$N2[i])
}

## S subs ##

NS_G9$S2 <- 0

for (i in 1:nrow(NS_G9))
{Nnse$value[i] <- NS_G9 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NS_G9$tree.nr[i]) %>% 
  .$S
NS_G9$S2[i] <- abs(NS_G9$S[i] - Nnse$value[i])
print(NS_G9$S2[i])
}


write.table(NS_G9,file="NS/NS_G9_distances_2.txt")




######## G10 ###########


## Data import ##
NS_G10 <- read.table("NS/NS_G10_distances.txt", header = T)
## Build redudant column ##
NS_G10 <- NS_G10[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NS_G10$tree.nr)
Nnse$value <- 0
NS_G10$N2 <- 0

for (i in 1:nrow(NS_G10))
{Nnse$value[i] <- NS_G10 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NS_G10$tree.nr[i]) %>% 
  .$N
NS_G10$N2[i] <- abs(NS_G10$N[i] - Nnse$value[i])
print(NS_G10$N2[i])
}

## S subs ##

NS_G10$S2 <- 0

for (i in 1:nrow(NS_G10))
{Nnse$value[i] <- NS_G10 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NS_G10$tree.nr[i]) %>% 
  .$S
NS_G10$S2[i] <- abs(NS_G10$S[i] - Nnse$value[i])
print(NS_G10$S2[i])
}


write.table(NS_G10,file="NS/NS_G10_distances_2.txt")




######## G11 ###########


## Data import ##
NS_G11 <- read.table("NS/NS_G11_distances.txt", header = T)
## Build redudant column ##
NS_G11 <- NS_G11[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NS_G11$tree.nr)
Nnse$value <- 0
NS_G11$N2 <- 0

for (i in 1:nrow(NS_G11))
{Nnse$value[i] <- NS_G11 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NS_G11$tree.nr[i]) %>% 
  .$N
NS_G11$N2[i] <- abs(NS_G11$N[i] - Nnse$value[i])
print(NS_G11$N2[i])
}

## S subs ##

NS_G11$S2 <- 0

for (i in 1:nrow(NS_G11))
{Nnse$value[i] <- NS_G11 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NS_G11$tree.nr[i]) %>% 
  .$S
NS_G11$S2[i] <- abs(NS_G11$S[i] - Nnse$value[i])
print(NS_G11$S2[i])
}


write.table(NS_G11,file="NS/NS_G11_distances_2.txt")




######## G12 ###########


## Data import ##
NS_G12 <- read.table("NS/NS_G12_distances.txt", header = T)
## Build redudant column ##
NS_G12 <- NS_G12[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NS_G12$tree.nr)
Nnse$value <- 0
NS_G12$N2 <- 0

for (i in 1:nrow(NS_G12))
{Nnse$value[i] <- NS_G12 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NS_G12$tree.nr[i]) %>% 
  .$N
NS_G12$N2[i] <- abs(NS_G12$N[i] - Nnse$value[i])
print(NS_G12$N2[i])
}

## S subs ##

NS_G12$S2 <- 0

for (i in 1:nrow(NS_G12))
{Nnse$value[i] <- NS_G12 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NS_G12$tree.nr[i]) %>% 
  .$S
NS_G12$S2[i] <- abs(NS_G12$S[i] - Nnse$value[i])
print(NS_G12$S2[i])
}


write.table(NS_G12,file="NS/NS_G12_distances_2.txt")




######## G13 ###########


## Data import ##
NS_G13 <- read.table("NS/NS_G13_distances.txt", header = T)
## Build redudant column ##
NS_G13 <- NS_G13[,1:5]

## Extract the "Noise AA sub values ##

## N subs ##
Nnse <- as.data.frame(NS_G13$tree.nr)
Nnse$value <- 0
NS_G13$N2 <- 0

for (i in 1:nrow(NS_G13))
{Nnse$value[i] <- NS_G13 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NS_G13$tree.nr[i]) %>% 
  .$N
NS_G13$N2[i] <- abs(NS_G13$N[i] - Nnse$value[i])
print(NS_G13$N2[i])
}

## S subs ##

NS_G13$S2 <- 0

for (i in 1:nrow(NS_G13))
{Nnse$value[i] <- NS_G13 %>% 
  filter(vaccineStrain == compareStrain) %>% 
  filter(tree.nr == NS_G13$tree.nr[i]) %>% 
  .$S
NS_G13$S2[i] <- abs(NS_G13$S[i] - Nnse$value[i])
print(NS_G13$S2[i])
}


write.table(NS_G13,file="NS/NS_G13_distances_2.txt")
