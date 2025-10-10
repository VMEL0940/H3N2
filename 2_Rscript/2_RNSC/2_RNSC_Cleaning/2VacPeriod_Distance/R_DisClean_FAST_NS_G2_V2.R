### Influ H3N2 GenDist Proj ###
## Package ##
library(tidyverse)

## Set working directory ##
setwd("/01_GenDisFlu/")

########### NS !!!!!! #################

## Data import ##
NS_G1 <- read.table("Data/Genetic Distances/G2/NS/NS_G1_distances.txt", header = T)
NS_G1_old <- read.table("Data/Genetic Distances/0_Old Distance/V2/NS distances/NS_G1_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
NS_G1 <- NS_G1[,1:5]

## Set the error value per each tree
Error <- NS_G1_old %>% 
  filter(vaccineStrain == unique(NS_G1_old$vaccineStrain) & 
           compareStrain == unique(NS_G1_old$vaccineStrain))

NS_G1_V2 <- NS_G1 %>% 
  full_join(Error, by = c("tree.nr"))

NS_G1$N2 <- abs(NS_G1_V2$N.x - NS_G1_V2$N.y)
NS_G1$S2 <- abs(NS_G1_V2$S.x - NS_G1_V2$S.y)

write.table(NS_G1,file= "Data/Genetic Distances/G2/NS/NS_G1_distances_2.txt")




## Data import ##
NS_G2 <- read.table("Data/Genetic Distances/G2/NS/NS_G2_distances.txt", header = T)
NS_G2_old <- read.table("Data/Genetic Distances/0_Old Distance/V2/NS distances/NS_G2_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
NS_G2 <- NS_G2[,1:5]

## Set the error value per each tree
Error <- NS_G2_old %>% 
  filter(vaccineStrain == unique(NS_G2_old$vaccineStrain) & 
           compareStrain == unique(NS_G2_old$vaccineStrain))

NS_G2_V2 <- NS_G2 %>% 
  full_join(Error, by = c("tree.nr"))

NS_G2$N2 <- abs(NS_G2_V2$N.x - NS_G2_V2$N.y)
NS_G2$S2 <- abs(NS_G2_V2$S.x - NS_G2_V2$S.y)

write.table(NS_G2,file= "Data/Genetic Distances/G2/NS/NS_G2_distances_2.txt")





## Data import ##
NS_G3 <- read.table("Data/Genetic Distances/G2/NS/NS_G3_distances.txt", header = T)
NS_G3_old <- read.table("Data/Genetic Distances/0_Old Distance/V2/NS distances/NS_G3_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
NS_G3 <- NS_G3[,1:5]



## Set the error value per each tree
Error <- NS_G3_old %>% 
  filter(vaccineStrain == unique(NS_G3_old$vaccineStrain) & 
           compareStrain == unique(NS_G3_old$vaccineStrain))

NS_G3_V2 <- NS_G3 %>% 
  full_join(Error, by = c("tree.nr"))

NS_G3$N2 <- abs(NS_G3_V2$N.x - NS_G3_V2$N.y)
NS_G3$S2 <- abs(NS_G3_V2$S.x - NS_G3_V2$S.y)

write.table(NS_G3,file= "Data/Genetic Distances/G2/NS/NS_G3_distances_2.txt")




## Data import ##
NS_G4 <- read.table("Data/Genetic Distances/G2/NS/NS_G4_distances.txt", header = T)
NS_G4_old <- read.table("Data/Genetic Distances/0_Old Distance/V2/NS distances/NS_G4_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
NS_G4 <- NS_G4[,1:5]

## Set the error value per each tree
Error <- NS_G4_old %>% 
  filter(vaccineStrain == unique(NS_G4_old$vaccineStrain) & 
           compareStrain == unique(NS_G4_old$vaccineStrain))

NS_G4_V2 <- NS_G4 %>% 
  full_join(Error, by = c("tree.nr"))

NS_G4$N2 <- abs(NS_G4_V2$N.x - NS_G4_V2$N.y)
NS_G4$S2 <- abs(NS_G4_V2$S.x - NS_G4_V2$S.y)

write.table(NS_G4,file= "Data/Genetic Distances/G2/NS/NS_G4_distances_2.txt")



## Data import ##
NS_G5 <- read.table("Data/Genetic Distances/G2/NS/NS_G5_distances.txt", header = T)
NS_G5_old <- read.table("Data/Genetic Distances/0_Old Distance/V2/NS distances/NS_G5_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
NS_G5 <- NS_G5[,1:5]

## Set the error value per each tree
Error <- NS_G5_old %>% 
  filter(vaccineStrain == unique(NS_G5_old$vaccineStrain) & 
           compareStrain == unique(NS_G5_old$vaccineStrain))

NS_G5_V2 <- NS_G5 %>% 
  full_join(Error, by = c("tree.nr"))

NS_G5$N2 <- abs(NS_G5_V2$N.x - NS_G5_V2$N.y)
NS_G5$S2 <- abs(NS_G5_V2$S.x - NS_G5_V2$S.y)

write.table(NS_G5,file= "Data/Genetic Distances/G2/NS/NS_G5_distances_2.txt")





## Data import ##
NS_G6 <- read.table("Data/Genetic Distances/G2/NS/NS_G6_distances.txt", header = T)
NS_G6_old <- read.table("Data/Genetic Distances/0_Old Distance/V2/NS distances/NS_G6_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
NS_G6 <- NS_G6[,1:5]

## Set the error value per each tree
Error <- NS_G6_old %>% 
  filter(vaccineStrain == unique(NS_G6_old$vaccineStrain) & 
           compareStrain == unique(NS_G6_old$vaccineStrain))

NS_G6_V2 <- NS_G6 %>% 
  full_join(Error, by = c("tree.nr"))

NS_G6$N2 <- abs(NS_G6_V2$N.x - NS_G6_V2$N.y)
NS_G6$S2 <- abs(NS_G6_V2$S.x - NS_G6_V2$S.y)

write.table(NS_G6,file= "Data/Genetic Distances/G2/NS/NS_G6_distances_2.txt")





## Data import ##
NS_G7 <- read.table("Data/Genetic Distances/G2/NS/NS_G7_distances.txt", header = T)
NS_G7_old <- read.table("Data/Genetic Distances/0_Old Distance/V2/NS distances/NS_G7_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
NS_G7 <- NS_G7[,1:5]

## Set the error value per each tree
Error <- NS_G7_old %>% 
  filter(vaccineStrain == unique(NS_G7_old$vaccineStrain) & 
           compareStrain == unique(NS_G7_old$vaccineStrain))

NS_G7_V2 <- NS_G7 %>% 
  full_join(Error, by = c("tree.nr"))

NS_G7$N2 <- abs(NS_G7_V2$N.x - NS_G7_V2$N.y)
NS_G7$S2 <- abs(NS_G7_V2$S.x - NS_G7_V2$S.y)

write.table(NS_G7,file= "Data/Genetic Distances/G2/NS/NS_G7_distances_2.txt")




## Data import ##
NS_G8 <- read.table("Data/Genetic Distances/G2/NS/NS_G8_distances.txt", header = T)
NS_G8_old <- read.table("Data/Genetic Distances/0_Old Distance/V2/NS distances/NS_G8_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
NS_G8 <- NS_G8[,1:5]

## Set the error value per each tree
Error <- NS_G8_old %>% 
  filter(vaccineStrain == unique(NS_G8_old$vaccineStrain) & 
           compareStrain == unique(NS_G8_old$vaccineStrain))

NS_G8_V2 <- NS_G8 %>% 
  full_join(Error, by = c("tree.nr"))

NS_G8$N2 <- abs(NS_G8_V2$N.x - NS_G8_V2$N.y)
NS_G8$S2 <- abs(NS_G8_V2$S.x - NS_G8_V2$S.y)

write.table(NS_G8,file= "Data/Genetic Distances/G2/NS/NS_G8_distances_2.txt")





