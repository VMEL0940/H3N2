### Influ H3N2 GenDist Proj ###
## Package ##
library(tidyverse)

## Set working directory ##
setwd("/01_GenDisFlu/")

########### PB1 !!!!!! #################

## Data import ##
PB1_G1 <- read.table("Data/Genetic Distances/PB1/PB1_G1_distances.txt", header = T)
PB1_G1_old <- read.table("Data/Old Distance/PB1 distances/PB1_G1_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
PB1_G1 <- PB1_G1[,1:5]

## Set the error value per each tree
Error <- PB1_G1_old %>% 
  filter(vaccineStrain == unique(PB1_G1$vaccineStrain) & 
           compareStrain == unique(PB1_G1$vaccineStrain))

PB1_G1_V2 <- PB1_G1 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

PB1_G1$N2 <- abs(PB1_G1_V2$N.x - PB1_G1_V2$N.y)
PB1_G1$S2 <- abs(PB1_G1_V2$S.x - PB1_G1_V2$S.y)

write.table(PB1_G1,file= "Data/Genetic Distances/PB1/PB1_G1_distances_2.txt")




## Data import ##
PB1_G2 <- read.table("Data/Genetic Distances/PB1/PB1_G2_distances.txt", header = T)
PB1_G2_old <- read.table("Data/Old Distance/PB1 distances/PB1_G2_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
PB1_G2 <- PB1_G2[,1:5]

## Set the error value per each tree
Error <- PB1_G2_old %>% 
  filter(vaccineStrain == unique(PB1_G2$vaccineStrain) & 
           compareStrain == unique(PB1_G2$vaccineStrain))

PB1_G2_V2 <- PB1_G2 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

PB1_G2$N2 <- abs(PB1_G2_V2$N.x - PB1_G2_V2$N.y)
PB1_G2$S2 <- abs(PB1_G2_V2$S.x - PB1_G2_V2$S.y)

write.table(PB1_G2,file= "Data/Genetic Distances/PB1/PB1_G2_distances_2.txt")





## Data import ##
PB1_G3 <- read.table("Data/Genetic Distances/PB1/PB1_G3_distances.txt", header = T)
PB1_G3_old <- read.table("Data/Old Distance/PB1 distances/PB1_G3_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
PB1_G3 <- PB1_G3[,1:5]



## Set the error value per each tree
Error <- PB1_G3_old %>% 
  filter(vaccineStrain == unique(PB1_G3$vaccineStrain) & 
           compareStrain == unique(PB1_G3$vaccineStrain))

PB1_G3_V2 <- PB1_G3 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

PB1_G3$N2 <- abs(PB1_G3_V2$N.x - PB1_G3_V2$N.y)
PB1_G3$S2 <- abs(PB1_G3_V2$S.x - PB1_G3_V2$S.y)

write.table(PB1_G3,file= "Data/Genetic Distances/PB1/PB1_G3_distances_2.txt")




## Data import ##
PB1_G4 <- read.table("Data/Genetic Distances/PB1/PB1_G4_distances.txt", header = T)
PB1_G4_old <- read.table("Data/Old Distance/PB1 distances/PB1_G4_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
PB1_G4 <- PB1_G4[,1:5]

## Set the error value per each tree
Error <- PB1_G4_old %>% 
  filter(vaccineStrain == unique(PB1_G4$vaccineStrain) & 
           compareStrain == unique(PB1_G4$vaccineStrain))

PB1_G4_V2 <- PB1_G4 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

PB1_G4$N2 <- abs(PB1_G4_V2$N.x - PB1_G4_V2$N.y)
PB1_G4$S2 <- abs(PB1_G4_V2$S.x - PB1_G4_V2$S.y)

write.table(PB1_G4,file= "Data/Genetic Distances/PB1/PB1_G4_distances_2.txt")



## Data import ##
PB1_G5 <- read.table("Data/Genetic Distances/PB1/PB1_G5_distances.txt", header = T)
PB1_G5_old <- read.table("Data/Old Distance/PB1 distances/PB1_G5_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
PB1_G5 <- PB1_G5[,1:5]

## Set the error value per each tree
Error <- PB1_G5_old %>% 
  filter(vaccineStrain == unique(PB1_G5$vaccineStrain) & 
           compareStrain == unique(PB1_G5$vaccineStrain))

PB1_G5_V2 <- PB1_G5 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

PB1_G5$N2 <- abs(PB1_G5_V2$N.x - PB1_G5_V2$N.y)
PB1_G5$S2 <- abs(PB1_G5_V2$S.x - PB1_G5_V2$S.y)

write.table(PB1_G5,file= "Data/Genetic Distances/PB1/PB1_G5_distances_2.txt")





## Data import ##
PB1_G6 <- read.table("Data/Genetic Distances/PB1/PB1_G6_distances.txt", header = T)
PB1_G6_old <- read.table("Data/Old Distance/PB1 distances/PB1_G6_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
PB1_G6 <- PB1_G6[,1:5]

## Set the error value per each tree
Error <- PB1_G6_old %>% 
  filter(vaccineStrain == unique(PB1_G6$vaccineStrain) & 
           compareStrain == unique(PB1_G6$vaccineStrain))

PB1_G6_V2 <- PB1_G6 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

PB1_G6$N2 <- abs(PB1_G6_V2$N.x - PB1_G6_V2$N.y)
PB1_G6$S2 <- abs(PB1_G6_V2$S.x - PB1_G6_V2$S.y)

write.table(PB1_G6,file= "Data/Genetic Distances/PB1/PB1_G6_distances_2.txt")





## Data import ##
PB1_G7 <- read.table("Data/Genetic Distances/PB1/PB1_G7_distances.txt", header = T)
PB1_G7_old <- read.table("Data/Old Distance/PB1 distances/PB1_G7_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
PB1_G7 <- PB1_G7[,1:5]

## Set the error value per each tree
Error <- PB1_G7_old %>% 
  filter(vaccineStrain == unique(PB1_G7$vaccineStrain) & 
           compareStrain == unique(PB1_G7$vaccineStrain))

PB1_G7_V2 <- PB1_G7 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

PB1_G7$N2 <- abs(PB1_G7_V2$N.x - PB1_G7_V2$N.y)
PB1_G7$S2 <- abs(PB1_G7_V2$S.x - PB1_G7_V2$S.y)

write.table(PB1_G7,file= "Data/Genetic Distances/PB1/PB1_G7_distances_2.txt")




## Data import ##
PB1_G8 <- read.table("Data/Genetic Distances/PB1/PB1_G8_distances.txt", header = T)
PB1_G8_old <- read.table("Data/Old Distance/PB1 distances/PB1_G8_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
PB1_G8 <- PB1_G8[,1:5]

## Set the error value per each tree
Error <- PB1_G8_old %>% 
  filter(vaccineStrain == unique(PB1_G8$vaccineStrain) & 
           compareStrain == unique(PB1_G8$vaccineStrain))

PB1_G8_V2 <- PB1_G8 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

PB1_G8$N2 <- abs(PB1_G8_V2$N.x - PB1_G8_V2$N.y)
PB1_G8$S2 <- abs(PB1_G8_V2$S.x - PB1_G8_V2$S.y)

write.table(PB1_G8,file= "Data/Genetic Distances/PB1/PB1_G8_distances_2.txt")





## Data import ##
PB1_G9 <- read.table("Data/Genetic Distances/PB1/PB1_G9_distances.txt", header = T)
PB1_G9_old <- read.table("Data/Old Distance/PB1 distances/PB1_G9_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
PB1_G9 <- PB1_G9[,1:5]

## Set the error value per each tree
Error <- PB1_G9_old %>% 
  filter(vaccineStrain == unique(PB1_G9$vaccineStrain) & 
           compareStrain == unique(PB1_G9$vaccineStrain))

PB1_G9_V2 <- PB1_G9 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

PB1_G9$N2 <- abs(PB1_G9_V2$N.x - PB1_G9_V2$N.y)
PB1_G9$S2 <- abs(PB1_G9_V2$S.x - PB1_G9_V2$S.y)

write.table(PB1_G9,file= "Data/Genetic Distances/PB1/PB1_G9_distances_2.txt")


