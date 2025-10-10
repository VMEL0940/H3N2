### Influ H3N2 GenDist Proj ###
## Package ##
library(tidyverse)

## Set working directory ##
setwd("/01_GenDisFlu/")

########### M !!!!!! #################

## Data import ##
M_G1 <- read.table("Data/Genetic Distances/G2/M/M_G1_distances.txt", header = T)
M_G1_old <- read.table("Data/Genetic Distances/0_Old Distance/V2/M distances/M_G1_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
M_G1 <- M_G1[,1:5]

## Set the error value per each tree
Error <- M_G1_old %>% 
  filter(vaccineStrain == unique(M_G1_old$vaccineStrain) & 
           compareStrain == unique(M_G1_old$vaccineStrain))

M_G1_V2 <- M_G1 %>% 
  full_join(Error, by = c("tree.nr"))

M_G1$N2 <- abs(M_G1_V2$N.x - M_G1_V2$N.y)
M_G1$S2 <- abs(M_G1_V2$S.x - M_G1_V2$S.y)

write.table(M_G1,file= "Data/Genetic Distances/G2/M/M_G1_distances_2.txt")




## Data import ##
M_G2 <- read.table("Data/Genetic Distances/G2/M/M_G2_distances.txt", header = T)
M_G2_old <- read.table("Data/Genetic Distances/0_Old Distance/V2/M distances/M_G2_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
M_G2 <- M_G2[,1:5]

## Set the error value per each tree
Error <- M_G2_old %>% 
  filter(vaccineStrain == unique(M_G2_old$vaccineStrain) & 
           compareStrain == unique(M_G2_old$vaccineStrain))

M_G2_V2 <- M_G2 %>% 
  full_join(Error, by = c("tree.nr"))

M_G2$N2 <- abs(M_G2_V2$N.x - M_G2_V2$N.y)
M_G2$S2 <- abs(M_G2_V2$S.x - M_G2_V2$S.y)

write.table(M_G2,file= "Data/Genetic Distances/G2/M/M_G2_distances_2.txt")





## Data import ##
M_G3 <- read.table("Data/Genetic Distances/G2/M/M_G3_distances.txt", header = T)
M_G3_old <- read.table("Data/Genetic Distances/0_Old Distance/V2/M distances/M_G3_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
M_G3 <- M_G3[,1:5]



## Set the error value per each tree
Error <- M_G3_old %>% 
  filter(vaccineStrain == unique(M_G3_old$vaccineStrain) & 
           compareStrain == unique(M_G3_old$vaccineStrain))

M_G3_V2 <- M_G3 %>% 
  full_join(Error, by = c("tree.nr"))

M_G3$N2 <- abs(M_G3_V2$N.x - M_G3_V2$N.y)
M_G3$S2 <- abs(M_G3_V2$S.x - M_G3_V2$S.y)

write.table(M_G3,file= "Data/Genetic Distances/G2/M/M_G3_distances_2.txt")




## Data import ##
M_G4 <- read.table("Data/Genetic Distances/G2/M/M_G4_distances.txt", header = T)
M_G4_old <- read.table("Data/Genetic Distances/0_Old Distance/V2/M distances/M_G4_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
M_G4 <- M_G4[,1:5]

## Set the error value per each tree
Error <- M_G4_old %>% 
  filter(vaccineStrain == unique(M_G4_old$vaccineStrain) & 
           compareStrain == unique(M_G4_old$vaccineStrain))

M_G4_V2 <- M_G4 %>% 
  full_join(Error, by = c("tree.nr"))

M_G4$N2 <- abs(M_G4_V2$N.x - M_G4_V2$N.y)
M_G4$S2 <- abs(M_G4_V2$S.x - M_G4_V2$S.y)

write.table(M_G4,file= "Data/Genetic Distances/G2/M/M_G4_distances_2.txt")



## Data import ##
M_G5 <- read.table("Data/Genetic Distances/G2/M/M_G5_distances.txt", header = T)
M_G5_old <- read.table("Data/Genetic Distances/0_Old Distance/V2/M distances/M_G5_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
M_G5 <- M_G5[,1:5]

## Set the error value per each tree
Error <- M_G5_old %>% 
  filter(vaccineStrain == unique(M_G5_old$vaccineStrain) & 
           compareStrain == unique(M_G5_old$vaccineStrain))

M_G5_V2 <- M_G5 %>% 
  full_join(Error, by = c("tree.nr"))

M_G5$N2 <- abs(M_G5_V2$N.x - M_G5_V2$N.y)
M_G5$S2 <- abs(M_G5_V2$S.x - M_G5_V2$S.y)

write.table(M_G5,file= "Data/Genetic Distances/G2/M/M_G5_distances_2.txt")





## Data import ##
M_G6 <- read.table("Data/Genetic Distances/G2/M/M_G6_distances.txt", header = T)
M_G6_old <- read.table("Data/Genetic Distances/0_Old Distance/V2/M distances/M_G6_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
M_G6 <- M_G6[,1:5]

## Set the error value per each tree
Error <- M_G6_old %>% 
  filter(vaccineStrain == unique(M_G6_old$vaccineStrain) & 
           compareStrain == unique(M_G6_old$vaccineStrain))

M_G6_V2 <- M_G6 %>% 
  full_join(Error, by = c("tree.nr"))

M_G6$N2 <- abs(M_G6_V2$N.x - M_G6_V2$N.y)
M_G6$S2 <- abs(M_G6_V2$S.x - M_G6_V2$S.y)

write.table(M_G6,file= "Data/Genetic Distances/G2/M/M_G6_distances_2.txt")





## Data import ##
M_G7 <- read.table("Data/Genetic Distances/G2/M/M_G7_distances.txt", header = T)
M_G7_old <- read.table("Data/Genetic Distances/0_Old Distance/V2/M distances/M_G7_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
M_G7 <- M_G7[,1:5]

## Set the error value per each tree
Error <- M_G7_old %>% 
  filter(vaccineStrain == unique(M_G7_old$vaccineStrain) & 
           compareStrain == unique(M_G7_old$vaccineStrain))

M_G7_V2 <- M_G7 %>% 
  full_join(Error, by = c("tree.nr"))

M_G7$N2 <- abs(M_G7_V2$N.x - M_G7_V2$N.y)
M_G7$S2 <- abs(M_G7_V2$S.x - M_G7_V2$S.y)

write.table(M_G7,file= "Data/Genetic Distances/G2/M/M_G7_distances_2.txt")




## Data import ##
M_G8 <- read.table("Data/Genetic Distances/G2/M/M_G8_distances.txt", header = T)
M_G8_old <- read.table("Data/Genetic Distances/0_Old Distance/V2/M distances/M_G8_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
M_G8 <- M_G8[,1:5]

## Set the error value per each tree
Error <- M_G8_old %>% 
  filter(vaccineStrain == unique(M_G8_old$vaccineStrain) & 
           compareStrain == unique(M_G8_old$vaccineStrain))

M_G8_V2 <- M_G8 %>% 
  full_join(Error, by = c("tree.nr"))

M_G8$N2 <- abs(M_G8_V2$N.x - M_G8_V2$N.y)
M_G8$S2 <- abs(M_G8_V2$S.x - M_G8_V2$S.y)

write.table(M_G8,file= "Data/Genetic Distances/G2/M/M_G8_distances_2.txt")





