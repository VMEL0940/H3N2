### Influ H3N2 GenDist Proj ###
## Package ##
library(tidyverse)

## Set working directory ##
setwd("/01_GenDisFlu/")

########### NA !!!!!! #################

## Data import ##
NA_G1 <- read.table("Data/Genetic Distances/NA/NA_G1_distances.txt", header = T)
NA_G1_old <- read.table("Data/Old Distance/NA distances/NA_G1_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
NA_G1 <- NA_G1[,1:5]

## Set the error value per each tree
Error <- NA_G1_old %>% 
  filter(vaccineStrain == unique(NA_G1$vaccineStrain) & 
           compareStrain == unique(NA_G1$vaccineStrain))

NA_G1_V2 <- NA_G1 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

NA_G1$N2 <- abs(NA_G1_V2$N.x - NA_G1_V2$N.y)
NA_G1$S2 <- abs(NA_G1_V2$S.x - NA_G1_V2$S.y)

write.table(NA_G1,file= "Data/Genetic Distances/NA/NA_G1_distances_2.txt")




## Data import ##
NA_G2 <- read.table("Data/Genetic Distances/NA/NA_G2_distances.txt", header = T)
NA_G2_old <- read.table("Data/Old Distance/NA distances/NA_G2_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
NA_G2 <- NA_G2[,1:5]

## Set the error value per each tree
Error <- NA_G2_old %>% 
  filter(vaccineStrain == unique(NA_G2$vaccineStrain) & 
           compareStrain == unique(NA_G2$vaccineStrain))

NA_G2_V2 <- NA_G2 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

NA_G2$N2 <- abs(NA_G2_V2$N.x - NA_G2_V2$N.y)
NA_G2$S2 <- abs(NA_G2_V2$S.x - NA_G2_V2$S.y)

write.table(NA_G2,file= "Data/Genetic Distances/NA/NA_G2_distances_2.txt")





## Data import ##
NA_G3 <- read.table("Data/Genetic Distances/NA/NA_G3_distances.txt", header = T)
NA_G3_old <- read.table("Data/Old Distance/NA distances/NA_G3_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
NA_G3 <- NA_G3[,1:5]



## Set the error value per each tree
Error <- NA_G3_old %>% 
  filter(vaccineStrain == unique(NA_G3$vaccineStrain) & 
           compareStrain == unique(NA_G3$vaccineStrain))

NA_G3_V2 <- NA_G3 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

NA_G3$N2 <- abs(NA_G3_V2$N.x - NA_G3_V2$N.y)
NA_G3$S2 <- abs(NA_G3_V2$S.x - NA_G3_V2$S.y)

write.table(NA_G3,file= "Data/Genetic Distances/NA/NA_G3_distances_2.txt")




## Data import ##
NA_G4 <- read.table("Data/Genetic Distances/NA/NA_G4_distances.txt", header = T)
NA_G4_old <- read.table("Data/Old Distance/NA distances/NA_G4_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
NA_G4 <- NA_G4[,1:5]

## Set the error value per each tree
Error <- NA_G4_old %>% 
  filter(vaccineStrain == unique(NA_G4$vaccineStrain) & 
           compareStrain == unique(NA_G4$vaccineStrain))

NA_G4_V2 <- NA_G4 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

NA_G4$N2 <- abs(NA_G4_V2$N.x - NA_G4_V2$N.y)
NA_G4$S2 <- abs(NA_G4_V2$S.x - NA_G4_V2$S.y)

write.table(NA_G4,file= "Data/Genetic Distances/NA/NA_G4_distances_2.txt")



## Data import ##
NA_G5 <- read.table("Data/Genetic Distances/NA/NA_G5_distances.txt", header = T)
NA_G5_old <- read.table("Data/Old Distance/NA distances/NA_G5_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
NA_G5 <- NA_G5[,1:5]

## Set the error value per each tree
Error <- NA_G5_old %>% 
  filter(vaccineStrain == unique(NA_G5$vaccineStrain) & 
           compareStrain == unique(NA_G5$vaccineStrain))

NA_G5_V2 <- NA_G5 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

NA_G5$N2 <- abs(NA_G5_V2$N.x - NA_G5_V2$N.y)
NA_G5$S2 <- abs(NA_G5_V2$S.x - NA_G5_V2$S.y)

write.table(NA_G5,file= "Data/Genetic Distances/NA/NA_G5_distances_2.txt")





## Data import ##
NA_G6 <- read.table("Data/Genetic Distances/NA/NA_G6_distances.txt", header = T)
NA_G6_old <- read.table("Data/Old Distance/NA distances/NA_G6_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
NA_G6 <- NA_G6[,1:5]

## Set the error value per each tree
Error <- NA_G6_old %>% 
  filter(vaccineStrain == unique(NA_G6$vaccineStrain) & 
           compareStrain == unique(NA_G6$vaccineStrain))

NA_G6_V2 <- NA_G6 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

NA_G6$N2 <- abs(NA_G6_V2$N.x - NA_G6_V2$N.y)
NA_G6$S2 <- abs(NA_G6_V2$S.x - NA_G6_V2$S.y)

write.table(NA_G6,file= "Data/Genetic Distances/NA/NA_G6_distances_2.txt")





## Data import ##
NA_G7 <- read.table("Data/Genetic Distances/NA/NA_G7_distances.txt", header = T)
NA_G7_old <- read.table("Data/Old Distance/NA distances/NA_G7_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
NA_G7 <- NA_G7[,1:5]

## Set the error value per each tree
Error <- NA_G7_old %>% 
  filter(vaccineStrain == unique(NA_G7$vaccineStrain) & 
           compareStrain == unique(NA_G7$vaccineStrain))

NA_G7_V2 <- NA_G7 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

NA_G7$N2 <- abs(NA_G7_V2$N.x - NA_G7_V2$N.y)
NA_G7$S2 <- abs(NA_G7_V2$S.x - NA_G7_V2$S.y)

write.table(NA_G7,file= "Data/Genetic Distances/NA/NA_G7_distances_2.txt")




## Data import ##
NA_G8 <- read.table("Data/Genetic Distances/NA/NA_G8_distances.txt", header = T)
NA_G8_old <- read.table("Data/Old Distance/NA distances/NA_G8_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
NA_G8 <- NA_G8[,1:5]

## Set the error value per each tree
Error <- NA_G8_old %>% 
  filter(vaccineStrain == unique(NA_G8$vaccineStrain) & 
           compareStrain == unique(NA_G8$vaccineStrain))

NA_G8_V2 <- NA_G8 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

NA_G8$N2 <- abs(NA_G8_V2$N.x - NA_G8_V2$N.y)
NA_G8$S2 <- abs(NA_G8_V2$S.x - NA_G8_V2$S.y)

write.table(NA_G8,file= "Data/Genetic Distances/NA/NA_G8_distances_2.txt")





## Data import ##
NA_G9 <- read.table("Data/Genetic Distances/NA/NA_G9_distances.txt", header = T)
NA_G9_old <- read.table("Data/Old Distance/NA distances/NA_G9_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
NA_G9 <- NA_G9[,1:5]

## Set the error value per each tree
Error <- NA_G9_old %>% 
  filter(vaccineStrain == unique(NA_G9$vaccineStrain) & 
           compareStrain == unique(NA_G9$vaccineStrain))

NA_G9_V2 <- NA_G9 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

NA_G9$N2 <- abs(NA_G9_V2$N.x - NA_G9_V2$N.y)
NA_G9$S2 <- abs(NA_G9_V2$S.x - NA_G9_V2$S.y)

write.table(NA_G9,file= "Data/Genetic Distances/NA/NA_G9_distances_2.txt")


