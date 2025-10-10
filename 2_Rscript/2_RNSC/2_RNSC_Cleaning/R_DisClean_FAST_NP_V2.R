### Influ H3N2 GenDist Proj ###
## Package ##
library(tidyverse)

## Set working directory ##
setwd("/01_GenDisFlu/")

########### NP !!!!!! #################

## Data import ##
NP_G1 <- read.table("Data/Genetic Distances/NP/NP_G1_distances.txt", header = T)
NP_G1_old <- read.table("Data/Old Distance/NP distances/NP_G1_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
NP_G1 <- NP_G1[,1:5]

## Set the error value per each tree
Error <- NP_G1_old %>% 
  filter(vaccineStrain == unique(NP_G1$vaccineStrain) & 
           compareStrain == unique(NP_G1$vaccineStrain))

NP_G1_V2 <- NP_G1 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

NP_G1$N2 <- abs(NP_G1_V2$N.x - NP_G1_V2$N.y)
NP_G1$S2 <- abs(NP_G1_V2$S.x - NP_G1_V2$S.y)

write.table(NP_G1,file= "Data/Genetic Distances/NP/NP_G1_distances_2.txt")




## Data import ##
NP_G2 <- read.table("Data/Genetic Distances/NP/NP_G2_distances.txt", header = T)
NP_G2_old <- read.table("Data/Old Distance/NP distances/NP_G2_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
NP_G2 <- NP_G2[,1:5]

## Set the error value per each tree
Error <- NP_G2_old %>% 
  filter(vaccineStrain == unique(NP_G2$vaccineStrain) & 
           compareStrain == unique(NP_G2$vaccineStrain))

NP_G2_V2 <- NP_G2 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

NP_G2$N2 <- abs(NP_G2_V2$N.x - NP_G2_V2$N.y)
NP_G2$S2 <- abs(NP_G2_V2$S.x - NP_G2_V2$S.y)

write.table(NP_G2,file= "Data/Genetic Distances/NP/NP_G2_distances_2.txt")





## Data import ##
NP_G3 <- read.table("Data/Genetic Distances/NP/NP_G3_distances.txt", header = T)
NP_G3_old <- read.table("Data/Old Distance/NP distances/NP_G3_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
NP_G3 <- NP_G3[,1:5]



## Set the error value per each tree
Error <- NP_G3_old %>% 
  filter(vaccineStrain == unique(NP_G3$vaccineStrain) & 
           compareStrain == unique(NP_G3$vaccineStrain))

NP_G3_V2 <- NP_G3 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

NP_G3$N2 <- abs(NP_G3_V2$N.x - NP_G3_V2$N.y)
NP_G3$S2 <- abs(NP_G3_V2$S.x - NP_G3_V2$S.y)

write.table(NP_G3,file= "Data/Genetic Distances/NP/NP_G3_distances_2.txt")




## Data import ##
NP_G4 <- read.table("Data/Genetic Distances/NP/NP_G4_distances.txt", header = T)
NP_G4_old <- read.table("Data/Old Distance/NP distances/NP_G4_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
NP_G4 <- NP_G4[,1:5]

## Set the error value per each tree
Error <- NP_G4_old %>% 
  filter(vaccineStrain == unique(NP_G4$vaccineStrain) & 
           compareStrain == unique(NP_G4$vaccineStrain))

NP_G4_V2 <- NP_G4 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

NP_G4$N2 <- abs(NP_G4_V2$N.x - NP_G4_V2$N.y)
NP_G4$S2 <- abs(NP_G4_V2$S.x - NP_G4_V2$S.y)

write.table(NP_G4,file= "Data/Genetic Distances/NP/NP_G4_distances_2.txt")



## Data import ##
NP_G5 <- read.table("Data/Genetic Distances/NP/NP_G5_distances.txt", header = T)
NP_G5_old <- read.table("Data/Old Distance/NP distances/NP_G5_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
NP_G5 <- NP_G5[,1:5]

## Set the error value per each tree
Error <- NP_G5_old %>% 
  filter(vaccineStrain == unique(NP_G5$vaccineStrain) & 
           compareStrain == unique(NP_G5$vaccineStrain))

NP_G5_V2 <- NP_G5 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

NP_G5$N2 <- abs(NP_G5_V2$N.x - NP_G5_V2$N.y)
NP_G5$S2 <- abs(NP_G5_V2$S.x - NP_G5_V2$S.y)

write.table(NP_G5,file= "Data/Genetic Distances/NP/NP_G5_distances_2.txt")





## Data import ##
NP_G6 <- read.table("Data/Genetic Distances/NP/NP_G6_distances.txt", header = T)
NP_G6_old <- read.table("Data/Old Distance/NP distances/NP_G6_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
NP_G6 <- NP_G6[,1:5]

## Set the error value per each tree
Error <- NP_G6_old %>% 
  filter(vaccineStrain == unique(NP_G6$vaccineStrain) & 
           compareStrain == unique(NP_G6$vaccineStrain))

NP_G6_V2 <- NP_G6 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

NP_G6$N2 <- abs(NP_G6_V2$N.x - NP_G6_V2$N.y)
NP_G6$S2 <- abs(NP_G6_V2$S.x - NP_G6_V2$S.y)

write.table(NP_G6,file= "Data/Genetic Distances/NP/NP_G6_distances_2.txt")





## Data import ##
NP_G7 <- read.table("Data/Genetic Distances/NP/NP_G7_distances.txt", header = T)
NP_G7_old <- read.table("Data/Old Distance/NP distances/NP_G7_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
NP_G7 <- NP_G7[,1:5]

## Set the error value per each tree
Error <- NP_G7_old %>% 
  filter(vaccineStrain == unique(NP_G7$vaccineStrain) & 
           compareStrain == unique(NP_G7$vaccineStrain))

NP_G7_V2 <- NP_G7 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

NP_G7$N2 <- abs(NP_G7_V2$N.x - NP_G7_V2$N.y)
NP_G7$S2 <- abs(NP_G7_V2$S.x - NP_G7_V2$S.y)

write.table(NP_G7,file= "Data/Genetic Distances/NP/NP_G7_distances_2.txt")




## Data import ##
NP_G8 <- read.table("Data/Genetic Distances/NP/NP_G8_distances.txt", header = T)
NP_G8_old <- read.table("Data/Old Distance/NP distances/NP_G8_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
NP_G8 <- NP_G8[,1:5]

## Set the error value per each tree
Error <- NP_G8_old %>% 
  filter(vaccineStrain == unique(NP_G8$vaccineStrain) & 
           compareStrain == unique(NP_G8$vaccineStrain))

NP_G8_V2 <- NP_G8 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

NP_G8$N2 <- abs(NP_G8_V2$N.x - NP_G8_V2$N.y)
NP_G8$S2 <- abs(NP_G8_V2$S.x - NP_G8_V2$S.y)

write.table(NP_G8,file= "Data/Genetic Distances/NP/NP_G8_distances_2.txt")





## Data import ##
NP_G9 <- read.table("Data/Genetic Distances/NP/NP_G9_distances.txt", header = T)
NP_G9_old <- read.table("Data/Old Distance/NP distances/NP_G9_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
NP_G9 <- NP_G9[,1:5]

## Set the error value per each tree
Error <- NP_G9_old %>% 
  filter(vaccineStrain == unique(NP_G9$vaccineStrain) & 
           compareStrain == unique(NP_G9$vaccineStrain))

NP_G9_V2 <- NP_G9 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

NP_G9$N2 <- abs(NP_G9_V2$N.x - NP_G9_V2$N.y)
NP_G9$S2 <- abs(NP_G9_V2$S.x - NP_G9_V2$S.y)

write.table(NP_G9,file= "Data/Genetic Distances/NP/NP_G9_distances_2.txt")


