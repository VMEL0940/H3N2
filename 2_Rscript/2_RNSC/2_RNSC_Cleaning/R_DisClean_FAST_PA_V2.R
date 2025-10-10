### Influ H3N2 GenDist Proj ###
## Package ##
library(tidyverse)

## Set working directory ##
setwd("/01_GenDisFlu/")

########### PA !!!!!! #################

## Data import ##
PA_G1 <- read.table("Data/Genetic Distances/PA/PA_G1_distances.txt", header = T)
PA_G1_old <- read.table("Data/Old Distance/PA distances/PA_G1_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
PA_G1 <- PA_G1[,1:5]

## Set the error value per each tree
Error <- PA_G1_old %>% 
  filter(vaccineStrain == unique(PA_G1$vaccineStrain) & 
           compareStrain == unique(PA_G1$vaccineStrain))

PA_G1_V2 <- PA_G1 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

PA_G1$N2 <- abs(PA_G1_V2$N.x - PA_G1_V2$N.y)
PA_G1$S2 <- abs(PA_G1_V2$S.x - PA_G1_V2$S.y)

write.table(PA_G1,file= "Data/Genetic Distances/PA/PA_G1_distances_2.txt")




## Data import ##
PA_G2 <- read.table("Data/Genetic Distances/PA/PA_G2_distances.txt", header = T)
PA_G2_old <- read.table("Data/Old Distance/PA distances/PA_G2_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
PA_G2 <- PA_G2[,1:5]

## Set the error value per each tree
Error <- PA_G2_old %>% 
  filter(vaccineStrain == unique(PA_G2$vaccineStrain) & 
           compareStrain == unique(PA_G2$vaccineStrain))

PA_G2_V2 <- PA_G2 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

PA_G2$N2 <- abs(PA_G2_V2$N.x - PA_G2_V2$N.y)
PA_G2$S2 <- abs(PA_G2_V2$S.x - PA_G2_V2$S.y)

write.table(PA_G2,file= "Data/Genetic Distances/PA/PA_G2_distances_2.txt")





## Data import ##
PA_G3 <- read.table("Data/Genetic Distances/PA/PA_G3_distances.txt", header = T)
PA_G3_old <- read.table("Data/Old Distance/PA distances/PA_G3_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
PA_G3 <- PA_G3[,1:5]



## Set the error value per each tree
Error <- PA_G3_old %>% 
  filter(vaccineStrain == unique(PA_G3$vaccineStrain) & 
           compareStrain == unique(PA_G3$vaccineStrain))

PA_G3_V2 <- PA_G3 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

PA_G3$N2 <- abs(PA_G3_V2$N.x - PA_G3_V2$N.y)
PA_G3$S2 <- abs(PA_G3_V2$S.x - PA_G3_V2$S.y)

write.table(PA_G3,file= "Data/Genetic Distances/PA/PA_G3_distances_2.txt")




## Data import ##
PA_G4 <- read.table("Data/Genetic Distances/PA/PA_G4_distances.txt", header = T)
PA_G4_old <- read.table("Data/Old Distance/PA distances/PA_G4_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
PA_G4 <- PA_G4[,1:5]

## Set the error value per each tree
Error <- PA_G4_old %>% 
  filter(vaccineStrain == unique(PA_G4$vaccineStrain) & 
           compareStrain == unique(PA_G4$vaccineStrain))

PA_G4_V2 <- PA_G4 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

PA_G4$N2 <- abs(PA_G4_V2$N.x - PA_G4_V2$N.y)
PA_G4$S2 <- abs(PA_G4_V2$S.x - PA_G4_V2$S.y)

write.table(PA_G4,file= "Data/Genetic Distances/PA/PA_G4_distances_2.txt")



## Data import ##
PA_G5 <- read.table("Data/Genetic Distances/PA/PA_G5_distances.txt", header = T)
PA_G5_old <- read.table("Data/Old Distance/PA distances/PA_G5_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
PA_G5 <- PA_G5[,1:5]

## Set the error value per each tree
Error <- PA_G5_old %>% 
  filter(vaccineStrain == unique(PA_G5$vaccineStrain) & 
           compareStrain == unique(PA_G5$vaccineStrain))

PA_G5_V2 <- PA_G5 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

PA_G5$N2 <- abs(PA_G5_V2$N.x - PA_G5_V2$N.y)
PA_G5$S2 <- abs(PA_G5_V2$S.x - PA_G5_V2$S.y)

write.table(PA_G5,file= "Data/Genetic Distances/PA/PA_G5_distances_2.txt")





## Data import ##
PA_G6 <- read.table("Data/Genetic Distances/PA/PA_G6_distances.txt", header = T)
PA_G6_old <- read.table("Data/Old Distance/PA distances/PA_G6_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
PA_G6 <- PA_G6[,1:5]

## Set the error value per each tree
Error <- PA_G6_old %>% 
  filter(vaccineStrain == unique(PA_G6$vaccineStrain) & 
           compareStrain == unique(PA_G6$vaccineStrain))

PA_G6_V2 <- PA_G6 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

PA_G6$N2 <- abs(PA_G6_V2$N.x - PA_G6_V2$N.y)
PA_G6$S2 <- abs(PA_G6_V2$S.x - PA_G6_V2$S.y)

write.table(PA_G6,file= "Data/Genetic Distances/PA/PA_G6_distances_2.txt")





## Data import ##
PA_G7 <- read.table("Data/Genetic Distances/PA/PA_G7_distances.txt", header = T)
PA_G7_old <- read.table("Data/Old Distance/PA distances/PA_G7_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
PA_G7 <- PA_G7[,1:5]

## Set the error value per each tree
Error <- PA_G7_old %>% 
  filter(vaccineStrain == unique(PA_G7$vaccineStrain) & 
           compareStrain == unique(PA_G7$vaccineStrain))

PA_G7_V2 <- PA_G7 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

PA_G7$N2 <- abs(PA_G7_V2$N.x - PA_G7_V2$N.y)
PA_G7$S2 <- abs(PA_G7_V2$S.x - PA_G7_V2$S.y)

write.table(PA_G7,file= "Data/Genetic Distances/PA/PA_G7_distances_2.txt")




## Data import ##
PA_G8 <- read.table("Data/Genetic Distances/PA/PA_G8_distances.txt", header = T)
PA_G8_old <- read.table("Data/Old Distance/PA distances/PA_G8_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
PA_G8 <- PA_G8[,1:5]

## Set the error value per each tree
Error <- PA_G8_old %>% 
  filter(vaccineStrain == unique(PA_G8$vaccineStrain) & 
           compareStrain == unique(PA_G8$vaccineStrain))

PA_G8_V2 <- PA_G8 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

PA_G8$N2 <- abs(PA_G8_V2$N.x - PA_G8_V2$N.y)
PA_G8$S2 <- abs(PA_G8_V2$S.x - PA_G8_V2$S.y)

write.table(PA_G8,file= "Data/Genetic Distances/PA/PA_G8_distances_2.txt")





## Data import ##
PA_G9 <- read.table("Data/Genetic Distances/PA/PA_G9_distances.txt", header = T)
PA_G9_old <- read.table("Data/Old Distance/PA distances/PA_G9_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
PA_G9 <- PA_G9[,1:5]

## Set the error value per each tree
Error <- PA_G9_old %>% 
  filter(vaccineStrain == unique(PA_G9$vaccineStrain) & 
           compareStrain == unique(PA_G9$vaccineStrain))

PA_G9_V2 <- PA_G9 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

PA_G9$N2 <- abs(PA_G9_V2$N.x - PA_G9_V2$N.y)
PA_G9$S2 <- abs(PA_G9_V2$S.x - PA_G9_V2$S.y)

write.table(PA_G9,file= "Data/Genetic Distances/PA/PA_G9_distances_2.txt")


