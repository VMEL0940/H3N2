### Influ H3N2 GenDist Proj ###
## Package ##
library(tidyverse)

## Set working directory ##
setwd("/01_GenDisFlu/")

########### HA !!!!!! #################

## Data import ##
HA_G1 <- read.table("Data/Genetic Distances/HA/HA_G1_distances.txt", header = T)
HA_G1_old <- read.table("Data/Old Distance/HA distances/HA_G1_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
HA_G1 <- HA_G1[,1:5]

unique(HA_G1$vaccineStrain)

## Set the error value per each tree
Error <- HA_G1_old %>% 
  filter(vaccineStrain == "EPI103320_A_Moscow_10_1999_NA_NA_NA_NA" & 
           compareStrain == "EPI103320_A_Moscow_10_1999_NA_NA_NA_NA")

HA_G1_V2 <- HA_G1 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

HA_G1$N2 <- abs(HA_G1_V2$N.x - HA_G1_V2$N.y)
HA_G1$S2 <- abs(HA_G1_V2$S.x - HA_G1_V2$S.y)



write.table(HA_G1,file= "Data/Genetic Distances/HA/HA_G1_distances_2.txt")





## Data import ##
HA_G2 <- read.table("Data/Genetic Distances/HA/HA_G2_distances.txt", header = T)
HA_G2_old <- read.table("Data/Old Distance/HA distances/HA_G2_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
HA_G2 <- HA_G2[,1:5]

unique(HA_G2$vaccineStrain)

## Set the error value per each tree
Error <- HA_G2_old %>% 
  filter(vaccineStrain == "EPI358781_A_Fujian_411_2002_NA_NA_NA_NA" & 
           compareStrain == "EPI358781_A_Fujian_411_2002_NA_NA_NA_NA")

HA_G2_V2 <- HA_G2 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

HA_G2$N2 <- abs(HA_G2_V2$N.x - HA_G2_V2$N.y)
HA_G2$S2 <- abs(HA_G2_V2$S.x - HA_G2_V2$S.y)



write.table(HA_G2,file= "Data/Genetic Distances/HA/HA_G2_distances_2.txt")





## Data import ##
HA_G3 <- read.table("Data/Genetic Distances/HA/HA_G3_distances.txt", header = T)
HA_G3_old <- read.table("Data/Old Distance/HA distances/HA_G3_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
HA_G3 <- HA_G3[,1:5]

unique(HA_G3$vaccineStrain)

## Set the error value per each tree
Error <- HA_G3_old %>% 
  filter(vaccineStrain == "EPI367109_A_California_7_2004_NA_NA_NA_NA" & 
           compareStrain == "EPI367109_A_California_7_2004_NA_NA_NA_NA")

HA_G3_V2 <- HA_G3 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

HA_G3$N2 <- abs(HA_G3_V2$N.x - HA_G3_V2$N.y)
HA_G3$S2 <- abs(HA_G3_V2$S.x - HA_G3_V2$S.y)

write.table(HA_G3,file= "Data/Genetic Distances/HA/HA_G3_distances_2.txt")




## Data import ##
HA_G4 <- read.table("Data/Genetic Distances/HA/HA_G4_distances.txt", header = T)
HA_G4_old <- read.table("Data/Old Distance/HA distances/HA_G4_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
HA_G4 <- HA_G4[,1:5]

unique(HA_G4$vaccineStrain)

## Set the error value per each tree
Error <- HA_G4_old %>% 
  filter(vaccineStrain == unique(HA_G4$vaccineStrain) & 
           compareStrain == unique(HA_G4$vaccineStrain))

HA_G4_V2 <- HA_G4 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

HA_G4$N2 <- abs(HA_G4_V2$N.x - HA_G4_V2$N.y)
HA_G4$S2 <- abs(HA_G4_V2$S.x - HA_G4_V2$S.y)

write.table(HA_G4,file= "Data/Genetic Distances/HA/HA_G4_distances_2.txt")



## Data import ##
HA_G5 <- read.table("Data/Genetic Distances/HA/HA_G5_distances.txt", header = T)
HA_G5_old <- read.table("Data/Old Distance/HA distances/HA_G5_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
HA_G5 <- HA_G5[,1:5]

## Set the error value per each tree
Error <- HA_G5_old %>% 
  filter(vaccineStrain == unique(HA_G5$vaccineStrain) & 
           compareStrain == unique(HA_G5$vaccineStrain))

HA_G5_V2 <- HA_G5 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

HA_G5$N2 <- abs(HA_G5_V2$N.x - HA_G5_V2$N.y)
HA_G5$S2 <- abs(HA_G5_V2$S.x - HA_G5_V2$S.y)

write.table(HA_G5,file= "Data/Genetic Distances/HA/HA_G5_distances_2.txt")





## Data import ##
HA_G6 <- read.table("Data/Genetic Distances/HA/HA_G6_distances.txt", header = T)
HA_G6_old <- read.table("Data/Old Distance/HA distances/HA_G6_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
HA_G6 <- HA_G6[,1:5]

## Set the error value per each tree
Error <- HA_G6_old %>% 
  filter(vaccineStrain == unique(HA_G6$vaccineStrain) & 
           compareStrain == unique(HA_G6$vaccineStrain))

HA_G6_V2 <- HA_G6 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

HA_G6$N2 <- abs(HA_G6_V2$N.x - HA_G6_V2$N.y)
HA_G6$S2 <- abs(HA_G6_V2$S.x - HA_G6_V2$S.y)

write.table(HA_G6,file= "Data/Genetic Distances/HA/HA_G6_distances_2.txt")





## Data import ##
HA_G7 <- read.table("Data/Genetic Distances/HA/HA_G7_distances.txt", header = T)
HA_G7_old <- read.table("Data/Old Distance/HA distances/HA_G7_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
HA_G7 <- HA_G7[,1:5]

## Set the error value per each tree
Error <- HA_G7_old %>% 
  filter(vaccineStrain == unique(HA_G7$vaccineStrain) & 
           compareStrain == unique(HA_G7$vaccineStrain))

HA_G7_V2 <- HA_G7 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

HA_G7$N2 <- abs(HA_G7_V2$N.x - HA_G7_V2$N.y)
HA_G7$S2 <- abs(HA_G7_V2$S.x - HA_G7_V2$S.y)

write.table(HA_G7,file= "Data/Genetic Distances/HA/HA_G7_distances_2.txt")




## Data import ##
HA_G8 <- read.table("Data/Genetic Distances/HA/HA_G8_distances.txt", header = T)
HA_G8_old <- read.table("Data/Old Distance/HA distances/HA_G8_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
HA_G8 <- HA_G8[,1:5]

## Set the error value per each tree
Error <- HA_G8_old %>% 
  filter(vaccineStrain == unique(HA_G8$vaccineStrain) & 
           compareStrain == unique(HA_G8$vaccineStrain))

HA_G8_V2 <- HA_G8 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

HA_G8$N2 <- abs(HA_G8_V2$N.x - HA_G8_V2$N.y)
HA_G8$S2 <- abs(HA_G8_V2$S.x - HA_G8_V2$S.y)

write.table(HA_G8,file= "Data/Genetic Distances/HA/HA_G8_distances_2.txt")





## Data import ##
HA_G9 <- read.table("Data/Genetic Distances/HA/HA_G9_distances.txt", header = T)
HA_G9_old <- read.table("Data/Old Distance/HA distances/HA_G9_distances_2.txt", header = T)

## Extract the "Noise AA sub values ##
## Build redudant column ##
HA_G9 <- HA_G9[,1:5]

## Set the error value per each tree
Error <- HA_G9_old %>% 
  filter(vaccineStrain == unique(HA_G9$vaccineStrain) & 
           compareStrain == unique(HA_G9$vaccineStrain))

HA_G9_V2 <- HA_G9 %>% 
  full_join(Error, by = c("vaccineStrain", "tree.nr"))

HA_G9$N2 <- abs(HA_G9_V2$N.x - HA_G9_V2$N.y)
HA_G9$S2 <- abs(HA_G9_V2$S.x - HA_G9_V2$S.y)

write.table(HA_G9,file= "Data/Genetic Distances/HA/HA_G9_distances_2.txt")


