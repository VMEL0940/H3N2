### Influ H3N2 GenDist Proj -- EValuation ###
## Package ##
library(tidyverse)

## Set working directory ##
setwd("/01_GenDisFlu/Evaluation/")

## Data import ##
HA_G1 <- read.table("RNSC/HA_G1_distances.txt", header = T)
HA_G2 <- read.table("RNSC/HA_G2_distances.txt", header = T)
HA_G3 <- read.table("RNSC/HA_G3_distances.txt", header = T)
HA_G4 <- read.table("RNSC/HA_G4_distances.txt", header = T)
HA_G5 <- read.table("RNSC/HA_G5_distances.txt", header = T)
HA_G6 <- read.table("RNSC/HA_G6_distances.txt", header = T)
HA_G7 <- read.table("RNSC/HA_G7_distances.txt", header = T)
HA_G8 <- read.table("RNSC/HA_G8_distances.txt", header = T)
HA_G9 <- read.table("RNSC/HA_G9_distances.txt", header = T)
HA_G10 <- read.table("RNSC/HA_G10_distances.txt", header = T)
HA_G11 <- read.table("RNSC/HA_G11_distances.txt", header = T)
HA_G12 <- read.table("RNSC/HA_G12_distances.txt", header = T)
HA_G13 <- read.table("RNSC/HA_G13_distances.txt", header = T)
HA_G14 <- read.table("RNSC/HA_G14_distances.txt", header = T)
HA_G15 <- read.table("RNSC/HA_G15_distances.txt", header = T)
HA_G16 <- read.table("RNSC/HA_G16_distances.txt", header = T)
HA_G17 <- read.table("RNSC/HA_G17_distances.txt", header = T)

## Data cleaning  ##
HA_G1 <- HA_G1[1:5] %>% 
  filter(tree.nr != 1)
HA_G2 <- HA_G2[1:5] %>% 
  filter(tree.nr != 1)
HA_G3 <- HA_G3[1:5] %>% 
  filter(tree.nr != 1)
HA_G4 <- HA_G4[1:5] %>% 
  filter(tree.nr != 1)
HA_G5 <- HA_G5[1:5] %>% 
  filter(tree.nr != 1)
HA_G6 <- HA_G6[1:5] %>% 
  filter(tree.nr != 1)
HA_G7 <- HA_G7[1:5] %>% 
  filter(tree.nr != 1)
HA_G8 <- HA_G8[1:5] %>% 
  filter(tree.nr != 1)
HA_G9 <- HA_G9[1:5] %>% 
  filter(tree.nr != 1)
HA_G10 <- HA_G10[1:5] %>% 
  filter(tree.nr != 1)
HA_G11 <- HA_G11[1:5] %>% 
  filter(tree.nr != 1)
HA_G12 <- HA_G12[1:5] %>% 
  filter(tree.nr != 1)
HA_G13 <- HA_G13[1:5] %>% 
  filter(tree.nr != 1)
HA_G14 <- HA_G14[1:5] %>% 
  filter(tree.nr != 1)
HA_G15 <- HA_G15[1:5] %>% 
  filter(tree.nr != 1)
HA_G16 <- HA_G16[1:5] %>% 
  filter(tree.nr != 1)
HA_G17 <- HA_G17[1:5] %>% 
  filter(tree.nr != 1)

## Put the group note ##
HA_G1$Group <- "G01"
HA_G2$Group <- "G02"
HA_G3$Group <- "G03"
HA_G4$Group <- "G04"
HA_G5$Group <- "G05"
HA_G6$Group <- "G06"
HA_G7$Group <- "G07"
HA_G8$Group <- "G08"
HA_G9$Group <- "G09"
HA_G10$Group <- "G10"
HA_G11$Group <- "G11"
HA_G12$Group <- "G12"
HA_G13$Group <- "G13"
HA_G14$Group <- "G14"
HA_G15$Group <- "G15"
HA_G16$Group <- "G16"
HA_G17$Group <- "G17"


## DisCleaning 
## Set the error value per each tree
Error_G1 <- HA_G1 %>% 
  filter(vaccineStrain == compareStrain)

HA_G1_V2 <- HA_G1 %>% 
  full_join(Error_G1, by = c("vaccineStrain", "tree.nr"))

HA_G1$N2 <- abs(HA_G1_V2$N.x - HA_G1_V2$N.y)
HA_G1$S2 <- abs(HA_G1_V2$S.x - HA_G1_V2$S.y)


## Set the error value per each tree
Error_G2 <- HA_G2 %>% 
  filter(vaccineStrain == compareStrain)

HA_G2_V2 <- HA_G2 %>% 
  full_join(Error_G2, by = c("vaccineStrain", "tree.nr"))

HA_G2$N2 <- abs(HA_G2_V2$N.x - HA_G2_V2$N.y)
HA_G2$S2 <- abs(HA_G2_V2$S.x - HA_G2_V2$S.y)


## Set the error value per each tree
Error_G3 <- HA_G3 %>% 
  filter(vaccineStrain == compareStrain)

HA_G3_V2 <- HA_G3 %>% 
  full_join(Error_G3, by = c("vaccineStrain", "tree.nr"))

HA_G3$N2 <- abs(HA_G3_V2$N.x - HA_G3_V2$N.y)
HA_G3$S2 <- abs(HA_G3_V2$S.x - HA_G3_V2$S.y)


## Set the error value per each tree
Error_G4 <- HA_G4 %>% 
  filter(vaccineStrain == compareStrain)

HA_G4_V2 <- HA_G4 %>% 
  full_join(Error_G4, by = c("vaccineStrain", "tree.nr"))

HA_G4$N2 <- abs(HA_G4_V2$N.x - HA_G4_V2$N.y)
HA_G4$S2 <- abs(HA_G4_V2$S.x - HA_G4_V2$S.y)


## Set the error value per each tree
Error_G5 <- HA_G5 %>% 
  filter(vaccineStrain == compareStrain)

HA_G5_V2 <- HA_G5 %>% 
  full_join(Error_G5, by = c("vaccineStrain", "tree.nr"))

HA_G5$N2 <- abs(HA_G5_V2$N.x - HA_G5_V2$N.y)
HA_G5$S2 <- abs(HA_G5_V2$S.x - HA_G5_V2$S.y)


## Set the error value per each tree
Error_G6 <- HA_G6 %>% 
  filter(vaccineStrain == compareStrain)

HA_G6_V2 <- HA_G6 %>% 
  full_join(Error_G6, by = c("vaccineStrain", "tree.nr"))

HA_G6$N2 <- abs(HA_G6_V2$N.x - HA_G6_V2$N.y)
HA_G6$S2 <- abs(HA_G6_V2$S.x - HA_G6_V2$S.y)


## Set the error value per each tree
Error_G7 <- HA_G7 %>% 
  filter(vaccineStrain == compareStrain)

HA_G7_V2 <- HA_G7 %>% 
  full_join(Error_G7, by = c("vaccineStrain", "tree.nr"))

HA_G7$N2 <- abs(HA_G7_V2$N.x - HA_G7_V2$N.y)
HA_G7$S2 <- abs(HA_G7_V2$S.x - HA_G7_V2$S.y)

## Set the error value per each tree
Error_G8 <- HA_G8 %>% 
  filter(vaccineStrain == compareStrain)

HA_G8_V2 <- HA_G8 %>% 
  full_join(Error_G8, by = c("vaccineStrain", "tree.nr"))

HA_G8$N2 <- abs(HA_G8_V2$N.x - HA_G8_V2$N.y)
HA_G8$S2 <- abs(HA_G8_V2$S.x - HA_G8_V2$S.y)

## Set the error value per each tree
Error_G9 <- HA_G9 %>% 
  filter(vaccineStrain == compareStrain)

HA_G9_V2 <- HA_G9 %>% 
  full_join(Error_G9, by = c("vaccineStrain", "tree.nr"))

HA_G9$N2 <- abs(HA_G9_V2$N.x - HA_G9_V2$N.y)
HA_G9$S2 <- abs(HA_G9_V2$S.x - HA_G9_V2$S.y)

## Set the error value per each tree
Error_G10 <- HA_G10 %>% 
  filter(vaccineStrain == compareStrain)

HA_G10_V2 <- HA_G10 %>% 
  full_join(Error_G10, by = c("vaccineStrain", "tree.nr"))

HA_G10$N2 <- abs(HA_G10_V2$N.x - HA_G10_V2$N.y)
HA_G10$S2 <- abs(HA_G10_V2$S.x - HA_G10_V2$S.y)

## Set the error value per each tree
Error_G11 <- HA_G11 %>% 
  filter(vaccineStrain == compareStrain)

HA_G11_V2 <- HA_G11 %>% 
  full_join(Error_G11, by = c("vaccineStrain", "tree.nr"))

HA_G11$N2 <- abs(HA_G11_V2$N.x - HA_G11_V2$N.y)
HA_G11$S2 <- abs(HA_G11_V2$S.x - HA_G11_V2$S.y)

## Set the error value per each tree
Error_G12 <- HA_G12 %>% 
  filter(vaccineStrain == compareStrain)

HA_G12_V2 <- HA_G12 %>% 
  full_join(Error_G12, by = c("vaccineStrain", "tree.nr"))

HA_G12$N2 <- abs(HA_G12_V2$N.x - HA_G12_V2$N.y)
HA_G12$S2 <- abs(HA_G12_V2$S.x - HA_G12_V2$S.y)

## Set the error value per each tree
Error_G13 <- HA_G13 %>% 
  filter(vaccineStrain == compareStrain)

HA_G13_V2 <- HA_G13 %>% 
  full_join(Error_G13, by = c("vaccineStrain", "tree.nr"))

HA_G13$N2 <- abs(HA_G13_V2$N.x - HA_G13_V2$N.y)
HA_G13$S2 <- abs(HA_G13_V2$S.x - HA_G13_V2$S.y)

## Set the error value per each tree
Error_G14 <- HA_G14 %>% 
  filter(vaccineStrain == compareStrain)

HA_G14_V2 <- HA_G14 %>% 
  full_join(Error_G14, by = c("vaccineStrain", "tree.nr"))

HA_G14$N2 <- abs(HA_G14_V2$N.x - HA_G14_V2$N.y)
HA_G14$S2 <- abs(HA_G14_V2$S.x - HA_G14_V2$S.y)

## Set the error value per each tree
Error_G15 <- HA_G15 %>% 
  filter(vaccineStrain == compareStrain)

HA_G15_V2 <- HA_G15 %>% 
  full_join(Error_G15, by = c("vaccineStrain", "tree.nr"))

HA_G15$N2 <- abs(HA_G15_V2$N.x - HA_G15_V2$N.y)
HA_G15$S2 <- abs(HA_G15_V2$S.x - HA_G15_V2$S.y)

## Set the error value per each tree
Error_G16 <- HA_G16 %>% 
  filter(vaccineStrain == compareStrain)

HA_G16_V2 <- HA_G16 %>% 
  full_join(Error_G16, by = c("vaccineStrain", "tree.nr"))

HA_G16$N2 <- abs(HA_G16_V2$N.x - HA_G16_V2$N.y)
HA_G16$S2 <- abs(HA_G16_V2$S.x - HA_G16_V2$S.y)

## Set the error value per each tree
Error_G17 <- HA_G17 %>% 
  filter(vaccineStrain == compareStrain)

HA_G17_V2 <- HA_G17 %>% 
  full_join(Error_G17, by = c("vaccineStrain", "tree.nr"))

HA_G17$N2 <- abs(HA_G17_V2$N.x - HA_G17_V2$N.y)
HA_G17$S2 <- abs(HA_G17_V2$S.x - HA_G17_V2$S.y)

## Remove "_NA"
HA_G1$vaccineStrain <- gsub('_NA', '', HA_G1$vaccineStrain)
HA_G1$compareStrain <- gsub('_NA', '', HA_G1$compareStrain)

HA_G2$vaccineStrain <- gsub('_NA', '', HA_G2$vaccineStrain)
HA_G2$compareStrain <- gsub('_NA', '', HA_G2$compareStrain)

HA_G3$vaccineStrain <- gsub('_NA', '', HA_G3$vaccineStrain)
HA_G3$compareStrain <- gsub('_NA', '', HA_G3$compareStrain)

HA_G4$vaccineStrain <- gsub('_NA', '', HA_G4$vaccineStrain)
HA_G4$compareStrain <- gsub('_NA', '', HA_G4$compareStrain)

HA_G5$vaccineStrain <- gsub('_NA', '', HA_G5$vaccineStrain)
HA_G5$compareStrain <- gsub('_NA', '', HA_G5$compareStrain)

HA_G6$vaccineStrain <- gsub('_NA', '', HA_G6$vaccineStrain)
HA_G6$compareStrain <- gsub('_NA', '', HA_G6$compareStrain)

HA_G7$vaccineStrain <- gsub('_NA', '', HA_G7$vaccineStrain)
HA_G7$compareStrain <- gsub('_NA', '', HA_G7$compareStrain)

HA_G8$vaccineStrain <- gsub('_NA', '', HA_G8$vaccineStrain)
HA_G8$compareStrain <- gsub('_NA', '', HA_G8$compareStrain)

HA_G9$vaccineStrain <- gsub('_NA', '', HA_G9$vaccineStrain)
HA_G9$compareStrain <- gsub('_NA', '', HA_G9$compareStrain)

HA_G10$vaccineStrain <- gsub('_NA', '', HA_G10$vaccineStrain)
HA_G10$compareStrain <- gsub('_NA', '', HA_G10$compareStrain)

HA_G11$vaccineStrain <- gsub('_NA', '', HA_G11$vaccineStrain)
HA_G11$compareStrain <- gsub('_NA', '', HA_G11$compareStrain)

HA_G12$vaccineStrain <- gsub('_NA', '', HA_G12$vaccineStrain)
HA_G12$compareStrain <- gsub('_NA', '', HA_G12$compareStrain)

HA_G13$vaccineStrain <- gsub('_NA', '', HA_G13$vaccineStrain)
HA_G13$compareStrain <- gsub('_NA', '', HA_G13$compareStrain)

HA_G14$vaccineStrain <- gsub('_NA', '', HA_G14$vaccineStrain)
HA_G14$compareStrain <- gsub('_NA', '', HA_G14$compareStrain)

HA_G15$vaccineStrain <- gsub('_NA', '', HA_G15$vaccineStrain)
HA_G15$compareStrain <- gsub('_NA', '', HA_G15$compareStrain)

HA_G16$vaccineStrain <- gsub('_NA', '', HA_G16$vaccineStrain)
HA_G16$compareStrain <- gsub('_NA', '', HA_G16$compareStrain)

HA_G17$vaccineStrain <- gsub('_NA', '', HA_G17$vaccineStrain)
HA_G17$compareStrain <- gsub('_NA', '', HA_G17$compareStrain)

## Data filter -- Remove Vaccine == Compare ##
HA_G1 <- HA_G1 %>% 
  filter(vaccineStrain != compareStrain)

HA_G2 <- HA_G2 %>% 
  filter(vaccineStrain != compareStrain)

HA_G3 <- HA_G3 %>% 
  filter(vaccineStrain != compareStrain)

HA_G4 <- HA_G4 %>% 
  filter(vaccineStrain != compareStrain)

HA_G5 <- HA_G5 %>% 
  filter(vaccineStrain != compareStrain)

HA_G6 <- HA_G6 %>% 
  filter(vaccineStrain != compareStrain)

HA_G7 <- HA_G7 %>% 
  filter(vaccineStrain != compareStrain)

HA_G8 <- HA_G8 %>% 
  filter(vaccineStrain != compareStrain)

HA_G9 <- HA_G9 %>% 
  filter(vaccineStrain != compareStrain)

HA_G10 <- HA_G10 %>% 
  filter(vaccineStrain != compareStrain)

HA_G11 <- HA_G11 %>% 
  filter(vaccineStrain != compareStrain)

HA_G12 <- HA_G12 %>% 
  filter(vaccineStrain != compareStrain)

HA_G13 <- HA_G13 %>% 
  filter(vaccineStrain != compareStrain)

HA_G14 <- HA_G14 %>% 
  filter(vaccineStrain != compareStrain)

HA_G15 <- HA_G15 %>% 
  filter(vaccineStrain != compareStrain)

HA_G16 <- HA_G16 %>% 
  filter(vaccineStrain != compareStrain)

HA_G17 <- HA_G17 %>% 
  filter(vaccineStrain != compareStrain)


## Summarize the 901 tree values
G1_dist <- HA_G1 %>% 
  group_by(vaccineStrain, compareStrain, Group) %>% 
  summarise(Nmean = mean(N2), Nmedian = median(N2), Nerror = mean(N2)-median(N2),
            Smean = mean(S2), Smedian = median(S2), Serror = mean(S2)-median(S2))

G2_dist <- HA_G2 %>% 
  group_by(vaccineStrain, compareStrain, Group) %>% 
  summarise(Nmean = mean(N2), Nmedian = median(N2), Nerror = mean(N2)-median(N2),
            Smean = mean(S2), Smedian = median(S2), Serror = mean(S2)-median(S2))

G3_dist <- HA_G3 %>% 
  group_by(vaccineStrain, compareStrain, Group) %>% 
  summarise(Nmean = mean(N2), Nmedian = median(N2), Nerror = mean(N2)-median(N2),
            Smean = mean(S2), Smedian = median(S2), Serror = mean(S2)-median(S2))

G4_dist <- HA_G4 %>% 
  group_by(vaccineStrain, compareStrain, Group) %>% 
  summarise(Nmean = mean(N2), Nmedian = median(N2), Nerror = mean(N2)-median(N2),
            Smean = mean(S2), Smedian = median(S2), Serror = mean(S2)-median(S2))

G5_dist <- HA_G5 %>% 
  group_by(vaccineStrain, compareStrain, Group) %>% 
  summarise(Nmean = mean(N2), Nmedian = median(N2), Nerror = mean(N2)-median(N2),
            Smean = mean(S2), Smedian = median(S2), Serror = mean(S2)-median(S2))

G6_dist <- HA_G6 %>% 
  group_by(vaccineStrain, compareStrain, Group) %>% 
  summarise(Nmean = mean(N2), Nmedian = median(N2), Nerror = mean(N2)-median(N2),
            Smean = mean(S2), Smedian = median(S2), Serror = mean(S2)-median(S2))

G7_dist <- HA_G7 %>% 
  group_by(vaccineStrain, compareStrain, Group) %>% 
  summarise(Nmean = mean(N2), Nmedian = median(N2), Nerror = mean(N2)-median(N2),
            Smean = mean(S2), Smedian = median(S2), Serror = mean(S2)-median(S2))

G8_dist <- HA_G8 %>% 
  group_by(vaccineStrain, compareStrain, Group) %>% 
  summarise(Nmean = mean(N2), Nmedian = median(N2), Nerror = mean(N2)-median(N2),
            Smean = mean(S2), Smedian = median(S2), Serror = mean(S2)-median(S2))

G9_dist <- HA_G9 %>% 
  group_by(vaccineStrain, compareStrain, Group) %>% 
  summarise(Nmean = mean(N2), Nmedian = median(N2), Nerror = mean(N2)-median(N2),
            Smean = mean(S2), Smedian = median(S2), Serror = mean(S2)-median(S2))

G10_dist <- HA_G10 %>% 
  group_by(vaccineStrain, compareStrain, Group) %>% 
  summarise(Nmean = mean(N2), Nmedian = median(N2), Nerror = mean(N2)-median(N2),
            Smean = mean(S2), Smedian = median(S2), Serror = mean(S2)-median(S2))

G11_dist <- HA_G11 %>% 
  group_by(vaccineStrain, compareStrain, Group) %>% 
  summarise(Nmean = mean(N2), Nmedian = median(N2), Nerror = mean(N2)-median(N2),
            Smean = mean(S2), Smedian = median(S2), Serror = mean(S2)-median(S2))

G12_dist <- HA_G12 %>% 
  group_by(vaccineStrain, compareStrain, Group) %>% 
  summarise(Nmean = mean(N2), Nmedian = median(N2), Nerror = mean(N2)-median(N2),
            Smean = mean(S2), Smedian = median(S2), Serror = mean(S2)-median(S2))

G13_dist <- HA_G13 %>% 
  group_by(vaccineStrain, compareStrain, Group) %>% 
  summarise(Nmean = mean(N2), Nmedian = median(N2), Nerror = mean(N2)-median(N2),
            Smean = mean(S2), Smedian = median(S2), Serror = mean(S2)-median(S2))

G14_dist <- HA_G14 %>% 
  group_by(vaccineStrain, compareStrain, Group) %>% 
  summarise(Nmean = mean(N2), Nmedian = median(N2), Nerror = mean(N2)-median(N2),
            Smean = mean(S2), Smedian = median(S2), Serror = mean(S2)-median(S2))

G15_dist <- HA_G15 %>% 
  group_by(vaccineStrain, compareStrain, Group) %>% 
  summarise(Nmean = mean(N2), Nmedian = median(N2), Nerror = mean(N2)-median(N2),
            Smean = mean(S2), Smedian = median(S2), Serror = mean(S2)-median(S2))

G16_dist <- HA_G16 %>% 
  group_by(vaccineStrain, compareStrain, Group) %>% 
  summarise(Nmean = mean(N2), Nmedian = median(N2), Nerror = mean(N2)-median(N2),
            Smean = mean(S2), Smedian = median(S2), Serror = mean(S2)-median(S2))

G17_dist <- HA_G17 %>% 
  group_by(vaccineStrain, compareStrain, Group) %>% 
  summarise(Nmean = mean(N2), Nmedian = median(N2), Nerror = mean(N2)-median(N2),
            Smean = mean(S2), Smedian = median(S2), Serror = mean(S2)-median(S2))

HA242_NSdist <- rbind(G1_dist, G2_dist, G3_dist, G4_dist, G5_dist,
                      G6_dist, G7_dist, G8_dist, G9_dist, G10_dist,
                      G11_dist, G12_dist, G13_dist, G14_dist, G15_dist,
                      G16_dist, G17_dist) 


HA242_NSdist$Year <- as.numeric(str_sub(HA242_NSdist$compareStrain, start= -4))

##Export Data 
write.csv(HA242_NSdist,file="HA242_NSdist.csv")



