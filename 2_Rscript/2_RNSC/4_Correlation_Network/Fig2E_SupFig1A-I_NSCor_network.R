### Influ H3N2 GenDist Proj ###
## Package ##
library(tidyverse)
library(factoextra)
library(FactoMineR)
library(ggpubr)
library(corrplot)
library(igraph)
library(Hmisc)
library(dplyr)
library(tidyr)  
library(svglite)

## custom function from Corr matrix to flat df

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}


## Set working directory ##
setwd("~/H3N2/")

## Data import ##
G1_NS <- read.table("1_Data/2_RNSC/3_DistanceLog/1_Before2020_601/Gs/G1_NS_2.txt", header = T, na.strings = "")
G2_NS <- read.table("1_Data/2_RNSC/3_DistanceLog/1_Before2020_601/Gs/G2_NS_2.txt", header = T, na.strings = "")
G3_NS <- read.table("1_Data/2_RNSC/3_DistanceLog/1_Before2020_601/Gs/G3_NS_2.txt", header = T, na.strings = "")
G4_NS <- read.table("1_Data/2_RNSC/3_DistanceLog/1_Before2020_601/Gs/G4_NS_2.txt", header = T, na.strings = "")
G5_NS <- read.table("1_Data/2_RNSC/3_DistanceLog/1_Before2020_601/Gs/G5_NS_2.txt", header = T, na.strings = "")
G6_NS <- read.table("1_Data/2_RNSC/3_DistanceLog/1_Before2020_601/Gs/G6_NS_2.txt", header = T, na.strings = "")
G7_NS <- read.table("1_Data/2_RNSC/3_DistanceLog/1_Before2020_601/Gs/G7_NS_2.txt", header = T, na.strings = "")
G8_NS <- read.table("1_Data/2_RNSC/3_DistanceLog/1_Before2020_601/Gs/G8_NS_2.txt", header = T, na.strings = "")
G9_NS <- read.table("1_Data/2_RNSC/3_DistanceLog/1_Before2020_601/Gs/G9_NS_2.txt", header = T, na.strings = "")

G1_NS$Gene <- factor(G1_NS$Gene, levels = c("PB2", "PB1", "PA", "HA", 
                                               "NP", "NA", "M", "NS"))
G2_NS$Gene <- factor(G2_NS$Gene, levels = c("PB2", "PB1", "PA", "HA", 
                                            "NP", "NA", "M", "NS"))
G3_NS$Gene <- factor(G3_NS$Gene, levels = c("PB2", "PB1", "PA", "HA", 
                                            "NP", "NA", "M", "NS"))
G4_NS$Gene <- factor(G4_NS$Gene, levels = c("PB2", "PB1", "PA", "HA", 
                                            "NP", "NA", "M", "NS"))
G5_NS$Gene <- factor(G5_NS$Gene, levels = c("PB2", "PB1", "PA", "HA", 
                                            "NP", "NA", "M", "NS"))
G6_NS$Gene <- factor(G6_NS$Gene, levels = c("PB2", "PB1", "PA", "HA", 
                                            "NP", "NA", "M", "NS"))
G7_NS$Gene <- factor(G7_NS$Gene, levels = c("PB2", "PB1", "PA", "HA", 
                                            "NP", "NA", "M", "NS"))
G8_NS$Gene <- factor(G8_NS$Gene, levels = c("PB2", "PB1", "PA", "HA", 
                                            "NP", "NA", "M", "NS"))
G9_NS$Gene <- factor(G9_NS$Gene, levels = c("PB2", "PB1", "PA", "HA", 
                                            "NP", "NA", "M", "NS"))

GeneCol <- c("#C780E8", "#9D94FF","#59ADF6","#08CAD1","#42D6A4","#F8F38D","#FFB480","#FF6961")


output_dir <- "5_SFigure/Sfig1_network_svg"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

gene_order <- c("PB2","PB1","PA","HA","NP","NA","M","NS")

for (i in 1:9) {
  message("Processing G", i, " ...")
  df_name <- paste0("G", i, "_NS")
  if (!exists(df_name)) {
    warning(df_name, " No data ")
    next
  }
  df <- get(df_name)
  
  df <- df %>%
    group_by(compareStrain, Gene) %>%
    mutate(Nmedian = mean(N2, na.rm = TRUE)) %>%
    select(vaccineStrain, compareStrain, Gene, Nmedian) %>%
    distinct()
  
  N_dist <- tidyr::spread(df, key = Gene, value = Nmedian)
  
  keep_cols <- intersect(colnames(N_dist), gene_order)
  
  res2 <- rcorr(as.matrix(N_dist[, keep_cols]))
  CorMat <- flattenCorrMatrix(res2$r, res2$P)
  colnames(CorMat) <- c("Gene1", "Gene2", "Correlation", "p")
  
  CorMat <- CorMat %>%
    filter(Correlation != 1, p < 0.05)
  
  CorMat03 <- CorMat %>% filter(Correlation >= 0.3)
  
  actors <- data.frame(name = gene_order)
  relations <- data.frame(
    from = CorMat03$Gene1,
    to   = CorMat03$Gene2,
    weight = CorMat03$Correlation
  )
  
  g <- graph_from_data_frame(relations, directed = FALSE, vertices = actors)
  
  if (ecount(g) > 0) {
    E(g)$width <- E(g)$weight * 2.5^2
    E(g)[E(g)$weight >= 0.3 & E(g)$weight < 0.595]$color <- "#ABB2B9"
    E(g)[E(g)$weight >= 0.595]$color <- "#F5B041"
  }
  
  set.seed(123)
  lay <- layout_with_fr(g)
  
  svg_file <- file.path(output_dir, sprintf("G%d_network.svg", i))
  svglite::svglite(svg_file, width = 7, height = 5)  # svg() 대신 이 줄
  plot(
    g,
    layout = lay,
    vertex.size = 30,
    vertex.color = GeneCol,
    vertex.label.color = "black",
    vertex.label.family = "sans",
    vertex.frame.width = 0.75
  )
  dev.off()
  
  message("Saved: ", svg_file)
}


########## for G1 ###########

## Let's use Median!! ##
G1_NS <- G1_NS %>%
  group_by(compareStrain, Gene) %>%
  mutate(Nmedian = mean(N2, na.rm = TRUE))

G1_NS <- G1_NS %>% select(vaccineStrain, compareStrain, Gene, Nmedian) %>% unique

## Build Correlation matrix - G1##
N_G1_Ndist_spread <- spread(G1_NS, key = Gene, value = Nmedian) 

## Correlation matrix - r and p
res2 <- rcorr(as.matrix(N_G1_Ndist_spread[,c(3:10)]))
G1_NS_Cor <- flattenCorrMatrix(res2$r, res2$P)

colnames(G1_NS_Cor) <- c("Gene1", "Gene2", "Correlation", "p")

G1_NS_Cor <- G1_NS_Cor %>% 
  filter(Correlation != 1 & p < 0.05) 

G1_NS_Cor$GenePair <- paste(G1_NS_Cor$Gene1, "-", G1_NS_Cor$Gene2)

G1_NS_Cor %>% 
  ggplot(aes(x=reorder(GenePair, Correlation), y=Correlation)) + 
  geom_bar(stat = "identity") +
  geom_hline(yintercept= 0.3, linetype=4, color='gray', size=0.75) +
  geom_hline(yintercept= -0.3, linetype=4, color='gray', size=0.75) +
  geom_hline(yintercept= 0.6, linetype=3, color='orange2', size=1) +
  geom_hline(yintercept= -0.6, linetype=3, color='orange2', size=1) +
  ylim(-1,1) +
  xlab("Gene Pair") +
  coord_flip() +
  theme_bw()

## Filter Corr > 0.3 ##
G1_NS_Cor_03 <- G1_NS_Cor %>% 
  filter(Correlation >= 0.3) 

actors <- data.frame(name=c("PB2", "PB1", "PA", "HA", 
                            "NP", "NA", "M", "NS"))
relations <- data.frame(from= G1_NS_Cor_03$Gene1,
                        to= G1_NS_Cor_03$Gene2,
                        weight = G1_NS_Cor_03$Correlation)

g <- graph_from_data_frame(relations, directed=F, vertices=actors)

E(g)$width <- E(g)$weight*2.5^2

E(g)[E(g)$weight >= 0.3 & E(g)$weight < 0.595]$color <- "#ABB2B9" 
E(g)[E(g)$weight >= 0.595]$color <- "#F5B041" 


plot(g, vertex.size=30, vertex.color=GeneCol, vertex.label.color="black", vertex.label.family = "sans",
     vertex.frame.width = 0.75)

########## for G2 ###########

## Let's use Median!! ##
G2_NS <- G2_NS %>%
  group_by(compareStrain, Gene) %>%
  mutate(Nmedian = mean(N2, na.rm = TRUE))

G2_NS <- G2_NS %>% select(vaccineStrain, compareStrain, Gene, Nmedian) %>% unique

## Build Correlation matrix - G2##
N_G2_Ndist_spread <- spread(G2_NS, key = Gene, value = Nmedian) 

## Correlation matrix - r and p
res2 <- rcorr(as.matrix(N_G2_Ndist_spread[,c(3:10)]))
G2_NS_Cor <- flattenCorrMatrix(res2$r, res2$P)

colnames(G2_NS_Cor) <- c("Gene1", "Gene2", "Correlation", "p")

G2_NS_Cor <- G2_NS_Cor %>% 
  filter(Correlation != 1 & p < 0.05) 

G2_NS_Cor$GenePair <- paste(G2_NS_Cor$Gene1, "-", G2_NS_Cor$Gene2)

G2_NS_Cor %>% 
  ggplot(aes(x=reorder(GenePair, Correlation), y=Correlation)) + 
  geom_bar(stat = "identity") +
  geom_hline(yintercept= 0.3, linetype=4, color='gray', size=0.75) +
  geom_hline(yintercept= -0.3, linetype=4, color='gray', size=0.75) +
  geom_hline(yintercept= 0.6, linetype=3, color='orange2', size=1) +
  geom_hline(yintercept= -0.6, linetype=3, color='orange2', size=1) +
    ylim(-1,1) +
  xlab("Gene Pair") +
  coord_flip() +
  theme_bw()

## Filter Corr > 0.3 ##
G2_NS_Cor_03 <- G2_NS_Cor %>% 
  filter(Correlation >= 0.3) 

actors <- data.frame(name=c("PB2", "PB1", "PA", "HA", 
                            "NP", "NA", "M", "NS"))
relations <- data.frame(from= G2_NS_Cor_03$Gene1,
                        to= G2_NS_Cor_03$Gene2,
                        weight = G2_NS_Cor_03$Correlation)

g <- graph_from_data_frame(relations, directed=F, vertices=actors)

E(g)$width <- E(g)$weight*2.5^2

E(g)[E(g)$weight >= 0.3 & E(g)$weight < 0.595]$color <- "#ABB2B9" 
E(g)[E(g)$weight >= 0.595]$color <- "#F5B041" 


plot(g, vertex.size=30, vertex.color=GeneCol, vertex.label.color="black", vertex.label.family = "sans",
     vertex.frame.width = 0.75)

########## for G3 ###########

## Let's use Median!! ##
G3_NS <- G3_NS %>%
  group_by(compareStrain, Gene) %>%
  mutate(Nmedian = mean(N2, na.rm = TRUE))

G3_NS <- G3_NS %>% select(vaccineStrain, compareStrain, Gene, Nmedian) %>% unique

## Build Correlation matrix - G3##
N_G3_Ndist_spread <- spread(G3_NS, key = Gene, value = Nmedian) 

## Correlation matrix - r and p
res2 <- rcorr(as.matrix(N_G3_Ndist_spread[,c(3:10)]))

## Build Correlation matrix - G3##
N_G3_Ndist_spread <- spread(G3_NS, key = Gene, value = Nmedian) 

res2<-rcorr(as.matrix(N_G3_Ndist_spread[,c(3:10)]))
G3_NS_Cor <- flattenCorrMatrix(res2$r, res2$P)

colnames(G3_NS_Cor) <- c("Gene1", "Gene2", "Correlation", "p")

G3_NS_Cor <- G3_NS_Cor %>% 
  filter(Correlation != 1 & p < 0.05) 

G3_NS_Cor$GenePair <- paste(G3_NS_Cor$Gene1, "-", G3_NS_Cor$Gene2)


G3_NS_Cor %>% 
  ggplot(aes(x=reorder(GenePair, Correlation), y=Correlation)) + 
  geom_bar(stat = "identity") +
  geom_hline(yintercept= 0.3, linetype=4, color='gray', size=0.75) +
  geom_hline(yintercept= -0.3, linetype=4, color='gray', size=0.75) +
  geom_hline(yintercept= 0.6, linetype=3, color='orange2', size=1) +
  geom_hline(yintercept= -0.6, linetype=3, color='orange2', size=1) +
  ylim(-1,1) +
  xlab("Gene Pair") +
  coord_flip() +
  theme_bw()

## Filter Corr > 0.3 ##

G3_NS_Cor_03 <- G3_NS_Cor %>% 
  filter(Correlation >= 0.3) 

actors <- data.frame(name=c("PB2", "PB1", "PA", "HA", 
                            "NP", "NA", "M", "NS"))
relations <- data.frame(from= G3_NS_Cor_03$Gene1,
                        to= G3_NS_Cor_03$Gene2,
                        weight = G3_NS_Cor_03$Correlation)

g <- graph_from_data_frame(relations, directed=F, vertices=actors)

E(g)$width <- E(g)$weight*2.5^2
E(g)[E(g)$weight >= 0.3 & E(g)$weight < 0.595]$color <- "#ABB2B9" 
E(g)[E(g)$weight >= 0.595]$color <- "#F5B041"  

as <- authority_score(g, weights=E(g)$weight)$vector

plot(g, vertex.size=30, vertex.color=GeneCol, vertex.label.color="black", vertex.label.family = "sans",
     vertex.frame.width = 0.75)


########## for G4 ###########

## Let's use Median!! ##
G4_NS <- G4_NS %>%
  group_by(compareStrain, Gene) %>%
  mutate(Nmedian = mean(N2, na.rm = TRUE))

G4_NS <- G4_NS %>% select(vaccineStrain, compareStrain, Gene, Nmedian) %>% unique

## Build Correlation matrix - G4##
N_G4_Ndist_spread <- spread(G4_NS, key = Gene, value = Nmedian) 

res2<-rcorr(as.matrix(N_G4_Ndist_spread[,c(3:10)]))
G4_NS_Cor <- flattenCorrMatrix(res2$r, res2$P)

colnames(G4_NS_Cor) <- c("Gene1", "Gene2", "Correlation", "p")

G4_NS_Cor <- G4_NS_Cor %>% 
  filter(Correlation != 1 & p < 0.05) 

G4_NS_Cor$GenePair <- paste(G4_NS_Cor$Gene1, "-", G4_NS_Cor$Gene2)

G4_NS_Cor %>% 
  ggplot(aes(x=reorder(GenePair, Correlation), y=Correlation)) + 
  geom_bar(stat = "identity") +
  geom_hline(yintercept= 0.3, linetype=4, color='gray', size=0.75) +
  geom_hline(yintercept= -0.3, linetype=4, color='gray', size=0.75) +
  geom_hline(yintercept= 0.6, linetype=3, color='orange2', size=1) +
  geom_hline(yintercept= -0.6, linetype=3, color='orange2', size=1) +
  ylim(-1,1) +
  xlab("Gene Pair") +
  coord_flip() +
  theme_bw()

## Filter Corr > 0.3 ##

G4_NS_Cor_03 <- G4_NS_Cor %>% 
  filter(Correlation >= 0.3) 

actors <- data.frame(name=c("PB2", "PB1", "PA", "HA", 
                            "NP", "NA", "M", "NS"))
relations <- data.frame(from= G4_NS_Cor_03$Gene1,
                        to= G4_NS_Cor_03$Gene2,
                        weight = G4_NS_Cor_03$Correlation)

g <- graph_from_data_frame(relations, directed=F, vertices=actors)

E(g)$width <- E(g)$weight*2.5^2

E(g)[E(g)$weight >= 0.3 & E(g)$weight < 0.595]$color <- "#ABB2B9" 
E(g)[E(g)$weight >= 0.595]$color <- "#F5B041" 

as <- authority_score(g, weights=E(g)$weight)$vector

plot(g, vertex.size=30, vertex.color=GeneCol, vertex.label.color="black", vertex.label.family = "sans",
     vertex.frame.width = 0.75)


########## for G5 ###########

## Let's use Median!! ##
G5_NS <- G5_NS %>%
  group_by(compareStrain, Gene) %>%
  mutate(Nmedian = mean(N2, na.rm = TRUE))

G5_NS <- G5_NS %>% select(vaccineStrain, compareStrain, Gene, Nmedian) %>% unique

## Spread the table!! ##
N_G5_Ndist_spread <- spread(G5_NS, key = Gene, value = Nmedian) 

## Build Correlation matrix - G5##
res2<-rcorr(as.matrix(N_G5_Ndist_spread[,c(3:10)]))
G5_NS_Cor <- flattenCorrMatrix(res2$r, res2$P)

colnames(G5_NS_Cor) <- c("Gene1", "Gene2", "Correlation", "p")

G5_NS_Cor <- G5_NS_Cor %>% 
  filter(Correlation != 1 & p < 0.05) 

G5_NS_Cor$GenePair <- paste(G5_NS_Cor$Gene1, "-", G5_NS_Cor$Gene2)

G5_NS_Cor %>% 
  ggplot(aes(x=reorder(GenePair, Correlation), y=Correlation)) + 
  geom_bar(stat = "identity") +
  geom_hline(yintercept= 0.3, linetype=4, color='gray', size=0.75) +
  geom_hline(yintercept= -0.3, linetype=4, color='gray', size=0.75) +
  geom_hline(yintercept= 0.6, linetype=3, color='orange2', size=1) +
  geom_hline(yintercept= -0.6, linetype=3, color='orange2', size=1) +
  ylim(-1,1) +
  xlab("Gene Pair") +
  coord_flip() +
  theme_bw()

## Filter Corr > 0.3 ##

G5_NS_Cor_03 <- G5_NS_Cor %>% 
  filter(Correlation >= 0.3) 

actors <- data.frame(name=c("PB2", "PB1", "PA", "HA", 
                            "NP", "NA", "M", "NS"))
relations <- data.frame(from= G5_NS_Cor_03$Gene1,
                        to= G5_NS_Cor_03$Gene2,
                        weight = G5_NS_Cor_03$Correlation)

g <- graph_from_data_frame(relations, directed=F, vertices=actors)

E(g)$width <- E(g)$weight*2.5^2


E(g)[E(g)$weight >= 0.3 & E(g)$weight < 0.595]$color <- "#ABB2B9" 
E(g)[E(g)$weight >= 0.595]$color <- "#F5B041" 

as <- authority_score(g, weights=E(g)$weight)$vector

plot(g, vertex.size=30, vertex.color=GeneCol, vertex.label.color="black", vertex.label.family = "sans",
     vertex.frame.width = 0.75)


########## for G6 ###########

## Let's use Median!! ##
G6_NS <- G6_NS %>%
  group_by(compareStrain, Gene) %>%
  mutate(Nmedian = mean(N2, na.rm = TRUE))

G6_NS <- G6_NS %>% select(vaccineStrain, compareStrain, Gene, Nmedian) %>% unique

## Spread the table!! ##
N_G6_Ndist_spread <- spread(G6_NS, key = Gene, value = Nmedian) 

## Build Correlation matrix - G6##
res2<-rcorr(as.matrix(N_G6_Ndist_spread[,c(3:10)]))
G6_NS_Cor <- flattenCorrMatrix(res2$r, res2$P)

colnames(G6_NS_Cor) <- c("Gene1", "Gene2", "Correlation", "p")

G6_NS_Cor <- G6_NS_Cor %>% 
  filter(Correlation != 1 & p < 0.05) 

G6_NS_Cor$GenePair <- paste(G6_NS_Cor$Gene1, "-", G6_NS_Cor$Gene2)

G6_NS_Cor %>% 
  ggplot(aes(x=reorder(GenePair, Correlation), y=Correlation)) + 
  geom_bar(stat = "identity") +
  geom_hline(yintercept= 0.3, linetype=4, color='gray', size=0.75) +
  geom_hline(yintercept= -0.3, linetype=4, color='gray', size=0.75) +
  geom_hline(yintercept= 0.6, linetype=3, color='orange2', size=1) +
  geom_hline(yintercept= -0.6, linetype=3, color='orange2', size=1) +
  ylim(-1,1) +
  xlab("Gene Pair") +
  coord_flip() +
  theme_bw()

## Filter Corr > 0.3 ##

G6_NS_Cor_03 <- G6_NS_Cor %>% 
  filter(Correlation >= 0.3) 

actors <- data.frame(name=c("PB2", "PB1", "PA", "HA", 
                            "NP", "NA", "M", "NS"))
relations <- data.frame(from= G6_NS_Cor_03$Gene1,
                        to= G6_NS_Cor_03$Gene2,
                        weight = G6_NS_Cor_03$Correlation)

g <- graph_from_data_frame(relations, directed=F, vertices=actors)

E(g)$width <- E(g)$weight*2.5^2

E(g)[E(g)$weight >= 0.3 & E(g)$weight < 0.595]$color <- "#ABB2B9" 
E(g)[E(g)$weight >= 0.595]$color <- "#F5B041" 

as <- authority_score(g, weights=E(g)$weight)$vector

plot(g, vertex.size=30, vertex.color=GeneCol, vertex.label.color="black", vertex.label.family = "sans",
     vertex.frame.width = 0.75)


########## for G7 ###########

## Let's use Median!! ##
G7_NS <- G7_NS %>%
  group_by(compareStrain, Gene) %>%
  mutate(Nmedian = mean(N2, na.rm = TRUE))

G7_NS <- G7_NS %>% select(vaccineStrain, compareStrain, Gene, Nmedian) %>% unique

## Spread the table!! ##
N_G7_Ndist_spread <- spread(G7_NS, key = Gene, value = Nmedian) 

## Build Correlation matrix - G7##
res2<-rcorr(as.matrix(N_G7_Ndist_spread[,c(3:10)]))
G7_NS_Cor <- flattenCorrMatrix(res2$r, res2$P)

colnames(G7_NS_Cor) <- c("Gene1", "Gene2", "Correlation", "p")

G7_NS_Cor <- G7_NS_Cor %>% 
  filter(Correlation != 1 & p < 0.05) 

G7_NS_Cor$GenePair <- paste(G7_NS_Cor$Gene1, "-", G7_NS_Cor$Gene2)


G7_NS_Cor %>% 
  ggplot(aes(x=reorder(GenePair, Correlation), y=Correlation)) + 
  geom_bar(stat = "identity") +
  geom_hline(yintercept= 0.3, linetype=4, color='gray', size=0.75) +
  geom_hline(yintercept= -0.3, linetype=4, color='gray', size=0.75) +
  geom_hline(yintercept= 0.6, linetype=3, color='orange2', size=1) +
  geom_hline(yintercept= -0.6, linetype=3, color='orange2', size=1) +
  ylim(-1,1) +
  xlab("Gene Pair") +
  coord_flip() +
  theme_bw()

## Filter Corr > 0.3 ##

G7_NS_Cor_03 <- G7_NS_Cor %>% 
  filter(Correlation >= 0.3) 

actors <- data.frame(name=c("PB2", "PB1", "PA", "HA", 
                            "NP", "NA", "M", "NS"))
relations <- data.frame(from= G7_NS_Cor_03$Gene1,
                        to= G7_NS_Cor_03$Gene2,
                        weight = G7_NS_Cor_03$Correlation)

g <- graph_from_data_frame(relations, directed=F, vertices=actors)

E(g)$width <- E(g)$weight*2.5^2

E(g)[E(g)$weight >= 0.3 & E(g)$weight < 0.595]$color <- "#ABB2B9" 
E(g)[E(g)$weight >= 0.595]$color <- "#F5B041" 

as <- authority_score(g, weights=E(g)$weight)$vector

plot(g, vertex.size=30, vertex.color=GeneCol, vertex.label.color="black", vertex.label.family = "sans",
     vertex.frame.width = 0.75)





########## for G8 ###########

## Let's use Median!! ##
G8_NS <- G8_NS %>%
  group_by(compareStrain, Gene) %>%
  mutate(Nmedian = mean(N2, na.rm = TRUE))

G8_NS <- G8_NS %>% select(vaccineStrain, compareStrain, Gene, Nmedian) %>% unique

## Spread the table!! ##
N_G8_Ndist_spread <- spread(G8_NS, key = Gene, value = Nmedian) 

## Build Correlation matrix - G8##
res2<-rcorr(as.matrix(N_G8_Ndist_spread[,c(3:10)]))
G8_NS_Cor <- flattenCorrMatrix(res2$r, res2$P)

colnames(G8_NS_Cor) <- c("Gene1", "Gene2", "Correlation", "p")

G8_NS_Cor <- G8_NS_Cor %>% 
  filter(Correlation != 1 & p < 0.05) 

G8_NS_Cor$GenePair <- paste(G8_NS_Cor$Gene1, "-", G8_NS_Cor$Gene2)


G8_NS_Cor %>% 
  ggplot(aes(x=reorder(GenePair, Correlation), y=Correlation)) + 
  geom_bar(stat = "identity") +
  geom_hline(yintercept= 0.3, linetype=4, color='gray', size=0.75) +
  geom_hline(yintercept= -0.3, linetype=4, color='gray', size=0.75) +
  geom_hline(yintercept= 0.6, linetype=3, color='orange2', size=1) +
  geom_hline(yintercept= -0.6, linetype=3, color='orange2', size=1) +
  ylim(-1,1) +
  xlab("Gene Pair") +
  coord_flip() +
  theme_bw()

## Filter Corr > 0.3 ##

G8_NS_Cor_03 <- G8_NS_Cor %>% 
  filter(Correlation >= 0.3) 

actors <- data.frame(name=c("PB2", "PB1", "PA", "HA", 
                            "NP", "NA", "M", "NS"))
relations <- data.frame(from= G8_NS_Cor_03$Gene1,
                        to= G8_NS_Cor_03$Gene2,
                        weight = G8_NS_Cor_03$Correlation)

g <- graph_from_data_frame(relations, directed=F, vertices=actors)

E(g)$width <- E(g)$weight*2.5^2

E(g)[E(g)$weight >= 0.3 & E(g)$weight < 0.595]$color <- "#ABB2B9" 
E(g)[E(g)$weight >= 0.595]$color <- "#F5B041" 

as <- authority_score(g, weights=E(g)$weight)$vector

plot(g, vertex.size=30, vertex.color=GeneCol, vertex.label.color="black", vertex.label.family = "sans",
     vertex.frame.width = 0.75)




########## for G9 ###########

## Let's use Median!! ##
G9_NS <- G9_NS %>%
  group_by(compareStrain, Gene) %>%
  mutate(Nmedian = mean(N2, na.rm = TRUE))

G9_NS <- G9_NS %>% select(vaccineStrain, compareStrain, Gene, Nmedian) %>% unique

## Spread the table!! ##
N_G9_Ndist_spread <- spread(G9_NS, key = Gene, value = Nmedian) 

## Build Correlation matrix - G9##
res2<-rcorr(as.matrix(N_G9_Ndist_spread[,c(3:10)]))
G9_NS_Cor <- flattenCorrMatrix(res2$r, res2$P)

colnames(G9_NS_Cor) <- c("Gene1", "Gene2", "Correlation", "p")

G9_NS_Cor <- G9_NS_Cor %>% 
  filter(Correlation != 1 & p < 0.05) 

G9_NS_Cor$GenePair <- paste(G9_NS_Cor$Gene1, "-", G9_NS_Cor$Gene2)


G9_NS_Cor %>% 
  ggplot(aes(x=reorder(GenePair, Correlation), y=Correlation)) + 
  geom_bar(stat = "identity") +
  geom_hline(yintercept= 0.3, linetype=4, color='gray', size=0.75) +
  geom_hline(yintercept= -0.3, linetype=4, color='gray', size=0.75) +
  geom_hline(yintercept= 0.6, linetype=3, color='orange2', size=1) +
  geom_hline(yintercept= -0.6, linetype=3, color='orange2', size=1) +
  ylim(-1,1) +
  xlab("Gene Pair") +
  coord_flip() +
  theme_bw()

## Filter Corr > 0.3 ##

G9_NS_Cor_03 <- G9_NS_Cor %>% 
  filter(Correlation >= 0.3) 

actors <- data.frame(name=c("PB2", "PB1", "PA", "HA", 
                            "NP", "NA", "M", "NS"))
relations <- data.frame(from= G9_NS_Cor_03$Gene1,
                        to= G9_NS_Cor_03$Gene2,
                        weight = G9_NS_Cor_03$Correlation)

g <- graph_from_data_frame(relations, directed=F, vertices=actors)

E(g)$width <- E(g)$weight*2.5^2

E(g)[E(g)$weight >= 0.3 & E(g)$weight < 0.595]$color <- "#ABB2B9" 
E(g)[E(g)$weight >= 0.595]$color <- "#F5B041" 

as <- authority_score(g, weights=E(g)$weight)$vector

plot(g, vertex.size=30, vertex.color=GeneCol, vertex.label.color="black", vertex.label.family = "sans",
     vertex.frame.width = 0.75)







############ Count Pairs + > 0.4 #########

## Filter Corr > 0.4 ##

G1_NS_Cor_04 <- G1_NS_Cor %>% 
  filter(Correlation >= 0.395) 
G2_NS_Cor_04 <- G2_NS_Cor %>% 
  filter(Correlation >= 0.395) 
G3_NS_Cor_04 <- G3_NS_Cor %>% 
  filter(Correlation >= 0.395) 
G4_NS_Cor_04 <- G5_NS_Cor %>% 
  filter(Correlation >= 0.395) 
G5_NS_Cor_04 <- G5_NS_Cor %>% 
  filter(Correlation >= 0.395) 
G6_NS_Cor_04 <- G6_NS_Cor %>% 
  filter(Correlation >= 0.395) 
G7_NS_Cor_04 <- G7_NS_Cor %>% 
  filter(Correlation >= 0.395) 
G8_NS_Cor_04 <- G8_NS_Cor %>% 
  filter(Correlation >= 0.395) 
G9_NS_Cor_04 <- G9_NS_Cor %>% 
  filter(Correlation >= 0.395) 

G1_NS_Cor_04$Group <- "G1"
G2_NS_Cor_04$Group <- "G2"
G3_NS_Cor_04$Group <- "G3"
G4_NS_Cor_04$Group <- "G4"
G5_NS_Cor_04$Group <- "G5"
G6_NS_Cor_04$Group <- "G6"
G7_NS_Cor_04$Group <- "G7"
G8_NS_Cor_04$Group <- "G8"
G9_NS_Cor_04$Group <- "G9"


AllGenepair_N_Cor_04 <- rbind(G1_NS_Cor_04, G2_NS_Cor_04,
                              G3_NS_Cor_04, G4_NS_Cor_04,
                              G5_NS_Cor_04, G6_NS_Cor_04,
                              G7_NS_Cor_04, G8_NS_Cor_04,
                              G9_NS_Cor_04)


AllGenepair_N_Cor_04 <- AllGenepair_N_Cor_04[complete.cases(AllGenepair_N_Cor_04), ]
AllGenepair_N_Cor_04[order(AllGenepair_N_Cor_04$GenePair),]


n <- AllGenepair_N_Cor_04 %>% 
  group_by(GenePair) %>% 
  summarise(n = n()) %>% 
  arrange(-n)


GenePairRank <- data.frame(n)
GenePairRank$Pair <- GenePairRank$GenePair

GenePairRank_2 <- GenePairRank %>% 
  separate(col = Pair, into = c("Gene1", "Gene2"), sep = " - ")


actors <- data.frame(name=c("PB2", "PB1", "PA", "HA", 
                            "NP", "NA", "M", "NS"))
relations <- data.frame(from= GenePairRank_2$Gene1,
                        to= GenePairRank_2$Gene2,
                        weight = GenePairRank_2$n)

g <- graph_from_data_frame(relations, directed=F, vertices=actors)

E(g)$width <- E(g)$weight^1.25

E(g)[E(g)$weight >= 5]$color <- "#873600" 
E(g)[E(g)$weight >= 3 & E(g)$weight < 5 ]$color <- "#F5B041" 
E(g)[E(g)$weight < 3]$color <- "#ABB2B9" 


w_rep      <- c(2, 4, 6)                 
lwd_legend <- w_rep^1.25                

col_legend <- c("#ABB2B9", "#F5B041", "#873600")
lab_legend <- c("< 3", "3–5", "≥ 5")


out_dir <- "3_Figures/Fig_2/Fig2E_GenPairCorr"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

svglite::svglite(
  filename = file.path(out_dir, "GenPairCorr.svg"),
  width = 350/25.4, height = 300/25.4
)

par(mar = c(4,4,2,2))
plot(
  g,
  vertex.size = 30,
  vertex.color = GeneCol,
  vertex.label.color = "black",
  vertex.label.family = "sans",
  vertex.frame.width = 0.75
)

legend("topleft",
       inset = 0.02,
       title = "Edge frequency",
       legend = lab_legend,
       col    = col_legend,
       lwd    = lwd_legend, 
       seg.len = 3,
       bty = "n")

dev.off()



################ For reference ################

############ Count Pairs + > 0.5 #########


## Filter Corr > 0.5 ##

G1_NS_Cor_05 <- G1_NS_Cor %>% 
  filter(Correlation >= 0.495) 
G2_NS_Cor_05 <- G2_NS_Cor %>% 
  filter(Correlation >= 0.495) 
G3_NS_Cor_05 <- G3_NS_Cor %>% 
  filter(Correlation >= 0.495) 
G4_NS_Cor_05 <- G5_NS_Cor %>% 
  filter(Correlation >= 0.495) 
G5_NS_Cor_05 <- G5_NS_Cor %>% 
  filter(Correlation >= 0.495) 
G6_NS_Cor_05 <- G6_NS_Cor %>% 
  filter(Correlation >= 0.495) 
G7_NS_Cor_05 <- G7_NS_Cor %>% 
  filter(Correlation >= 0.495) 
G8_NS_Cor_05 <- G8_NS_Cor %>% 
  filter(Correlation >= 0.495) 
G9_NS_Cor_05 <- G9_NS_Cor %>% 
  filter(Correlation >= 0.495) 

G1_NS_Cor_05$Group <- "G1"
G2_NS_Cor_05$Group <- "G2"
G3_NS_Cor_05$Group <- "G3"
G4_NS_Cor_05$Group <- "G4"
G5_NS_Cor_05$Group <- "G5"
G6_NS_Cor_05$Group <- "G6"
G7_NS_Cor_05$Group <- "G7"
G8_NS_Cor_05$Group <- "G8"
G9_NS_Cor_05$Group <- "G9"


AllGenepair_N_Cor_05 <- rbind(G1_NS_Cor_05, G2_NS_Cor_05,
                              G3_NS_Cor_05, G4_NS_Cor_05,
                              G5_NS_Cor_05, G6_NS_Cor_05,
                              G7_NS_Cor_05, G8_NS_Cor_05,
                              G9_NS_Cor_05)


AllGenepair_N_Cor_05 <- AllGenepair_N_Cor_05[complete.cases(AllGenepair_N_Cor_05), ]
AllGenepair_N_Cor_05[order(AllGenepair_N_Cor_05$GenePair),]


n <- AllGenepair_N_Cor_05 %>% 
  group_by(GenePair) %>% 
  summarise(n = n()) %>% 
  arrange(-n)


GenePairRank <- data.frame(n)
GenePairRank$Pair <- GenePairRank$GenePair

GenePairRank_2 <- GenePairRank %>% 
  separate(col = Pair, into = c("Gene1", "Gene2"), sep = " - ")

GenePairRank_2

actors <- data.frame(name=c("PB2", "PB1", "PA", "HA", 
                            "NP", "NA", "M", "NS"))
relations <- data.frame(from= GenePairRank_2$Gene1,
                        to= GenePairRank_2$Gene2,
                        weight = GenePairRank_2$n)

g <- graph_from_data_frame(relations, directed=F, vertices=actors)

E(g)$width <- E(g)$weight^1.5

E(g)[E(g)$weight == 2]$color <- "orange" 
E(g)[E(g)$weight == 1]$color <- "gray80" 

plot(g, vertex.size=30, vertex.color=GeneCol, vertex.label.color="black", vertex.label.family = "sans",
     vertex.frame.width = 0.5)


############ Count Pairs + > 0.3 #########

## Filter Corr > 0.3 ##

G1_NS_Cor_03 <- G1_NS_Cor %>% 
  filter(Correlation >= 0.295) 
G2_NS_Cor_03 <- G2_NS_Cor %>% 
  filter(Correlation >= 0.295) 
G3_NS_Cor_03 <- G3_NS_Cor %>% 
  filter(Correlation >= 0.295) 
G4_NS_Cor_03 <- G5_NS_Cor %>% 
  filter(Correlation >= 0.295) 
G5_NS_Cor_03 <- G5_NS_Cor %>% 
  filter(Correlation >= 0.295) 
G6_NS_Cor_03 <- G6_NS_Cor %>% 
  filter(Correlation >= 0.295) 
G7_NS_Cor_03 <- G7_NS_Cor %>% 
  filter(Correlation >= 0.295) 
G8_NS_Cor_03 <- G8_NS_Cor %>% 
  filter(Correlation >= 0.295) 
G9_NS_Cor_03 <- G9_NS_Cor %>% 
  filter(Correlation >= 0.295) 

G1_NS_Cor_03$Group <- "G1"
G2_NS_Cor_03$Group <- "G2"
G3_NS_Cor_03$Group <- "G3"
G4_NS_Cor_03$Group <- "G4"
G5_NS_Cor_03$Group <- "G5"
G6_NS_Cor_03$Group <- "G6"
G7_NS_Cor_03$Group <- "G7"
G8_NS_Cor_03$Group <- "G8"
G9_NS_Cor_03$Group <- "G9"


AllGenepair_N_Cor_03 <- rbind(G1_NS_Cor_03, G2_NS_Cor_03,
                              G3_NS_Cor_03, G4_NS_Cor_03,
                              G5_NS_Cor_03, G6_NS_Cor_03,
                              G7_NS_Cor_03, G8_NS_Cor_03,
                              G9_NS_Cor_03)


AllGenepair_N_Cor_03 <- AllGenepair_N_Cor_03[complete.cases(AllGenepair_N_Cor_03), ]
AllGenepair_N_Cor_03[order(AllGenepair_N_Cor_03$GenePair),]


n <- AllGenepair_N_Cor_03 %>% 
  group_by(GenePair) %>% 
  summarise(n = n()) %>% 
  arrange(-n)


GenePairRank <- data.frame(n)
GenePairRank$Pair <- GenePairRank$GenePair

GenePairRank_2 <- GenePairRank %>% 
  separate(col = Pair, into = c("Gene1", "Gene2"), sep = " - ")

GenePairRank_2

actors <- data.frame(name=c("PB2", "PB1", "PA", "HA", 
                            "NP", "NA", "M", "NS"))
relations <- data.frame(from= GenePairRank_2$Gene1,
                        to= GenePairRank_2$Gene2,
                        weight = GenePairRank_2$n)

g <- graph_from_data_frame(relations, directed=F, vertices=actors)

E(g)$width <- E(g)$weight^1.25

E(g)[E(g)$weight >= 5]$color <- "#873600" 
E(g)[E(g)$weight >= 3 & E(g)$weight < 5 ]$color <- "#F5B041" 
E(g)[E(g)$weight < 3]$color <- "#ABB2B9" 


plot(g, vertex.size=30, vertex.color=GeneCol, vertex.label.color="black", vertex.label.family = "sans",
     vertex.frame.width = 0.75)



##### by p and cor values



############ All N Mutation, Line by correlation #########
AllG_NS <- rbind(G1_NS, G2_NS, G3_NS, G4_NS, G5_NS, 
                 G6_NS, G7_NS, G8_NS, G9_NS)

## Build Correlation matrix - G2##
N_All_Ndist_spread <- spread(AllG_NS, key = Gene, value = Nmedian) 

## Correlation matrix - r and p
res2 <- rcorr(as.matrix(N_All_Ndist_spread[,c(3:10)]))
All_N_Cor <- flattenCorrMatrix(res2$r, res2$P)

colnames(All_N_Cor) <- c("Gene1", "Gene2", "Correlation", "p")

All_N_Cor <- All_N_Cor %>% 
  filter(Correlation != 1 & p < 0.05) 

All_N_Cor$GenePair <- paste(All_N_Cor$Gene1, "-", All_N_Cor$Gene2)

All_N_Cor %>% 
  ggplot(aes(x=reorder(GenePair, Correlation), y=Correlation)) + 
  geom_bar(stat = "identity") +
  geom_hline(yintercept= 0.3, linetype=4, color='gray', size=0.75) +
  geom_hline(yintercept= -0.3, linetype=4, color='gray', size=0.75) +
  geom_hline(yintercept= 0.6, linetype=3, color='orange2', size=1) +
  geom_hline(yintercept= -0.6, linetype=3, color='orange2', size=1) +
  ylim(-1,1) +
  xlab("Gene Pair") +
  coord_flip() +
  theme_bw()

## Filter Corr > 0.3 ##
actors <- data.frame(name=c("PB2", "PB1", "PA", "HA", 
                            "NP", "NA", "M", "NS"))
relations <- data.frame(from= All_N_Cor$Gene1,
                        to= All_N_Cor$Gene2,
                        weight = All_N_Cor$Correlation)

g <- graph_from_data_frame(relations, directed=F, vertices=actors)

E(g)$width <- E(g)$weight*1.08^30

E(g)[E(g)$weight < 0.45]$color <- "#ABB2B9" 
E(g)[E(g)$weight >= 0.45 & E(g)$weight < 0.55]$color <- "#F5B041" 
E(g)[E(g)$weight >= 0.55]$color <- "#873600" 


plot(g, vertex.size=30, vertex.color=GeneCol, vertex.label.color="black", vertex.label.family = "sans",
     vertex.frame.width = 0.75)

## Filter p <0.05 ##

G1_NS_CorP <- G1_NS_Cor %>% 
  filter(p < 0.05) 
G2_NS_CorP <- G2_NS_Cor %>% 
  filter(p < 0.05) 
G3_NS_CorP <- G3_NS_Cor %>% 
  filter(p < 0.05) 
G4_NS_CorP <- G5_NS_Cor %>% 
  filter(p < 0.05) 
G5_NS_CorP <- G5_NS_Cor %>% 
  filter(p < 0.05) 
G6_NS_CorP <- G6_NS_Cor %>% 
  filter(p < 0.05) 
G7_NS_CorP <- G7_NS_Cor %>% 
  filter(p < 0.05) 
G8_NS_CorP <- G8_NS_Cor %>% 
  filter(p < 0.05) 
G9_NS_CorP <- G9_NS_Cor %>% 
  filter(p < 0.05) 

G1_NS_CorP$Group <- "G1"
G2_NS_CorP$Group <- "G2"
G3_NS_CorP$Group <- "G3"
G4_NS_CorP$Group <- "G4"
G5_NS_CorP$Group <- "G5"
G6_NS_CorP$Group <- "G6"
G7_NS_CorP$Group <- "G7"
G8_NS_CorP$Group <- "G8"
G9_NS_CorP$Group <- "G9"

AllGenepair_N_CorP <- rbind(G1_NS_CorP, G2_NS_CorP,
                              G3_NS_CorP, G4_NS_CorP,
                              G5_NS_CorP, G6_NS_CorP,
                              G7_NS_CorP, G8_NS_CorP,
                              G9_NS_CorP)


AllGenepair_N_CorP <- AllGenepair_N_CorP[complete.cases(AllGenepair_N_CorP), ]
AllGenepair_N_CorP[order(AllGenepair_N_CorP$GenePair),]


n <- AllGenepair_N_CorP %>% 
  group_by(GenePair) %>% 
  summarise(n = n()) %>% 
  arrange(-n)


GenePairRank <- data.frame(n)
GenePairRank$Pair <- GenePairRank$GenePair

GenePairRank_2 <- GenePairRank %>% 
  separate(col = Pair, into = c("Gene1", "Gene2"), sep = " - ")

GenePairRank_2

actors <- data.frame(name=c("PB2", "PB1", "PA", "HA", 
                            "NP", "NA", "M", "NS"))
relations <- data.frame(from= GenePairRank_2$Gene1,
                        to= GenePairRank_2$Gene2,
                        weight = GenePairRank_2$n)

g <- graph_from_data_frame(relations, directed=F, vertices=actors)

E(g)$width <- E(g)$weight^1.25

E(g)[E(g)$weight >= 6]$color <- "#873600" 
E(g)[E(g)$weight >= 4 & E(g)$weight < 6 ]$color <- "#F5B041" 
E(g)[E(g)$weight < 4]$color <- "#ABB2B9" 

plot(g, vertex.size=30, vertex.color=GeneCol, vertex.label.color="black", vertex.label.family = "sans",
     vertex.frame.width = 0.75)

