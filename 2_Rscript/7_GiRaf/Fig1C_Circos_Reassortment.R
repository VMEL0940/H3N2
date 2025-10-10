## Package ##
library(tidyverse)
library(circlize)
library(ggraph)
library(igraph)
library(RColorBrewer)
library(edgebundleR)
library(grDevices)


## Set working directory ##
setwd("~/H3N2/")

### import "GiRaf" data 
GiRaf_Summary <- read.csv("1_Data/7_GIRAF/GiRaf_summary.csv", header = TRUE, na.string = "") # 

##### 1. Tile Plot #####

GiRaf_summary_GeneIdx <- GiRaf_Summary %>% 
  group_by(GenePair) %>% 
  dplyr::summarize(n = n()) %>% 
  arrange(-n)

## split by gene
table2 <- separate(data = GiRaf_summary_GeneIdx, col = GenePair, into = c("Gene1", "Gene2"), sep = "-")

table2$Gene1 <- factor(table2$Gene1, levels = c("PB2","PB1","PA","HA",
                                                "NP", "NA", "M", "NS"))

table2$Gene2 <- factor(table2$Gene2, levels = c("PB2","PB1","PA","HA",
                                                "NP", "NA", "M", "NS"))

table2 <- table2 %>% 
  arrange(Gene2, Gene1)

table2$Gene1[26] <- "NS"
table2$Gene2[26] <- "NP"

table2$Gene1[22] <- "NA"
table2$Gene2[22] <- "HA"

table2$Gene1[25] <- "NS"
table2$Gene2[25] <- "HA"

table2$Gene1[24] <- "M"
table2$Gene2[24] <- "HA"

table2$Gene1[19] <- "NP"
table2$Gene2[19] <- "HA"

table2$Gene1[27] <- "NS"
table2$Gene2[27] <- "NA"

table2$Gene1[28] <- "NS"
table2$Gene2[28] <- "M"

table2 <- table2 %>% 
  arrange(Gene2, Gene1)

table2 <- table2 %>% 
  mutate(n_s = scale(n))



##### 2. Chord gram -- by ggraph #####
actors <- data.frame(name=c("PB2", "PB1", "PA", "HA", 
                            "NP", "NA", "M", "NS"))

grid.col = c("PB2" = "#C780E8", "PB1" = "#9D94FF", "PA" = "#59ADF6", "HA" = "#08CAD1",
             "NP" = "#42D6A4", "NA" = "#F8F38D", "M" = "#FFB480", "NS" = "#FF6961")

relations <- data.frame(from= table2$Gene1,
                        to= table2$Gene2,
                        weight = table2$n)

g <- graph_from_data_frame(relations, directed=F, vertices=actors)

V(g)$degree <- strength(g)

GeneCol <- c("#C780E8", "#9D94FF","#59ADF6","#08CAD1","#42D6A4","#F8F38D","#FFB480","#FF6961")



p <- ggraph(g, layout = 'linear', circular = TRUE) + 
    geom_edge_arc(aes(color = weight, edge_width = weight)) +  # 두께 매핑 추가
    geom_node_text(aes(label = name, angle = node_angle(0, 0)), hjust = -0.9) +
    geom_node_point(aes(fill = name, size = degree), shape = 21) +
    theme_graph() +
    scale_edge_colour_distiller(palette = 1, trans = "reverse", name = "Frequency") + 
    scale_edge_width(range = c(0.5, 2), guide = "none") +     # 두께 스케일 설정
    scale_fill_manual(values = grid.col, guide="none") +
    coord_fixed(xlim = c(-1.4, 1.4), ylim = c(-1.4, 1.4)) +
    scale_size(range = c(2,12), guide="none")


ggsave("Circos_Reassort.svg", plot = p, device = "svg", 
       path = "3_Figures/Fig_1/Fig1C_Reassort",
       width = 350, height = 300, units = "mm", dpi = 720)

################ For reference ################
Chord_Table <- table2[,c(1:3)]

chordDiagram(Chord_Table)
chordDiagram(Chord_Table, order = c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"))

chordDiagram(Chord_Table, order = c("NS","M","NA","NP","HA","PA", "PB1", "PB2"), 
             grid.col = grid.col, transparency = 0.2)
