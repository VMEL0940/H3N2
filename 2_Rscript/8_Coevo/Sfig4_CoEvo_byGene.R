### Influ H3N2 CoEvo Position ###
## Package ##
library(tidyverse)

## Set working directory ##
setwd("~/H3N2/")

## Data import ##
HA_M <- read.csv("1_Data/8_CoEvo/Epistasis Analysis/Epistasis Analysis/Datamonkey Results/M_datamonkey_table.csv", header = T, na.strings = "")
HA_NA <- read.csv("1_Data/8_CoEvo/Epistasis Analysis/Epistasis Analysis/Datamonkey Results/NA_datamonkey_table.csv", header = T, na.strings = "")
HA_PA <- read.csv("1_Data/8_CoEvo/Epistasis Analysis/Epistasis Analysis/Datamonkey Results/PA_datamonkey_table.csv", header = T, na.strings = "")
HA_NS <- read.csv("1_Data/8_CoEvo/Epistasis Analysis/Epistasis Analysis/Datamonkey Results/NS_datamonkey_table.csv", header = T, na.strings = "")
HA_NP <- read.csv("1_Data/8_CoEvo/Epistasis Analysis/Epistasis Analysis/Datamonkey Results/NP_datamonkey_table.csv", header = T, na.strings = "")
HA_PB1 <- read.csv("1_Data/8_CoEvo/Epistasis Analysis/Epistasis Analysis/Datamonkey Results/PB1_datamonkey_table.csv", header = T, na.strings = "")
HA_PB2 <- read.csv("1_Data/8_CoEvo/Epistasis Analysis/Epistasis Analysis/Datamonkey Results/PB2_datamonkey_table.csv", header = T, na.strings = "")

## Data Cleaning
colnames(HA_M)[1] <- "HA"
colnames(HA_M)[2] <- "Position"
HA_M$Gene <- "M"

colnames(HA_NA)[1] <- "HA"
colnames(HA_NA)[2] <- "Position"
HA_NA$Gene <- "NA"

colnames(HA_PA)[1] <- "HA"
colnames(HA_PA)[2] <- "Position"
HA_PA$Gene <- "PA"

colnames(HA_NS)[1] <- "HA"
colnames(HA_NS)[2] <- "Position"
HA_NS$Gene <- "NS"

colnames(HA_NP)[1] <- "HA"
colnames(HA_NP)[2] <- "Position"
HA_NP$Gene <- "NP"

colnames(HA_PB1)[1] <- "HA"
colnames(HA_PB1)[2] <- "Position"
HA_PB1$Gene <- "PB1"

colnames(HA_PB2)[1] <- "HA"
colnames(HA_PB2)[2] <- "Position"
HA_PB2$Gene <- "PB2"

## Binding all data 
CoE_Gene <- rbind(HA_M, HA_NA, HA_NP, HA_PA, HA_PB1, HA_PB2, HA_NS)

## Subsetting non-gene gene interaction region and align position
CoE_Gene_sb <- CoE_Gene %>% 
  filter(HA < 567 & Position > 566)

CoE_Gene_sb$Position <- CoE_Gene_sb$Position - 566  ## Subtract HA seq Number 


#### ADVANCED!!!! RBS Information!!! !! ! ! ! #####

sum(CoE_Gene_sb$HA %in% c(149, 151, 161, 171, 172, 174, 175, 205, 209)) ##0
sum(CoE_Gene_sb$HA %in% c(114, 150, 151, 152, 153, 169, 171, 199, 206, 210, 241, 242, 244))  ##2
sum(CoE_Gene_sb$HA %in% c(112, 113, 115, 116, 147, 148, 149, 154, 155, 156, 157, 158, 159, 
                      160, 161, 162, 163, 164, 167, 170, 172, 200, 201, 202, 203, 204, 
                      205, 207, 208, 209, 211, 234, 235, 236, 237, 238, 239, 240, 243, 
                      245, 246, 266, 267, 268)) ##20
sum(CoE_Gene_sb$HA %in% c(138,140,142,146,147,148,149,151,153,154,156,157,158,159,160,161,162,166,168,184,
                          144,145,171,172,173,174,175,176,179,180,181,202,203,204,205,206,208,209,210,212,213,214,215,
                          60,61,62,63,64,66,67,69,70,289,291,292,294,295,296,310,313,315,316,320,321,323,324,325,326,327,328,
                          112,118,119,133,137,183,186,187,188,189,190,191,192,193,195,198,217,219,223,224,225,228,229,230,231,232,233,
                          234,235,238,239,241,242,243,244,245,246,249,254,256,258,260,262,263,264,
                          73,75,78,79,83,91,94,96,97,98,99,102,103,104,107,108,110,125,276,277,278,281))  ##82

CoE_Gene_sb$HA
CoE_Gene_sb$RBD15A <- ifelse(CoE_Gene_sb$HA %in% c(112, 113, 115, 116, 147, 148, 149, 154, 155, 156, 157, 158, 159, 
                                            160, 161, 162, 163, 164, 167, 170, 172, 200, 201, 202, 203, 204, 
                                            205, 207, 208, 209, 211, 234, 235, 236, 237, 238, 239, 240, 243, 
                                            245, 246, 266, 267, 268), "RBD-15A", 
                             ifelse(CoE_Gene_sb$HA %in% c(114, 150, 151, 152, 153, 169, 171, 199, 206, 210, 241, 242, 244), "RBD", "NO"))

CoE_Gene_sb$AgBABCDE <- ifelse(CoE_Gene_sb$HA %in% c(138,140,142,146,	147,148,149,151,153,154,156,157,158,159,160,161,162,166,168,184,
                                                     144,145,171,172,173,174,175,176,179,180,181,202,203,204,205,206,208,209,210,212,213,214,215,
                                                     60,61,62,63,64,66,67,69,70,289,291,292,294,295,296,310,313,315,316,320,321,323,324,325,326,327,328,
                                                     112,118,119,133,137,183,186,187,188,189,190,191,192,193,195,198,217,219,223,224,225,228,229,230,231,232,233,
                                                     234,235,238,239,241,242,243,244,245,246,249,254,256,258,260,262,263,264,
                                                     73,75,78,79,83,91,94,96,97,98,99,102,103,104,107,108,110,125,276,277,278,281), "AgBABCDE", "NO")
                                                    


## Figures ##
p <- CoE_Gene_sb %>% 
    filter(Gene == "M") %>%
    ggplot() +
    geom_segment(aes(x=1, y = 1, xend = 1, yend = 566), color = '#08CAD1', size = 5) +
    geom_segment(aes(x=2, y = 1, xend = 2, yend = 349), color = '#FFB480', size = 5) +  ## M == 349
    scale_colour_manual(name = "Receptor binding sites",values = c("black", "red4", "red1")) +
    geom_point(aes(x= 1, y = HA, color = RBD15A), size=2.5) + 
    geom_point(aes(x= 2, y = Position, color = RBD15A), size=2.5) +
    geom_text(aes(x= 1, y = HA, label= HA, color = RBD15A), nudge_x=-0.05, size=2) +
    geom_text(aes(x= 2, y = Position, label= Position, color = RBD15A), nudge_x=0.05, size=2) +
    geom_segment(aes(x=1, y= HA, xend=2, yend= Position, color = RBD15A), size = 0.2) + 
    annotate('text', x=1, y=-20, label='HA', size=5, color="black") +
    annotate('text', x=2, y=-20, label='M', size=5, color="black") +
    scale_y_reverse() +
    theme_void() +
    theme(legend.position = "none") 

ggsave("Sfig_4_H_M_Coevolution.svg", plot = p, device = "svg",
       path = "5_SFigure/Sfig4_Coevo",
       width = 350, height = 700, units = "mm", dpi = 720)


p <- CoE_Gene_sb %>% 
    filter(Gene == "NA") %>%
    ggplot() +
    geom_segment(aes(x=1, y = 1, xend = 1, yend = 566), color = '#08CAD1', size = 5) +
    geom_segment(aes(x=2, y = 1, xend = 2, yend = 469), color = '#F8F38D', size = 5) + ## NA == 469
    scale_colour_manual(name = "Receptor binding sites",values = c("black", "red1")) +
    geom_point(aes(x= 1, y = HA, color = RBD15A), size=2.5) + 
    geom_point(aes(x= 2, y = Position, color = RBD15A), size=2.5) +
    geom_text(aes(x= 1, y = HA, label= HA, color = RBD15A), nudge_x=-0.05, size=2) +
    geom_text(aes(x= 2, y = Position, label= Position, color = RBD15A), nudge_x=0.05, size=2) +
    geom_segment(aes(x=1, y= HA, xend=2, yend= Position, color = RBD15A), size = 0.2) + 
    annotate('text', x=1, y=-20, label='HA', size=5, color="black") +
    annotate('text', x=2, y=-20, label='NA', size=5, color="black") +
    scale_y_reverse() +
    theme_void() +
    theme(legend.position = "none") 


ggsave("Sfig_4_B_HA_NA.svg", plot = p, device = "svg",
       path = "5_SFigure/Sfig4_Coevo",
       width = 350, height = 700, units = "mm", dpi = 720)


p <- CoE_Gene_sb %>% 
    filter(Gene == "NP") %>%
    ggplot() +
    geom_segment(aes(x=1, y = 1, xend = 1, yend = 566), color = '#08CAD1', size = 5) +
    geom_segment(aes(x=2, y = 1, xend = 2, yend = 498), color = '#42D6A4', size = 5) +  ## NP == 498
    scale_colour_manual(name = "Receptor binding sites",values = c("black", "red1")) +
    geom_point(aes(x= 1, y = HA, color = RBD15A), size=2.5) + 
    geom_point(aes(x= 2, y = Position, color = RBD15A), size=2.5) +
    geom_text(aes(x= 1, y = HA, label= HA, color = RBD15A), nudge_x=-0.05, size=2) +
    geom_text(aes(x= 2, y = Position, label= Position, color = RBD15A), nudge_x=0.05, size=2) +
    geom_segment(aes(x=1, y= HA, xend=2, yend= Position, color = RBD15A), size = 0.2) + 
    annotate('text', x=1, y=-20, label='HA', size=5, color="black") +
    annotate('text', x=2, y=-20, label='NP', size=5, color="black") +
    scale_y_reverse() +
    theme_void() +
    theme(legend.position = "none") 


ggsave("Sfig_4_G_NP_Coevolution.svg", plot = p, device = "svg",
       path = "5_SFigure/Sfig4_Coevo",
       width = 350, height = 700, units = "mm", dpi = 720)

p <- CoE_Gene_sb %>% 
    filter(Gene == "NS") %>%
    ggplot() +
    geom_segment(aes(x=1, y = 1, xend = 1, yend = 566), color = '#08CAD1', size = 5) +
    geom_segment(aes(x=2, y = 1, xend = 2, yend = 351), color = '#FF6961', size = 5) +  ## NS == 351
    scale_colour_manual(name = "Receptor binding sites",values = c("black", "red4", "red1")) +
    geom_point(aes(x= 1, y = HA, color = RBD15A), size=2.5) + 
    geom_point(aes(x= 2, y = Position, color = RBD15A), size=2.5) +
    geom_text(aes(x= 1, y = HA, label= HA, color = RBD15A), nudge_x=-0.05, size=2) +
    geom_text(aes(x= 2, y = Position, label= Position, color = RBD15A), nudge_x=0.05, size=2) +
    geom_segment(aes(x=1, y= HA, xend=2, yend= Position, color = RBD15A), size = 0.2) + 
    annotate('text', x=1, y=-20, label='HA', size=5, color="black") +
    annotate('text', x=2, y=-20, label='NS', size=5, color="black") +
    scale_y_reverse() +
    theme_void() +
    theme(legend.position = "none") 


ggsave("Sfig_4_C_HA_NS.svg", plot = p, device = "svg",
       path = "5_SFigure/Sfig4_Coevo",
       width = 350, height = 700, units = "mm", dpi = 720)

p <- CoE_Gene_sb %>% 
    filter(Gene == "PA") %>%
    ggplot() +
    geom_segment(aes(x=1, y = 1, xend = 1, yend = 566), color = '#08CAD1', size = 5) +
    geom_segment(aes(x=2, y = 1, xend = 2, yend = 716), color = '#59ADF6', size = 5) +  ## PA == 716
    scale_colour_manual(name = "Receptor binding sites",values = c("black", "red1")) +
    geom_point(aes(x= 1, y = HA, color = RBD15A), size=2.5) + 
    geom_point(aes(x= 2, y = Position, color = RBD15A), size=2.5) +
    geom_text(aes(x= 1, y = HA, label= HA, color = RBD15A), nudge_x=-0.05, size=2) +
    geom_text(aes(x= 2, y = Position, label= Position, color = RBD15A), nudge_x=0.05, size=2) +
    geom_segment(aes(x=1, y= HA, xend=2, yend= Position, color = RBD15A), size = 0.2) + 
    annotate('text', x=1, y=-20, label='HA', size=5, color="black") +
    annotate('text', x=2, y=-20, label='PA', size=5, color="black") +
    scale_y_reverse() +
    theme_void() +
    theme(legend.position = "none") 

ggsave("Sfig_4_F_PA_Coevolution.svg", plot = p, device = "svg",
       path = "5_SFigure/Sfig4_Coevo",
       width = 350, height = 700, units = "mm", dpi = 720)


p <- CoE_Gene_sb %>% 
    filter(Gene == "PB1") %>%
    ggplot() +
    geom_segment(aes(x=1, y = 1, xend = 1, yend = 566), color = '#08CAD1', size = 5) +
    geom_segment(aes(x=2, y = 1, xend = 2, yend = 757), color = '#9D94FF', size = 5) +  ## Pb1 == 757
    scale_colour_manual(name = "Receptor binding sites",values = c("black", "red1")) +
    geom_point(aes(x= 1, y = HA, color = RBD15A), size=2.5) + 
    geom_point(aes(x= 2, y = Position, color = RBD15A), size=2.5) +
    geom_text(aes(x= 1, y = HA, label= HA, color = RBD15A), nudge_x=-0.05, size=2) +
    geom_text(aes(x= 2, y = Position, label= Position, color = RBD15A), nudge_x=0.05, size=2) +
    geom_segment(aes(x=1, y= HA, xend=2, yend= Position, color = RBD15A), size = 0.2) + 
    annotate('text', x=1, y=-20, label='HA', size=5, color="black") +
    annotate('text', x=2, y=-20, label='PB1', size=5, color="black") +
    scale_y_reverse() +
    theme_void() +
    theme(legend.position = "none") 

ggsave("Sfig_4_E_PB1_Coevolution.svg", plot = p, device = "svg",
       path = "5_SFigure/Sfig4_Coevo",
       width = 350, height = 700, units = "mm", dpi = 720)


p <- CoE_Gene_sb %>% 
    filter(Gene == "PB2") %>%
    ggplot() +
    geom_segment(aes(x=1, y = 1, xend = 1, yend = 566), color = '#08CAD1', size = 5) +
    geom_segment(aes(x=2, y = 1, xend = 2, yend = 759), color = '#C780e8', size = 5) +  ## PB2 == 759
    scale_colour_manual(name = "Receptor binding sites",values = c("black", "red")) +
    geom_point(aes(x= 1, y = HA, color = RBD15A), size=2.5) + 
    geom_point(aes(x= 2, y = Position, color = RBD15A), size=2.5) +
    geom_text(aes(x= 1, y = HA, label= HA, color = RBD15A), nudge_x=-0.05, size=2) +
    geom_text(aes(x= 2, y = Position, label= Position, color = RBD15A), nudge_x=0.05, size=2) +
    geom_segment(aes(x=1, y= HA, xend=2, yend= Position, color = RBD15A), size = 0.2) + 
    annotate('text', x=1, y=-20, label='HA', size=5, color="black") +
    annotate('text', x=2, y=-20, label='PB2', size=5, color="black") +
    scale_y_reverse() +
    theme_void() +
    theme(legend.position = "none") 


ggsave("Sfig_4_D_PB2_Coevolution.svg", plot = p, device = "svg",
       path = "5_SFigure/Sfig4_Coevo",
       width = 350, height = 700, units = "mm", dpi = 720)



## Data import ##
HA_HA <- read.csv("1_Data/8_CoEvo/Epistasis Analysis/Epistasis Analysis/Datamonkey Results/HA_datamonkey_table.csv", header = T, na.strings = "")

## Data Cleaning
colnames(HA_HA)[1] <- "Position1"
colnames(HA_HA)[2] <- "Position2"

HA_HA$RBD15A_1 <- ifelse(HA_HA$Position1 %in% c(112, 113, 115, 116, 147, 148, 149, 154, 155, 156, 157, 158, 159, 
                                               160, 161, 162, 163, 164, 167, 170, 172, 200, 201, 202, 203, 204, 
                                               205, 207, 208, 209, 211, 234, 235, 236, 237, 238, 239, 240, 243, 
                                               245, 246, 266, 267, 268), "RBD-15A", 
                         ifelse(HA_HA$Position1 %in% c(114, 150, 151, 152, 153, 169, 171, 199, 206, 210, 241, 242, 244), "RBD", "NO"))
HA_HA$RBD15A_2 <- ifelse(HA_HA$Position2 %in% c(112, 113, 115, 116, 147, 148, 149, 154, 155, 156, 157, 158, 159, 
                                               160, 161, 162, 163, 164, 167, 170, 172, 200, 201, 202, 203, 204, 
                                               205, 207, 208, 209, 211, 234, 235, 236, 237, 238, 239, 240, 243, 
                                               245, 246, 266, 267, 268), "RBD-15A", 
                         ifelse(HA_HA$Position2 %in% c(114, 150, 151, 152, 153, 169, 171, 199, 206, 210, 241, 242, 244), "RBD", "NO"))

HA_HA$RBD15A_3 <- ifelse(HA_HA$RBD15A_1 == "RBD-15A" | HA_HA$RBD15A_2 == "RBD-15A", "RBD-15A" , 
                         ifelse(HA_HA$RBD15A_1 == "RBD" | HA_HA$RBD15A_2 == "RBD", "RBD" , "NO"))

summary(HA_HA$RBD15A_1)


p <- HA_HA %>% 
    ggplot() +
    geom_segment(aes(x=1, y = 1, xend = 1, yend = 566), color = '#08CAD1', size = 5) +
    scale_colour_manual(name = "Receptor binding site", values = c("black", "red4", "red1")) +
    geom_point(aes(x= 1, y = Position1, color = RBD15A_1), size=2.5) + 
    geom_point(aes(x= 1, y = Position2, color = RBD15A_2), size=2.5) + 
    geom_curve(aes(x = 1, y = Position1, xend = 1, yend = Position2, color = RBD15A_3), curvature = 0.5, size = 0.2) +
    geom_text(aes(x= 1, y = Position1, label= Position1, color = RBD15A_1), nudge_x=0.005, size=2) +
    geom_text(aes(x= 1, y = Position2, label= Position2, color = RBD15A_2), nudge_x=0.005, size=2) +
    annotate('text', x=1, y=-20, label='HA', size=7, color="black") +
    xlim(0.95, 1.05) +
    scale_y_reverse() +
    theme_void()

ggsave("Sfig_4_A_HA_HA.svg", plot = p, device = "svg",
       path = "5_SFigure/Sfig4_Coevo",
       width = 300, height = 450, units = "mm", dpi = 720)

