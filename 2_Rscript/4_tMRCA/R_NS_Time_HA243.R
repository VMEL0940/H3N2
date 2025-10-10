### Influ H3N2 GenDist Proj ###
## Package ##
library(tidyverse)

## Set working directory ##
setwd("/01_GenDisFlu/Evaluation/")

## Data import ##
AllG_NS <- read.csv("HA242_NSdist.csv", header = T, na.strings = "")

AllG_NS$Vaccine_code <- factor(AllG_NS$Group,levels = c("G01","G02","G03","G04","G05","G06","G07","G08","G09",
                                                        "G10","G11","G12","G13","G14","G15","G16","G17"),
                                             labels = c("1.HK68","2.EN72", "3.PC73", "4.VI75","5.TE77","6.PH82","7.LE86","8.SI87","9.SH87",
                                                        "10.GU89","11.BE89","12.BE92","13.SD93","14.JO94","15.WU95","16.SY97","17.MW99"))

jpeg(filename = "/01_GenDisFlu/Fig/NS_Time/HA243_Eval_N_Time.jpeg", width = 300, height = 150, units = "mm", res = 520)

AllG_NS %>%
  ggplot() +
  geom_jitter(aes(x = Year, y = Nmedian, fill = Vaccine_code), size = 2.5, shape = 21, color = "gray20",
              height = 0, width = 0.3) +
  geom_vline(xintercept = c(1972.5,1974.5,1976.5,1978.5,1983.5,1987.5,
                            1988.5,1990.5,1991.5,1993.5,1994.5,
                            1995.5,1996.5,1998.5,2000.5), color = "red", lty = 2, lwd = 0.4) +
  ylab("Non-synonymous RNSC of the HA gene") +
  scale_fill_discrete(name = "Vaccine Strains") +
  scale_y_continuous(breaks = seq(0, 45, by = 5)) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9))

dev.off()
  

jpeg(filename = "/01_GenDisFlu/Fig/NS_Time/HA243_Eval_S_Time.jpeg", width = 300, height = 150, units = "mm", res = 520)

AllG_NS %>%
  ggplot() +
  geom_jitter(aes(x = Year, y = Smedian, fill = Vaccine_code), size = 2.5, shape = 21, color = "gray20",
              height = 0, width = 0.3) +
  geom_vline(xintercept = c(1972.5,1974.5,1976.5,1978.5,1983.5,1987.5,
                            1988.5,1990.5,1991.5,1993.5,1994.5,
                            1995.5,1996.5,1998.5,2000.5), color = "red", lty = 2, lwd = 0.4) +
  ylab("Synonymous RNSC of the HA gene") +
  scale_fill_discrete(name = "Vaccine Strains") +
  scale_x_continuous(breaks = seq(1965, 2005, by = 5)) +
  scale_y_continuous(breaks = seq(0, 45, by = 5)) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9))

dev.off()


#### By Domiance -- N value #### 

AllG_NS$Dominance_F <- factor(AllG_NS$Dominance, levels = c(0, 1), labels = c("Extinct", "Dominant"))

jpeg(filename = "/01_GenDisFlu/Fig/NS_Time/HA243_Eval_N_Dom_Time.jpeg", width = 300, height = 150, units = "mm", res = 520)

AllG_NS %>%
  ggplot() +
  geom_jitter(aes(x = Year, y = Nmedian, fill = Dominance_F), size = 2.5, shape = 21, color = "gray20",
                height = 0, width = 0.3) +
  geom_vline(xintercept = c(1972.5,1974.5,1976.5,1978.5,1983.5,1987.5,
                            1988.5,1990.5,1991.5,1993.5,1994.5,
                            1995.5,1996.5,1998.5,2000.5), color = "red", lty = 2, lwd = 0.4) +
  ylab("Non-synonymous RNSC of the HA gene") +
  scale_fill_discrete(name = "Vaccine Strains") +
  scale_y_continuous(breaks = seq(0, 45, by = 5)) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9))


dev.off()

  


jpeg(filename = "/01_GenDisFlu/Fig/NS_Time/HA243_Eval_S_Dom_Time.jpeg", width = 300, height = 150, units = "mm", res = 520)

AllG_NS %>%
  ggplot() +
  geom_jitter(aes(x = Year, y = Smedian, fill = Dominance_F), size = 2.5, shape = 21, color = "gray20",
              height = 0, width = 0.3) +
  geom_vline(xintercept = c(1972.5,1974.5,1976.5,1978.5,1983.5,1987.5,
                            1988.5,1990.5,1991.5,1993.5,1994.5,
                            1995.5,1996.5,1998.5,2000.5), color = "red", lty = 2, lwd = 0.4) +
  ylab("Synonymous RNSC of the HA gene") +
  scale_fill_discrete(name = "Vaccine Strains") +
  scale_y_continuous(breaks = seq(0, 45, by = 5)) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9))


dev.off()

