setwd("path")

library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
data=read.csv("./tax_cal.csv",sep = "\t")
p=ggscatter(data, 
            x = "s__Clostridium_sulfidigenes", y = "Fecal_Calprotectin",
            color = "black",
            add = "reg.line", 
            conf.int = TRUE,
            add.params = list( color = "red",fill = "grey",fill = "grey"),
            cor.coef = T,
            size = 1,
            cor.coef.size = 2.5,
            cor.coeff.args = list(),
            cor.method = "spearman",
            repel=T)+
  labs(y="Fecal Calprotectin",x="Abundance of C.sulfidigenes")+
  theme_bw(base_size = 8)+
  theme(
    panel.border = element_blank(),
    axis.line = element_line(colour = "black",size = 0.3),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5))

p

ggsave("Cs_calprotectin_correlation.pdf",
       p, width = 4/2.53, height = 4/2.53)



p1=ggscatter(data, 
            x = "BSH", y = "Fecal_Calprotectin",
            color = "black",
            add = "reg.line", 
            conf.int = TRUE,
            add.params = list( color = "red",fill = "grey",fill = "grey"),
            cor.coef = T,
            size = 1,
            cor.coef.size = 2.5,
            cor.coeff.args = list(),
            cor.method = "spearman",
            repel=T)+
  labs(y="Fecal Calprotectin",x="Abundance of BSH (mus)")+
  theme_bw(base_size = 8)+
  theme(
    panel.border = element_blank(),
    axis.line = element_line(colour = "black",size = 0.3),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5))

p1

ggsave("BSH_calprotectin_correlation.pdf",
       p1, width = 4/2.53, height = 4/2.53)
