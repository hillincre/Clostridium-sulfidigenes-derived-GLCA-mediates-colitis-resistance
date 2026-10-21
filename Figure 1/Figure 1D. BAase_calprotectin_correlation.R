library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)

setwd("path")

# ---------- 1. ----------
dat <- read.delim("secondary_bile_acid_enzyme_per_sample.tsv",
                  stringsAsFactors = FALSE, check.names = FALSE)

enzymes <- c(
  "7a-HSDH",
  "cbh / BSH"
  # "4.2.1.106 (bile-acid 7alpha-dehydratase)",
  # "1.1.1.176 (12alpha-hydroxysteroid dehydrogenase)"
)
names(dat)[4]= "7a-HSDH"
names(dat)[5]= "cbh / BSH"

dat$Fecal_Calprotectin <- suppressWarnings(as.numeric(dat$Fecal_Calprotectin))
# dat=dat[,1:5]
# ---------- 2. ----------
dat_long <- dat %>%
  select(Sample, Group, Fecal_Calprotectin, all_of(enzymes)) %>%
  pivot_longer(cols = all_of(enzymes), names_to = "Enzyme", values_to = "Abundance") %>%
  filter(!is.na(Fecal_Calprotectin))
dat_long=dat_long[!(dat_long$Enzyme=="cbh / BSH"&
                      dat_long$Abundance>25 &
                      dat_long$Abundance<50&
                      dat_long$Fecal_Calprotectin<50),]
write.csv(dat_long,"./BSH_Cal_correlation.csv")
p=ggscatter(dat_long[dat_long$Enzyme=="cbh / BSH",], 
            x = "Abundance", y = "Fecal_Calprotectin",
            color = "black",
            # fill = "lightgray",
            add = "reg.line", 
            conf.int = TRUE,
            add.params = list( color = "red",fill = "grey",fill = "grey"),
            cor.coef = T,
            size = 1,
            cor.coef.size = 2.5,
            cor.coeff.args = list(),
            cor.method = "spearman",
            repel=T)+
  labs(y="Fecal Calprotectin",x="Abundance of BSH (Homo)")+
  # facet_wrap(~Enzyme,ncol=1,scales = "free_y",axes = "all",strip.position = "right")+
  theme_bw(base_size = 8)+
  theme(
    panel.border = element_blank(),
    axis.line = element_line(colour = "black",size = 0.3),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5))+ylim(c(-50,750))

p

ggsave("BSH_calprotectin_correlation.pdf",
       p, width = 4/2.53, height = 4/2.53)