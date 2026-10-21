setwd("path")

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

obj=readRDS("./anchored_data_processed.rds")
immu=subset(obj,!celltype%in%c("Endothelial cells","Epithelial cells","Fibroblasts","Dendritic cells"))

DotPlot(immu,features = c("Gprc5a","Gpr157","Cxcr5"),scale = F)+
  coord_flip()+
  scale_color_gradientn(colors=c("#1E3A78","lightgrey","#B22"))+
  theme_bw(base_size = 6)+
  labs(x="",y="")+
  theme(legend.key.size = unit(0.3,"cm"),
        legend.position = "top",
        legend.title.position = "top",
        axis.text.x = element_text(angle = -45,hjust =0,vjust = 0.5))
ggsave("./GPCRs_expression.pdf",width = 4.8/2.53,height = 4/2.53)