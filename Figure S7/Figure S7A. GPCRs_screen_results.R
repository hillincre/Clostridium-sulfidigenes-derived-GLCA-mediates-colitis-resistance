setwd("path")
data=xlsx::read.xlsx("GPCRs_screen.xlsx",1,header = T)
data=rbind(data,0)
data$rank=c(nrow(data):1)
library(ggplot2)
library(ggrepel)
p=ggplot(data,aes(x=rank,y=mean))+
  geom_point(data = data[-179,],shape=21,aes(fill=mean,size=mean),color="white")+
  scale_fill_gradient(low = "grey",high = "#4393C4")+
  scale_size_continuous(range = c(1,8))+
  geom_text_repel(data = head(data,20),aes(label=SYMBOL),max.overlaps = 100,size=2)+
  scale_y_continuous(expand = c(0,0.01))+
  labs(x="Ranks",y="Normalized value",title = "Expression of GPCRs")+
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        text = element_text(size = 8),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,size = 8))
p

ggsave("./GPCRs_Rank.pdf",units = "cm",width = 10,height = 9)

library(ggpie)
library(ggplot2)
library(patchwork)
temp=data.frame(Expresssion=c("Yes","No"),count=c(179,113))
p1 <- ggdonut(temp,
              group_key = "Expresssion",
              donut.label = T,
              donut.label.color = "black",
              donut.label.size = 2,
              label_size = 2,
              count_type = "count",
              border_size = 1,
              fill_color = c("#4393C4","grey"),
              r0 = 1.5,
              r1 = 3) +
  theme(text = element_text(size = 8),
        legend.key.size = unit(0.3,"cm"),
        legend.position = "right",
        legend.title = element_blank()) 
pdf("./genecount_pie.pdf",width = 7/2.53,height = 5/2.53)
p1
dev.off()