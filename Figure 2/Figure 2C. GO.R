setwd("path")
data=read.csv("./GO_sec_stat.csv")
data=data[-1,]
data=data[order(data$Genome,decreasing = T),]
data$X.GO_classify1=factor(data$X.GO_classify1,levels = unique(data$X.GO_classify1))
data=data[order(data$X.GO_classify1),]
data$GO_classify2=factor(data$GO_classify2,levels = unique(data$GO_classify2))

library(ggplot2)

ggplot(data,aes(x=GO_classify2,y=Genome,))+
  geom_col(aes(fill=X.GO_classify1))+
  scale_fill_manual(values = c("#E1C092","#A3C9DC","#C6B0D1"))+
  labs(x="",y="Count")+
  geom_text(
    data = data,
    aes(x = GO_classify2, y = 0, label = GO_classify2),
    angle = 90,  
    hjust = 0,   
    vjust = 0.5,
    size = 1.5 
  ) +scale_y_continuous(expand = c(0,10))+
  theme(text = element_text(size = 6),
        axis.text.x = element_blank(),axis.line = element_line(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position="inside",
        legend.position.inside = c(0.8,0.9),
        legend.title = element_blank(),
        legend.key.size = unit(0.2,"cm"))
ggsave("./GO.pdf",width = 13/2.54,height = 2.5/2.54)
