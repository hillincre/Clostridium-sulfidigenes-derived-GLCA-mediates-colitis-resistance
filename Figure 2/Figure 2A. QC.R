setwd("path")
data=read.table("./GCvsDepth.txt",header = T)

library(ggplot2)

p1=  ggplot(data, aes(x = GC., y = Med_DEPTH)) +
    geom_hex(bins = 30) +
    scale_fill_gradientn(colors = c("darkblue", "cyan", "yellow", "red", "darkred")) +
    labs(x = "GC content[%]", y = "Median Depth") +
    theme(text = element_text(size = 6),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_rect(fill = NULL),
          legend.key.size = unit(0.2,"cm"),
          legend.position = "right")

p2=ggplot(data,aes(x=GC.))+
  geom_density(fill="#89bbd9")+
  labs(x = "GC content[%]")+
  theme(text = element_text(size = 6),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NULL),
        legend.key.size = unit(0.5,"cm"))

p2

p3=ggplot(data,aes(x=Med_DEPTH))+
  geom_density(fill="#89bbd9")+
  labs(x = "Median Depth")+
  theme(text = element_text(size = 6),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NULL),
        legend.key.size = unit(0.5,"cm"))

p3

library(patchwork)
p1/p2/p3

ggsave("./QC.pdf",width = 6/2.54,height = 8/2.54)
