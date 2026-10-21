setwd("E:/博士2023-2025/课题/Mettl3/10.C.sulf/全基因组测序/customer_backup/function_annotation/01.general_database_annotation/03.kegg_annotation")
data=read.csv("./kegg_pathway_anno.csv")

data=data[order(data$Gene_number,decreasing = F),]
data$KEGG_pathway.level.2=factor(data$KEGG_pathway.level.2,levels = unique(data$KEGG_pathway.level.2))
data=data[order(data$KEGG_pathway.level.2),]
data$X.pathway=factor(data$X.pathway,levels = unique(data$X.pathway))
data=data[data$Gene_number>10,]

library(ggplot2)
library(RColorBrewer)
ggplot(data,aes(y=X.pathway,x=Gene_number,))+
  geom_col(aes(fill=KEGG_pathway.level.2))+
  scale_fill_manual(values = colorRampPalette(brewer.pal(12,"Set3"))(19)#c("#E1C092","#A3C9DC","#C6B0D1")
                    )+
  labs(x="Count",y="KEGG pathway")+
  geom_text(
    data = data,
    aes(y = X.pathway, x = 0, label = X.pathway),
    # angle = 90,          # 文字垂直旋转
    hjust = 0,            # 对齐方式，使文字从底部开始
    vjust = 0.5,
    size = 1.5            # 调整字号以适应图大小
  ) +
  scale_x_continuous(expand = c(0,1))+
  theme(text = element_text(size = 6),
        axis.text.y = element_blank(),axis.line = element_line(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        # panel.border = element_rect(fill = NULL),
        legend.position="inside",
        legend.position.inside = c(0.8,0.2),
        legend.title = element_blank(),
        legend.key.size = unit(0.2,"cm"))
ggsave("E:/博士2023-2025/课题/Mettl3/10.C.sulf/全基因组测序/KEGG.pdf",width = 8/2.54,height = 16/2.54)
