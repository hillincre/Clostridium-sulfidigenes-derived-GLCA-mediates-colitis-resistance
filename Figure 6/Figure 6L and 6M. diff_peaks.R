setwd("path")

data=read.csv("./DIFF_peaks_all.csv")
data$Regulate="no"
data$Regulate[data$Fold>0 & data$p.value<0.05]="up"
data$Regulate[data$Fold<0 & data$p.value<0.05]="down"

library(ggplot2)
library(ggrepel)
data <- data[order(data$Regulate == "no", decreasing = TRUE), ]
label=data[data$X%in%c(15574,15575),,drop=F]
label$X[label$X==15574]="Rorc promoter-TSS"
label$X[label$X==15575]="Rorc intron\n(intron 3 of 9)"

up <- sum(data$Regulate == "up"); down <- sum(data$Regulate == "down")
ggplot(data,aes(log10(Conc_ctrl),log10(Conc_treat)))+
  geom_point(aes(color=Regulate),shape=20,size=0.05)+
  theme_classic()+
  scale_color_manual(values = c("down"="#b4c6e2","up"= "#e4c8ca"))+
  xlim(0.4,1)+ylim(0.5,1)+
  labs(x="Log10(PBS mean)",
       y="Log10(GLCA mean)",
       title = paste0("Con vs GLCA\nDown: ", down, "  Up: ", up))+
  geom_point(data = label,aes(log10(Conc_ctrl),log10(Conc_treat)),shape=21,size=0.05,color="black")+
  geom_text_repel(
    data = label,
    aes(label = X),
    size = 1,
    hjust = 0,
    vjust=-1,
    color = "black",force_pull = 3,
    show.legend = FALSE
  )+
  theme(title=element_text(size=6),plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NULL),
        text=element_text(size = 6),
        axis.line = element_blank(),
        axis.title.x = element_text(size = 6),
        axis.title.y = element_text(size = 6),
        legend.title = element_blank(),
        legend.position = "none")
ggsave("./diff_peaks_volcanol.pdf",width = 3.5/2.53,height = 4/2.53)





########KEGG enrichment#########
res=read.csv("./ATAC_UP_KEGG.csv")
res$BgCount <- as.numeric(sub("/.*", "", res$BgRatio))
res$RichFactor <- res$Count / res$BgCount
res$GeneTotal <- as.numeric(sub(".*/", "", res$GeneRatio))
res$BgTotal   <- as.numeric(sub(".*/", "", res$BgRatio))

res$FoldEnrichment <- (res$Count / res$GeneTotal) / (res$BgCount / res$BgTotal)

head(res[, c("ID", "Description", "Count", "BgCount", 
             "RichFactor", "FoldEnrichment", "p.adjust", "qvalue")])

res=res[order(res$RichFactor,decreasing = F),]
res$Description=factor(res$Description,levels = res$Description)
p= ggplot(res, aes(x = RichFactor, y = Description)) +
  geom_col(fill = "#e4c8ca",
  )+
  labs(x="RichFactor",y="KEGG pathway")+
  geom_text(
    data = res,
    aes(y = Description, x = 0, label = Description),
    hjust = 0,
    vjust = 0.5,
    size = 2.5
  ) +
  scale_x_continuous(expand = c(0,0))+
  theme(text = element_text(size = 8),
        axis.text.y = element_blank(),axis.line = element_line(),
        axis.line.y = element_blank(),axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position="inside",
        legend.position.inside = c(0.8,0.2),
        legend.title = element_blank(),
        legend.key.size = unit(0.2,"cm"))
p
ggsave("./kegg_bar.pdf",plot = p,width = 4/2.53,height = 5.5/2.53)