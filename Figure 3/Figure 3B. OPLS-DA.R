rm(list=ls())
setwd("path")
library(ropls)
library(ggplot2)
data=read.csv("./04_AllMet.csv",row.names = 4)
data=data[,-c(1:3,grep("QC",names(data)))]
group=data.frame(group=gsub("\\_.*","",names(data)),row.names = names(data))
group$group=factor(group$group,levels = c("CON","CS"))
group=group[order(group$group),,drop=F]
data=data[,match(rownames(group),names(data))]
data=data.frame(t(data))

df1_oplsda <- opls(data, group$group, predI = 12, orthoI = 0)

data <- as.data.frame(df1_oplsda@scoreMN)
o1 <- df1_oplsda@orthoScoreMN[,1]
data$o1 <- o1
data$group = group$group
data$samples = rownames(data)
x_lab <- df1_oplsda@modelDF[1, "R2X"] * 100
y_lab <- df1_oplsda@modelDF[2, "R2X"] * 100
col=c("#4e62ab", "#d6404e")
p1 <- ggplot(data,aes(x=p1,y=p2,color=group))+
  theme_bw()+
  geom_point(size=3)+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed",color="grey")+
  geom_hline(yintercept = 0,lty="dashed",color="grey")+
  labs(x=paste0("P1 (",x_lab,"%)"),
       y=paste0("P2(",y_lab,"%)"))+
  stat_ellipse(data=data,
               geom = "path",level = 0.95,
               linetype = 2,size=0.5,
               alpha=0.2,
               show.legend = T)+
  scale_color_manual(values = col)+
  scale_fill_manual(values = col)+
  theme(axis.title.x=element_text(size=8),
        axis.title.y=element_text(size=8,angle=90),
        axis.text.y=element_text(size=8),
        axis.text.x=element_text(size=8),
        panel.grid=element_blank(),
        plot.title = element_text(size = 10, hjust = 0.5),
        legend.title = element_blank()) +
  ggtitle("PLS-DA")

pdf("./OPLS-DA.pdf",width = 6/2.53,height = 4.5/2.53)
p1
dev.off()