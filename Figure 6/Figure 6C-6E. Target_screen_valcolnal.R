setwd("path")
data=read.csv("./DE_proteins_ALL.csv")

data$Regulate="no"
data$Regulate[data$FC.GLCA.PBS.>1 & data$P_value.GLCA.PBS.<0.05]="up"
data$Regulate[data$FC.GLCA.PBS.<1 & data$P_value.GLCA.PBS.<0.05]="down"

library(ggplot2)
library(ggrepel)
label=data[c(grep("Guanine nucleotide-binding protein",data$Protein.Name),
             grep("Arrb",data$Protein.Name)),,drop=F]
label=label[label$Regulate!="no",]
data <- data[order(data$Regulate == "no", decreasing = TRUE), ]
data=data[data$GLCA!=0.1 & data$PBS!=0.1,]
up <- sum(data$Regulate == "up"); down <- sum(data$Regulate == "down")
ggplot(data,aes(log10(PBS),log10(GLCA)))+
  geom_point(aes(color=Regulate),shape=20,size=1)+
  theme_classic()+
  scale_color_manual(values = c("down"="#b4c6e2","up"= "#e4c8ca"))+
  labs(x="Log10(PBS mean)",y="Log10(GLCA mean)",
       title = paste0("Con vs GLCA\nDown: ", down, "  Up: ", up)
  )+
  geom_point(data = label,aes(log10(PBS),log10(GLCA)),color="black",
             shape=21,size=1)+
  geom_text_repel(
    data = label,
    aes(label = Symbol),force = 1,
    size = 1,
    hjust = 0,
    vjust=-1,
    color = "black",
    show.legend = FALSE
  ) +
  theme(title=element_text(size=6),plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NULL),
        text=element_text(size = 6),
        axis.line = element_blank(),
        axis.title.x = element_text(size = 6),
        axis.title.y = element_text(size = 6),
        legend.title = element_blank(),
        legend.position = "none")
ggsave("./volcanol.pdf",width = 3/2.53,height =3.8/2.53)

data=read.csv("./DE_proteins_ALL.csv")

data$Regulate="no"
data$Regulate[data$FC.GLCA.PBS.>1 & data$P_value.GLCA.PBS.<0.05]="up"
data$Regulate[data$FC.GLCA.PBS.<1 & data$P_value.GLCA.PBS.<0.05]="down"

library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(patchwork)
genes=data$Symbol[data$Regulate=="up"]
trans=bitr(genes,"SYMBOL","ENTREZID",OrgDb = org.Mm.eg.db)
kegg=enrichKEGG(trans$ENTREZID,organism = "mmu",
                pvalueCutoff = 1,
                pAdjustMethod = "BH",
                qvalueCutoff = 1)

res=kegg@result
res=res[res$subcategory=="Signal transduction",]
res=na.omit(res)
res=res[order(res$RichFactor,decreasing = F),]
res$Description=factor(res$Description,levels = res$Description)
write.csv(res,"./KEGG_up_protein.csv")

df=tail(res,10)
p= ggplot(df, aes(x = RichFactor, y = Description)) +
  geom_col(fill = "#e4c8ca")+
  labs(x="RichFactor",y="KEGG pathway")+
  geom_text(
    data = df,
    aes(y = Description, x = 0, label = Description),
    hjust = 0,
    vjust = 0.5,
    size = 1.5
  ) +
  scale_x_continuous(expand = c(0,0))+
  theme(text = element_text(size = 6),
        axis.text.y = element_blank(),axis.line = element_line(),
        axis.line.y = element_blank(),axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position="inside",
        legend.position.inside = c(0.8,0.2),
        legend.title = element_blank(),
        legend.key.size = unit(0.2,"cm"))
p
ggsave("./kegg_bar.pdf",plot = p,width = 3/2.53,height = 3/2.53)

pro1=res$geneID[grep("cAMP",res$Description)]
pro2=res$geneID[grep("cGMP",res$Description)]

pro1 <- unlist(strsplit(pro1, "/", fixed = TRUE))
pro2 <- unlist(strsplit(pro2, "/", fixed = TRUE))

pro1=trans[trans$ENTREZID%in%pro1,]
pro2=trans[trans$ENTREZID%in%pro2,]

df1=data[data$Symbol %in% pro1$SYMBOL,]
df2=data[data$Symbol %in% pro2$SYMBOL,]
df1=df1[order(df1$FC.GLCA.PBS.,decreasing = F),]
df1$Symbol=factor(df1$Symbol,levels = df1$Symbol)

df2=df2[order(df2$FC.GLCA.PBS.,decreasing = F),]
df2$Symbol=factor(df2$Symbol,levels = df2$Symbol)

p1= ggplot(df1, aes(x = FC.GLCA.PBS., y = Symbol)) +
  geom_col(fill = "#e4c8ca")+
  labs(x="Fold Change",y="")+
  geom_text(
    data = df1,
    aes(y = Symbol, x = 0, label = Symbol),
    hjust = 0,
    vjust = 0.5,
    size = 1.5
  ) +labs(title="cGMP-PKG signaling pathway")+
  scale_x_continuous(expand = c(0,0.1))+
  theme(text = element_text(size = 6),
        axis.text.y = element_blank(),axis.line = element_line(),
        panel.background = element_blank(),
        panel.grid = element_blank(),panel.border = element_rect(fill = NULL),
        legend.position="inside",
        legend.position.inside = c(0.8,0.2),
        legend.title = element_blank(),
        legend.key.size = unit(0.2,"cm"))
p1

p2= ggplot(df2, aes(x = FC.GLCA.PBS., y = Symbol)) +
  geom_col(fill = "#e4c8ca")+
  labs(x="Fold Change",y="")+
  geom_text(
    data = df2,
    aes(y = Symbol, x = 0, label = Symbol),
    hjust = 0,
    vjust = 0.5,
    size = 1.5
  ) +labs(title ="cAMP signaling pathway",x="")+
  scale_x_continuous(expand = c(0,0.1))+
  theme(text = element_text(size = 6),
        axis.text.y = element_blank(),axis.line = element_line(),
        panel.background = element_blank(),
        panel.grid = element_blank(),panel.border = element_rect(fill = NULL),
        legend.position="inside",
        legend.position.inside = c(0.8,0.2),
        legend.title = element_blank(),
        legend.key.size = unit(0.2,"cm"))
p2

p2/p1
ggsave("./int_path_pro_abun.pdf",p2/p1,width = 3/2.53,height = 5/2.53)