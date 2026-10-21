setwd("path")

data=read.table("./GSE225392_merged_gene_counts.txt",header = T)
data=data[,c(2,grep("DMSO",colnames(data)))]
data=data[!duplicated(data$gene_name),]
rownames(data)=data$gene_name
data=data[,-1]

group=data.frame(group=sub(".*(KODMSO|WTDMSO)$", "\\1", colnames(data)),
                 row.names = colnames(data))

library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(AnnotationDbi)
library(ggplot2)
library(patchwork)

cts <- as.matrix(data)
coldata <- data.frame(group = factor(group$group))

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ group)
dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)

res <- results(dds, contrast = c("group", "KODMSO", "WTDMSO"))
res <- lfcShrink(dds, contrast = c("group", "KODMSO", "WTDMSO"), res = res, type = "ashr")

de_markers=as.data.frame(res)
de_markers$regulation="no"
de_markers$regulation[de_markers$log2FoldChange > 0.3785116&de_markers$padj<0.05]="up"
de_markers$regulation[de_markers$log2FoldChange < -0.3785116&de_markers$padj<0.05]="down"

de_markers=read.csv("./DEGs.csv",row.names = 1)
de_markers$regulation="no"

up_count=table(de_markers$regulation)["up"]
down_count=table(de_markers$regulation)["down"]

de_markers$gene=rownames(de_markers)

p1=ggplot(de_markers,aes(log2FoldChange,-log2(pvalue)))+
  geom_point(aes(color=regulation),alpha=0.5,size=0.5)+
  labs(x="log2(Fold Change)",
       y="-log2(padj)",
       title = paste0("DEGs induced by Gprc5a KO \n","Down: ",down_count," - ","Up: ",up_count)
  )+
  scale_color_manual(values = c("#66c2a5", "grey", "#f46d43"))+
  theme_bw(base_size = 8)+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(linewidth = 0.5))
p1

ggsave("./DEGs_volcano.png",p1,units = "cm",width = 5.5,height = 6)

library(tidyr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggsignif)
norm <- counts(dds, normalized = TRUE)
exp=norm[grep("GNAI3",rownames(norm)),,drop=F]
exp=as.data.frame(exp)
exp=pivot_longer(exp,cols = names(exp),names_to = "sample")
exp$group=sub(".*(KODMSO|WTDMSO)$", "\\1", exp$sample)
exp$group=gsub("DMSO","",exp$group)
mytheme= 
  theme(text=element_text(size = 8),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black",linewidth = 0.5),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",legend.margin = margin(0,0,0,0),
        plot.title = element_text(hjust = 0.5)
        
  )
exp$group=factor(exp$group,levels = c("WT","KO"))
p=ggplot(exp, aes(x = group, y = value, color = group)) +
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(width = 0.5, outlier.color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1) +
  scale_color_manual(values = c("WT" = "#4e62ab", "KO" = "#d6404e")) +
  geom_signif(
    comparisons = list(c("WT","KO")),
    map_signif_level = function(p) sprintf("p = %.2g", p),
    test = "t.test",
    tip_length = 0,
    textsize = 2.5
  ) +
  labs(x="",y="Normolized value")+
  mytheme +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))
p
ggsave("GNAI3_expression.pdf", p, width = 4/2.53, height = 4/2.53)

library(ggplot2)
library(scales)

heat_color = c("#0f2c7a","#345aa6","#5c88d1","#8bace9","#c8d8f6","white",
               "#f6c8bc","#e98b7a","#d15c4b","#a63439","#7a0f29")

int = c("ARRB1","ARRB2",
        "GNAQ","GNA11","GNA14","GNA15",
        "GNAI1","GNAI2","GNAI3","GNAO1",
        "GNA12","GNA13",
        "LATS1","YAP1",
        "RORC","IL17D","IL22RA1")

heatplot = function(int, data){
  df = data[data$gene %in% int, ]
  df$gene = factor(df$gene, levels = int)
  
  df$log2FoldChange = as.numeric(df$log2FoldChange)
  df$pvalue         = as.numeric(df$pvalue)
  
  df$label = ifelse(
    is.na(df$pvalue), "",
    ifelse(df$pvalue < 0.0001, "***",
           ifelse(df$pvalue < 0.01,   "**",
                  ifelse(df$pvalue < 0.05,   "*", "")))
  )
  
  val = max(abs(df$log2FoldChange), na.rm = TRUE)
  if (!is.finite(val) || val == 0) val = 0.1
  
  df$fill_norm = abs(df$log2FoldChange) / val
  df$text_color = ifelse(df$fill_norm > 0.45, "white", "black")
  
  df$y = 1
  
  p1 = ggplot(df, aes(x = gene, y = y)) +
    geom_tile(aes(fill = log2FoldChange), color = "gray90", linewidth = 0.2) +
    scale_fill_gradientn(
      colors = heat_color, na.value = "white",
      limits = c(-val, val),
      breaks = c(-round(val, 1), 0, round(val, 1)),
      values = c(0, 0.5, 1),
      name   = "log2FC"
    ) +
    geom_text(aes(y = 1, label = label, color = text_color),
              size = 3, show.legend = FALSE) +
    scale_color_identity() +
    scale_y_continuous(
      limits = c(0.5, 1.5), expand = c(0, 0)
    ) +
    labs(y = "", x = "") +
    theme_bw(base_size = 8) +
    theme(
      text          = element_text(size = 8),
      axis.text      = element_text(size = 8),
      axis.text.x    = element_text(size = 8, angle = 45, hjust = 1),
      axis.text.y    = element_blank(),
      axis.ticks     = element_blank(),
      axis.title     = element_text(size = 8),
      legend.text    = element_text(size = 6),
      legend.title   = element_text(size = 8),
      legend.title.position = "top",
      legend.position = "top",
      legend.key.size = unit(0.3, "cm"),
      plot.title     = element_text(size = 8, hjust = 0.5),
      panel.grid     = element_blank(),
      panel.border   = element_blank()
    )
  return(p1)
}

G = heatplot(int, de_markers)
G

ggsave(plot = G,"./int_gene_exp.pdf",width = 10/2.53,height = 4/2.53)

enrich=function(gene_list){
  gene_entrez <- mapIds(org.Hs.eg.db,
                        keys=gene_list, keytype="SYMBOL",
                        column="ENTREZID", multiVals="first")
  gene_entrez <- na.omit(gene_entrez)
  
  kk2 <- enrichKEGG(
    gene = gene_entrez,
    organism = "hsa",
    pvalueCutoff = 0.05
  )
  kegg_result=as.data.frame(kk2@result)
  return(kegg_result)
}

deg1 <- de_markers %>%
  filter(regulation == "up")
deg2 <- de_markers %>%
  filter(regulation == "down")

up=enrich(rownames(deg1))
down=enrich(rownames(deg2))

up$regulation="up"
down$regulation="down"

down$RichFactor=-down$RichFactor
kegg=rbind(up,down)
write.csv(kegg,"./kegg_enrich.csv")

kegg=kegg[kegg$category=="Environmental Information Processing",]
kegg=na.omit(kegg)
kegg=kegg[order(kegg$RichFactor,decreasing = F),]

kegg2=rbind(head(kegg,10),tail(kegg,10))
kegg2=kegg2[!duplicated(kegg2$Description),]

kegg2$Description=factor(kegg2$Description,levels = kegg2$Description)
kegg2$regulation=factor(kegg2$regulation,levels = c("down","up"))

p2= ggplot(kegg2, aes(x = RichFactor, y = Description, fill = regulation)) +
  geom_bar(stat = "identity", width = 0.8, color = "white") +
  scale_fill_manual(values = c("down" = "#66c2a5", "up" = "#f46d43"),
                    name = NULL, labels = c("down",  "up")) +
  scale_x_continuous( limits = c(-0.2,0.2),breaks = seq(-0.2,0.2, 0.1),
                      labels = abs(seq(-0.2,0.2, 0.1))) +
  labs(
    title = "",
    y = ""
  ) +
  geom_text(data = kegg2[kegg2$regulation=="down",],aes(x=0.01,label=Description),hjust=0,size.unit = "pt",size=8)+
  geom_text(data = kegg2[kegg2$regulation=="up",],aes(x=-0.01,label=Description),hjust=1,size.unit = "pt",size=8)+
  theme_minimal() +
  coord_cartesian(expand = T,clip = "off") +
  theme(
    text=element_text(size = 8),
    axis.text.y = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.5),
    axis.line.x = element_line(linewidth = 0.5),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.key.size = unit(0.5, "cm"),
    plot.margin = margin(r=10,l = 20,t=0)
  )
p2

ggsave("./kegg_bar.pdf",plot = p2,width = 10/2.53,height = 11.6/2.53)

library(ggplot2)
library(patchwork)
library(clusterProfiler)

GKO=read.csv("./DEGs.csv")
ILC3s=read.csv("./DEGs_in_ILC3s_m3KO.csv")
ILC3s$X=toupper(ILC3s$X)

UC_ILC3=read.csv("./DEGs_in_ILC3s_UC.csv")
integrate=data.frame(genes=intersect(GKO$X,ILC3s$X))
integrate=left_join(integrate,GKO,by = join_by(genes==X))
integrate=left_join(integrate,ILC3s,by = join_by(genes==X))
integrate$regulation="no"
integrate$regulation[integrate$regulation.x != "no" | integrate$regulation.y!= "no"]="significant"
integrate$regulation[integrate$regulation.x == "down"&integrate$regulation.y=="up"]="overlap"

label=integrate[integrate$genes=="ZNRF3",]
p1=ggplot(integrate,aes(log2FoldChange,avg_log2FC))+
  geom_point(aes(color=regulation),alpha=1,size=0.5)+
  labs(x="log2(Fold Change)\n(Gprc5a KO / Control)",
       y="log2(Fold Change)\n(mettl3DIE / mettl3f/f)",
  )+
  scale_color_manual(values = c("significant"="#999999",
                                "no"="grey", "overlap"="#f46d43"))+
  geom_point(data=label,aes(log2FoldChange,avg_log2FC),shape=21,size=0.5,color="black")+
  geom_text(data = label,aes(label=gene),vjust=1,size = 5,size.unit = "pt")+
  theme_bw(base_size = 8)+
  theme(legend.position = "top",legend.title = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),margins = margin(r=20,l=10),
        panel.border = element_rect(linewidth = 0.5))
p1

ggsave("./DEGs_volcano.pdf",p1,units = "cm",width = 5.5,height = 5)

enrich=function(gene_list){
  gene_entrez <- mapIds(org.Hs.eg.db,
                        keys=gene_list, keytype="SYMBOL",
                        column="ENTREZID", multiVals="first")
  gene_entrez <- na.omit(gene_entrez)
  
  kk2 <- enrichKEGG(
    gene = gene_entrez,
    organism = "hsa",
    pvalueCutoff = 0.05
  )
  kegg_result=as.data.frame(kk2@result)
  return(kegg_result)
}

kegg=enrich(lap)

go <- enrichGO(
  gene = lap,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05
)
go_result <- as.data.frame(go)

library(ggplot2)
library(scales)
data=read.csv("./DEGs.csv")
int1=c(
  "GNAQ","GNA11","GNA14","GNA15",
  "GNAI1","GNAI2","GNAI3","GNAO1",
  "GNA12","GNA13")
int2=c(
  "NFKB1","NFKB2","RELA","RELB","REL"
)
int3=c( "RORC","TCF7","RORA","ID2","NFIL3","TOX")

heat_color=c("#0f2c7a","#345aa6","#5c88d1","#8bace9","#c8d8f6","white","#f6c8bc","#e98b7a","#d15c4b","#a63439","#7a0f29")
heatplot=function(int){
  df=data[data$gene %in% int,]
  df$gene=factor(df$gene,levels = int)
  df$NAME="expression"
  df$label=ifelse(df$pvalue<0.05,ifelse(df$pvalue>0.01,"*",ifelse(df$pvalue>0.0001,"**","***")),"")
  val=max(abs(df$log2FoldChange))
  p1=ggplot(df,aes(x=NAME,y=gene))+
    geom_tile(aes(fill = log2FoldChange), color = "gray90", linewidth = 0.2) + 
    scale_fill_gradientn(colors = heat_color,na.value = "white",
                         limits = c(-val, val),breaks=c(-round(val,digits = 1),0, round(val,digits = 1)),
                         values = rescale(c(-val, 0, val)) )+
    labs(y="",x="")+
    geom_text(aes(label=label),color="white")+
    theme_bw(base_size = 8) +
    theme(
      text = element_text(size = 8),
      axis.text = element_text(size = 8),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title = element_text(size = 8),
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 8),
      legend.title.position = "top",
      legend.position = "top",legend.key.size = unit(0.3,"cm"),
      plot.title = element_text(size = 8, hjust = 0.5),
      panel.grid = element_blank(),panel.border = element_blank()
    )
  return(p1)
}
G=heatplot(int1)
NF=heatplot(int2)
ILC3=heatplot(int3)
library(patchwork)
G+NF+ILC3
ggsave("./int_gene_exp_3class.pdf",width = 7.5/2.53,height = 6/2.53)