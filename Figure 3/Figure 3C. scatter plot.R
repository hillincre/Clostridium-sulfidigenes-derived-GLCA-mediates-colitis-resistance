rm(list=ls())
setwd("path")

data=read.csv("./OPLSDA_VIP.csv")

data$sig="insig"
data$sig[data$VIP>1 & data$Corr.Coeffs.>0]="up"
data$sig[data$VIP>1 & data$Corr.Coeffs.<0]="down"

genes_to_label <- data$Metabolite[data$VIP>1]
label_data <- data[data$Metabolite %in% genes_to_label,]
library(ggplot2)
p1=ggplot(data,aes(Corr.Coeffs.,VIP))+
  geom_point(aes(fill=sig),size=2,shape=21)+
  theme_classic()+
  scale_fill_manual(values = c("#66c2a5", "grey", "#f46d43"))+
  theme(title=element_text(size=10),text=element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),legend.position = "none",
        legend.title = element_blank())+
  labs(x="Correlation Coefficients",y="VIP")
p1

library(ggrepel)
p2 = p1 +
  geom_text_repel(
    data = label_data,
    aes(label = Metabolite),
    size = 1.5,max.overlaps = 20,
    color = "black",
    show.legend = FALSE
  )
p2

pdf("./volcanol.pdf",width = 4/2.53,height = 3.5/2.53)
p2
dev.off()

library(psych)
library(tidyr)
library(dplyr)

data=read.csv("./04_AllMet.csv",row.names = 4)
data=data[,-c(1:3,grep("QC",names(data)))]
group=data.frame(group=gsub("\\_.*","",names(data)),row.names = names(data))
group$group=factor(group$group,levels = c("flox","KO"))
group=group[order(group$group),,drop=F]
data=data[,match(rownames(group),names(data))]
data=data[rownames(data)%in%label_data$Metabolite,]
data=data.frame(t(data))    
rownames(data)=gsub(".*_","",rownames(data))

ILC3_ra=read.csv("./ILC3_ratio.csv",row.names = 1)
ILC3_ra=ILC3_ra[rownames(ILC3_ra) %in% rownames(data),,drop=F]

data=data[rownames(data) %in% rownames(ILC3_ra),]
data=data[match(rownames(data),rownames(ILC3_ra)),]

res=corr.test(data,ILC3_ra,method = "pearson")
cor_df <- as.data.frame(res$r)
pvalue_df <- as.data.frame(res$p)

long_cor <- cor_df %>%
  as_tibble(rownames = "factor_2") %>%
  mutate(factor_2 = rownames(cor_df)) %>%
  pivot_longer(cols = -factor_2, names_to = "factor_1", values_to = "correlation")

long_pvalue <- pvalue_df %>%
  as_tibble(rownames = "factor_2") %>%
  mutate(factor_2 = rownames(pvalue_df)) %>%
  pivot_longer(cols = -factor_2, names_to = "factor_1", values_to = "p_value")

df=left_join(long_cor,long_pvalue[,c(1,3)],by = "factor_2")
df$factor_2=gsub("X","",df$factor_2)
df$factor_2=gsub("\\.","-",df$factor_2)
df$factor_1="ILC3 ratio"
df$stars=ifelse(abs(df$correlation)<0.4,"",
                ifelse(df$p_value < 0.01, "**",
                       ifelse(df$p_value < 0.05, "*", "")
                )
)

p1=ggplot(df, aes(factor_1,factor_2)) +
  geom_point(aes(fill=correlation,size=abs(correlation)),
             color="black",shape=21) +
  geom_text(aes(label=stars), color="black", size=2) + 
  scale_fill_gradient2(low='#4671e0', high='#e04146',mid = 'white', limit=c(-1,1)) +
  labs(x=NULL,y=NULL) + theme_bw(base_size = 8)+
  theme(panel.border=element_rect(fill = NULL,linewidth = 0.5),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.key.size = unit(0.3,"cm"),
        legend.position = "right" )
p1

ggsave(plot = p1,"./metabo_ILC3_cor.pdf",width = 5/2.53,height = 8/2.53)