setwd("path")

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)


obj=readRDS("./anchored_data_processed.rds")
group_cols <- c("Con" = "#7b95c6", "GLCA" = "#f0486f")
color=c("#7b95c6","#a1d8e8","#67a583","#d0e2c0",
        "#ffcc1a","#f47254","#f0486f","#c960f0")
n=length(unique(Idents(seurat)))
color=colorRampPalette(color)(n)
names(color)=levels(Idents(seurat))
p_umap <- DimPlot(obj, reduction = "umap", group.by = "celltype",cols=color,
                  split.by = "group", 
                  label = TRUE, repel = TRUE,
                  label.size = 3, raster = FALSE) +
  NoLegend() + labs(title="")+
  theme(text = element_text(size=10),
        legend.key.size = unit(0.35,"cm"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        aspect.ratio = 1
  )
p_umap

markers <- FindAllMarkers(obj, assay = "RNA", slot = "data",
                          only.pos = TRUE
)

features <- c(
  "Cd19","Ms4a1","Cd79a","Cd79b",     
  "Cd3e","Cd3d","Cd4",                
  "Cd3e","Cd8a","Cd8b","Gzmb","Nkg7", 
  "Flt3","H2-Aa",      
  "Cdh5","Pecam1","Kdr","Vwf",        
  "Epcam","Krt8","Krt18","Cdh1",      
  "Col1a1","Col1a2","Dcn","Pdgfra",   
  "Ncr1","Klrb1c","Ifng","Tbx21",     
  "Gata3","Il1rl1","Il13","Rora",     
  "Rorc","Il23r","Il17a","Kit",       
  "Cd68","Adgre1","Lyz2","Csf1r",     
  "Ncr1","Nkg7","Klrd1","Prf1"        
)

cells <- c("B cells","Cd4+ T","Cd8+ T","Dendritic cells","Endothelial cells",
           "Epithelial cells","Fibroblasts","ILC1","ILC2","ILC3","Macrophage","NK cells")

df <- markers[markers$gene %in% features & markers$avg_log2FC > 0, ]
df$cluster <- factor(as.character(df$cluster), levels = cells)
df$gene    <- factor(df$gene, levels = unique(features))

p_marker <- ggplot(df, aes(x = gene, y = cluster)) +
  geom_point(aes(size = pct.1, fill = avg_log2FC), shape = 21) +
  scale_fill_gradientn(colors = c("#31479F","#1CB4C5","#85C047","#F7BA07","#FF564C"),
                       name = "avg_log2FC") +
  scale_size_continuous(name = "Percent", range = c(1, 5)) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5),
        legend.key.size = unit(0.3, "cm")) +
  labs(x = "", y = "")
p_marker
ggsave("./marker_pointplot.pdf", p_marker, device = cairo_pdf, width = 17, height = 6, units = "cm")

immu=subset(obj,!celltype%in%c("Endothelial cells","Epithelial cells","Fibroblasts","Dendritic cells"))
Ratio <- immu@meta.data %>%
  group_by(group, celltype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(group) %>%
  mutate(count = round(n / sum(n) * 100,2))

library(ggpie)
col=color[!names(color)%in%c("Endothelial cells","Epithelial cells","Fibroblasts","Dendritic cells")]
df_con  <- Ratio[Ratio$group == "Con", ]
df_glca <- Ratio[Ratio$group == "GLCA", ]

p1 <- ggdonut(df_con,
              group_key = "celltype",
              donut.label = F,
              label_size = 2, 
              label_pos = "out",
              label_gap = 0.1,
              count_type = "count",
              border_size = 0,
              fill_color = col,
              r0 = 1.5,
              r1 = 3) +
  theme(legend.position = "none") +
  geom_text(aes(x=0,y=0,label="Con"),size=8,size.unit = "pt")
p1
p2 <- ggdonut(df_glca,
              group_key = "celltype",
              donut.label = F,
              label_size = 2, 
              label_pos = "out",
              label_gap = 0.1,
              count_type = "count",
              border_size = 0,
              fill_color = col,
              r0 = 1.5,
              r1 = 3) +
  theme(legend.position = "none") +
  geom_text(aes(x=0,y=0,label="GLCA"),size=8,size.unit = "pt")

library(patchwork)
p_donut=p1 / p2 + plot_layout(guides = "collect") &
  theme(legend.position='right',legend.key.size = unit(0.3,"cm"),
        text = element_text(size = 8),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.3))
ggsave("./celltype_ratio.pdf", p_donut, device = cairo_pdf, width = 8, height = 6, units = "cm")




plot_volcano <- function(de) {
  de$regulation <- "no"
  de$regulation[de$avg_log2FC >  0 & de$p_val < 0.05] <- "up"
  de$regulation[de$avg_log2FC < 0 & de$p_val < 0.05] <- "down"
  
  up <- sum(de$regulation == "up"); down <- sum(de$regulation == "down")
  ggplot(de, aes(avg_log2FC, -log10(p_val))) +
    geom_point(aes(color = regulation), alpha = 0.5, size = 0.5) +
    scale_color_manual(values = c("down" = "#66c2a5", "no" = "grey", "up" = "#f46d43")) +
    labs(x = "log2(FC)", y = "-log10(p)",
         title = paste0("Con vs GLCA\nDown: ", down, "  Up: ", up)) +
    theme_bw(base_size = 8) +
    theme(legend.position = "none", panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5))
}

library(clusterProfiler)
library(org.Mm.eg.db)
do_go <- function(de) {
  genes <- de$gene[de$avg_log2FC > 0 & de$p_val < 0.05]
  if (length(genes) < 5) return(NULL)
  df_id <- bitr(genes, "SYMBOL", "ENTREZID", OrgDb = org.Mm.eg.db)
  enrichKEGG(gene = df_id$ENTREZID, "mmu",
             pvalueCutoff = 0.05, keyType = "kegg")
  
}

for (ct in names(de_list)) {
  de <- de_list[[ct]]
  tryCatch(
    ggsave(paste0("./DEG_volcano_", gsub("[ +]", "_", ct), ".pdf"), plot_volcano(de),
           width = 6, height = 6, units = "cm"),
    error = function(e) cat("skip volcano", ct, ":", conditionMessage(e), "\n")
  )
  go <- tryCatch(do_go(de), error = function(e) NULL)
  if (!is.null(go) && nrow(go@result) > 0)
    write.csv(go@result, paste0("./DEG_KEGG2_", gsub("[ +]", "_", ct), ".csv"))
}


p1=plot_volcano(de)
kegg=read.csv("./kegg_int.csv",row.names = 1)
kegg=kegg[order(kegg$RichFactor,decreasing = F),]
kegg$Description=factor(kegg$Description,levels = kegg$Description)
p2=ggplot(kegg,aes(x=RichFactor,y=Description))+
  geom_col(fill="grey",color="black")+theme_bw(base_size = 8)+
  geom_text(aes(x=0,label=Description),hjust=0,size=8,size.unit = "pt")+
  labs(title = "KEGG")+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),plot.title = element_text(hjust = 0.5))
p2
library(patchwork)
p=p1+p2+
  plot_layout(ncol=2,widths = c(1.5,2))
p
ggsave(plot=p,"./ILC3_DEGs.pdf",units = "cm",width = 10,height = 6)
