setwd("path")

library(CellChat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ComplexHeatmap)
library(Seurat)
obj=readRDS("./anchored_data_processed.rds")

cellchat_con  <- readRDS("./cellchat_con.rds")
cellchat_glca <- readRDS("./cellchat_glca.rds")

immune_cts <- c("B cells","Cd4+ T","Cd8+ T","Dendritic cells",
                "ILC1","ILC2","ILC3","Macrophage")
epi <- "Epithelial cells"
sub_cells <- c(immune_cts, epi)

group_cols <- c("Con" = "#B6CAE7", "GLCA" = "#E7CACC")
col=c( "#93BFDB","#6CA98C","#D8DDA1","#F6D038",  "#FAA334",  "#F36E56",  "#F15367",  "#C960F0","grey" )
names(col)=c("Macrophage","B cells","ILC1","ILC2","Cd8+ T","ILC3","Cd4+ T","NK cells","Epithelial cells")

plot_immune_epi_net <- function(cc, measure = c("count", "weight"), title = "") {
  measure <- match.arg(measure)
  mat <- cc@net[[measure]][sub_cells, sub_cells, drop = FALSE]
  mat[immune_cts, immune_cts] <- 0
  size <- as.numeric(table(cc@idents)[sub_cells])
  netVisual_circle(mat, vertex.weight = size,edge.label.cex = 0.5,vertex.label.cex = 0.5,color.use = col, weight.scale = TRUE,
                   label.edge = FALSE, title.name = title)
}

pdf("./cellchat_immune_epi_network.pdf", width = 3, height = 3)
plot_immune_epi_net(cellchat_con,  "count",  "Con  - Number of interactions")
plot_immune_epi_net(cellchat_glca, "count",  "GLCA - Number of interactions")
plot_immune_epi_net(cellchat_con,  "weight", "Con  - Interaction strength")
plot_immune_epi_net(cellchat_glca, "weight", "GLCA - Interaction strength")
dev.off()

get_epi <- function(cc) {
  cm <- cc@net$count; wm <- cc@net$weight
  do.call(rbind, lapply(immune_cts, function(ct) {
    data.frame(celltype = ct,
               count  = cm[ct, epi] + cm[epi, ct],
               weight = wm[ct, epi] + wm[epi, ct])
  }))
}

df <- rbind(transform(get_epi(cellchat_con),  group = "Con"),
            transform(get_epi(cellchat_glca), group = "GLCA"))
df$group    <- factor(df$group, levels = c("Con", "GLCA"))
df$celltype <- factor(df$celltype, levels = immune_cts)

bar_theme <- theme_bw(base_size = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top", legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))

p1 <- ggplot(df, aes(x = celltype, y = count, fill = group)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  scale_fill_manual(values = group_cols) +
  labs(x = "", y = "Number of interactions",
       title = "") + bar_theme+theme(axis.text.x = element_blank())

p2 <- ggplot(df, aes(x = celltype, y = weight, fill = group)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  scale_fill_manual(values = group_cols) +
  labs(x = "", y = "Interaction strength",
       title = "") + bar_theme

p_bar <- p1 / p2 + plot_layout(guides = "collect")&theme(legend.key.size = unit(0.3,"cm"),legend.position = "right")
p_bar
ggsave("./cellchat_immune_epi_bar.pdf", p_bar, device = cairo_pdf,
       width = 6, height = 6, units = "cm")

ht_con  <- netVisual_heatmap(cellchat_con,  measure = "count", title.name = "")
ht_glca <- netVisual_heatmap(cellchat_glca, measure = "count", title.name = "")
pdf("./cellchat_heatmap_compare.pdf", width = 10/2.53, height = 15/2.53)
ComplexHeatmap::draw(ht_con %v% ht_glca)
dev.off()

