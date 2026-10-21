################################################################################
## ILC3 <-> 上皮 CellChat 结果的交叉验证 + GLCA 优势信号堆叠柱状图
##
##   1) 用 DEG_all_celltypes.csv 里 ILC3 的差异基因（GLCA vs Con）验证 CellChat
##      推断的配受体 / 效应因子方向是否一致；
##   2) 取 GLCA 中上调（delta>0）的信号通路，画堆叠柱状图（颜色=组别）。
################################################################################

setwd("E:/博士2023-2025/课题/Mettl3/7.GLCA/单细胞/rds")
suppressMessages({ library(ggplot2); library(dplyr); library(tidyr) })

group_cols <- c("Con" = "#B6CAE7", "GLCA" = "#E7CACC")

## ========================== 1. 交叉验证 ==========================
deg <- read.csv("DEG_all_celltypes.csv", stringsAsFactors = FALSE, check.names = FALSE)
ilc3 <- deg[deg$celltype == "ILC3", ]
ilc3 <- ilc3[order(-ilc3$avg_log2FC), ]
rownames(ilc3) <- ilc3$gene

cat(">>> ILC3 差异基因总数：", nrow(ilc3), "（GLCA vs Con，avg_log2FC>0 = Up）\n\n", sep = "")

## 关键基因（按功能分组）
key <- c(
  ## ILC3 效应 / 身份
  "Il22", "Il17a", "Il17f", "Csf2", "Rorc", "Il23r",
  ## CEACAM 轴
  "Ceacam1", "Ceacam2",
  ## MHC-I / 抗原递呈
  "H2-D1", "H2-K1", "H2-T23", "B2m", "Cd1d1",
  ## MHC-II（ILC3 经典抗原递呈）
  "H2-Aa", "H2-Ab1", "H2-Eb1", "Cd74",
  ## CD8
  "Cd8a", "Cd8b1",
  ## 主要差异配受体（ILC3<->Epi 中显著）
  "Cldn3", "Guca2a", "Gucy2c", "App", "Sorl1", "Nr1h4", "Cdh1",
  "F11r", "Ocln", "Nectin1", "Nectin3", "Bsg", "Tnfrsf21", "Htr4",
  "Lamb3", "Lgals9", "Hbegf", "Dsc2", "Dsg2", "P4hb", "Pparg",
  ## Up通路相关基因
  "Ephb2", "Ephb3", "Ephb4", "Efnb1", "Efnb2", "Kit", "Kitl", "Grn"
)

hit  <- key[key %in% rownames(ilc3)]
miss <- key[!key %in% rownames(ilc3)]

cv <- ilc3[hit, c("gene", "avg_log2FC", "p_val", "p_val_adj", "pct.1", "pct.2")]
cv <- cv[order(cv$avg_log2FC), ]
cv$dir <- ifelse(cv$avg_log2FC > 0, "Up",
          ifelse(cv$avg_log2FC < 0, "Down", "不变"))

# write.csv(cv, "./cellchat_ilc3_epi_deg_crossval.csv", row.names = FALSE)

cat(">>> 命中基因（", length(hit), " 个）：\n", sep = "")
cat("    gene\tavg_log2FC\tp_val_adj\tpct.1(G)\tpct.2(C)\n")
for (i in seq_len(nrow(cv))) {
  cat(sprintf("    %-10s %+8.3f  %8.2e  %5.3f  %5.3f\n",
              cv$gene[i], cv$avg_log2FC[i], cv$p_val_adj[i], cv$pct.1[i], cv$pct.2[i]))
}
cat(">>> 未在 ILC3 差异表命中（可能低表达/未检出）：",
    paste(miss, collapse = ", "), "\n\n")

## 交叉验证柱状图：log2FC（红=Up，蓝=下调）
cv$gene <- factor(cv$gene, levels = cv$gene)
p_cv <- ggplot(cv, aes(x = gene, y = avg_log2FC, fill = dir)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = c("Up" = "#E7CACC", "Down" = "#B6CAE7")) +
  coord_flip() +
  labs(x = "", y = "avg_log2FC", fill = "",
       title = "") +
  theme_bw(base_size = 6) +
  theme(legend.position = "top", legend.key.size = unit(0.3,"cm"),
        plot.title = element_text(hjust = 0.5))
p_cv
ggsave("./cellchat_ilc3_epi_deg_crossval.pdf", p_cv,
       device = cairo_pdf, width = 4, height = 8, units = "cm")
cat(">>> 交叉验证图已保存：cellchat_ilc3_epi_deg_crossval.pdf\n")

## ========================== 2. GLCA 优势信号堆叠柱状图 ==========================
pd <- read.csv("cellchat_ilc3_epi_pathway_delta.csv", stringsAsFactors = FALSE)

## 保留 GLCA 中上调（delta>0）的通路，两个方向求和得该通路总通讯概率
up <- pd %>%
  # filter(delta > 0) %>%
  group_by(pathway) %>%
  summarise(Con  = sum(prob_con),
            GLCA = sum(prob_glca),
            delta = sum(delta), .groups = "drop") %>%
  arrange(Con + GLCA)

up_long <- up %>%
  select(pathway, Con, GLCA) %>%
  pivot_longer(c(Con, GLCA), names_to = "group", values_to = "prob") %>%
  mutate(group = factor(group, levels = c("Con", "GLCA")),
         pathway = factor(pathway, levels = up$pathway))


p_stack <- ggplot(up_long, aes(x = pathway, y = prob, fill = group)) +
  geom_col(position = "fill", width = 0.6) +
  coord_flip() +                                     # 交换 x/y 轴
  scale_fill_manual(values = group_cols) +
  scale_y_continuous(labels = scales::percent,expand = c(0,0.01)) +
  labs(x = "", y = "Relative information flow", fill = "",
       title = "") +
  theme_bw(base_size = 6) +
  theme(legend.position = "top",legend.key.size = unit(0.3,"cm"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))
p_stack
ggsave("./cellchat_ilc3_epi_up_stack.pdf", p_stack,
       device = cairo_pdf, width = 12, height = 8, units = "cm")
cat(">>> 优势信号堆叠柱状图已保存：cellchat_ilc3_epi_up_stack.pdf\n")

cat("\n>>> 全部完成。\n")
