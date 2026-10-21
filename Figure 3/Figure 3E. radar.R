setwd("path")

library(ggplot2)
library(dplyr)
library(patchwork)

data_path <- "bile_acid_ase_abun.csv"
df <- read.csv(data_path, header = TRUE, stringsAsFactors = FALSE)

df <- df[, c("Bile_Acid", "Sample", "Abundance", "group")]
df$group <- factor(df$group, levels = c("CON", "CS"))

ba_order <- df %>%
  group_by(Bile_Acid) %>%
  summarise(mean_abun = mean(Abundance, na.rm = TRUE)) %>%
  arrange(desc(mean_abun)) %>%
  pull(Bile_Acid)

ba_levels <- rev(ba_order)

summ <- df %>%
  group_by(Bile_Acid, group) %>%
  summarise(
    mean = mean(Abundance, na.rm = TRUE),
    sem  = sd(Abundance, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(Bile_Acid = factor(Bile_Acid, levels = ba_levels))

or_df <- do.call(rbind, lapply(split(df, df$Bile_Acid), function(d) {
  d$z  <- as.numeric(scale(d$Abundance))
  d$cs <- as.integer(d$group == "CS")
  fit <- glm(cs ~ z, data = d, family = binomial())
  s   <- summary(fit)$coefficients
  beta <- s["z", "Estimate"]
  se   <- s["z", "Std. Error"]
  data.frame(
    Bile_Acid = d$Bile_Acid[1],
    OR    = exp(beta),
    lower = exp(beta - 1.96 * se),
    upper = exp(beta + 1.96 * se)
  )
}))
or_df$Bile_Acid <- factor(or_df$Bile_Acid, levels = ba_levels)

pval_df <- df %>%
  group_by(Bile_Acid) %>%
  summarise(p = wilcox.test(Abundance ~ group, exact = FALSE)$p.value,
            .groups = "drop") %>%
  mutate(
    Bile_Acid = factor(Bile_Acid, levels = ba_levels),
    neglog10p = -log10(p)
  )

col_pbs <- "#B6CAE7"
col_cs  <- "#E7CACC"

use_log_scale <- TRUE

p1 <- ggplot(summ, aes(x = Bile_Acid, y = mean, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72) +
  geom_errorbar(
    aes(ymin = mean - sem, ymax = mean + sem),
    position = position_dodge(width = 0.8),
    width = 0.25
  ) +
  coord_flip() +
  scale_fill_manual(values = c("CON" = col_pbs, "CS" = col_cs)) +
  labs(x = NULL, y = "Log10(Abundance)", title = "") +
  theme_minimal(base_size = 6) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "left",legend.key.size = unit(0.3,"cm"),
    legend.title       = element_blank(),
    panel.border = element_rect(fill = NULL)
  )
if (use_log_scale) p1 <- p1 + scale_y_log10()
p1

p2 <- ggplot(or_df, aes(x = OR, y = Bile_Acid)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.3, colour = "grey30") +
  geom_point(size = 2.2, colour = col_pbs) +
  scale_x_log10() +
  labs(x = "OR (95% CI)", y = NULL, title = "") +
  theme_minimal(base_size = 6) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        panel.border = element_rect(fill = NULL))
p2

p3 <- ggplot(pval_df, aes(x = neglog10p, y = Bile_Acid)) +
  geom_col(fill = col_pbs) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  labs(x = "-log10(pvalue)", y = NULL, title = "") +
  theme_minimal(base_size = 6) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        panel.border = element_rect(fill = NULL))

p2 <- p2 + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
p3 <- p3 + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

p_comb <- p1 + p2 + p3 +
  plot_layout(ncol = 3, widths = c(1.4, 1, 1))
p_comb
ggsave("bile_acid_3panel.pdf", p_comb, width = 7/2.53, height = 6/2.53)

fc_data <- long_data %>%
  group_by(Bile_Acid, group) %>%
  summarise(MeanAbun = median(Abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = MeanAbun) %>%
  mutate(FoldChange = CS / CON) %>%
  select(Bile_Acid, FoldChange)

fc_data=fc_data[fc_data$Bile_Acid %in% c("3-DHCA","GDCA","GLCA"),]
fc_data=fc_data[order(fc_data$FoldChange,decreasing = F),]
fc_data$Bile_Acid=factor(fc_data$Bile_Acid,levels = fc_data$Bile_Acid)
p_fc=ggplot(fc_data,aes(y=Bile_Acid,x=log2(FoldChange)))+
  geom_col(fill="#5fe87f")+
  labs(x="log2(FC)",y="")+
  scale_x_continuous(expand = expansion(mult = c(0,0.1)),limits = c(0,1),
                     breaks = c(0,0.5,1))+
  theme_bw(base_size = 8)+
  theme(panel.grid = element_blank(),axis.line = element_line(),
        panel.border = element_rect(fill=NULL),
        legend.position = "none")

library(tidyverse)
library(randomForest)
library(pROC)
library(caret)
library(ggpubr)

wide_data <- long_data %>%
  filter(Bile_Acid %in% c("3-DHCA","GDCA","GLCA")) %>%
  select(Bile_Acid, Sample, Abundance, group) %>%
  pivot_wider(names_from = Bile_Acid, values_from = Abundance, values_fill = 0) %>%
  column_to_rownames("Sample")

X <- wide_data %>% select(-group)
y <- factor(wide_data$group)

keep <- apply(X, 2, var) > 0
X <- X[, keep, drop = FALSE]

cat("用于建模的胆汁酸数量：", ncol(X), "\n")

set.seed(123)

rf_model <- randomForest(
  x = X,
  y = y,
  ntree      = 1000,
  mtry       = floor(sqrt(ncol(X))),
  importance = TRUE
)

print(rf_model)

imp <- importance(rf_model, scale = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column("Bile_Acid")

imp$Bile_Acid=factor(imp$Bile_Acid,levels = levels(fc_data$Bile_Acid))
p_rf <- imp %>%
  slice_max(MeanDecreaseAccuracy , n = 20) %>%
  ggplot(aes(MeanDecreaseAccuracy , Bile_Acid )) +
  geom_col( fill = "#9ECAE1") +
  labs(x = "MeanDecreaseAccuracy", y = NULL,
       title = "") +
  scale_x_continuous(expand = expansion(mult = c(0,0.1)))+
  theme_bw(base_size = 8)+
  theme(panel.grid = element_blank(),
        axis.line = element_line(),        
        axis.text.y = element_blank(),
        panel.border = element_rect(fill=NULL),
        legend.position = "none")

print(p_rf)

library(psych)
library(tidyr)
library(dplyr)
cal=read.csv("./Fecal Calprotectin.csv",row.names = 1)
ba=wide_data[rownames(wide_data) %in% rownames(cal),-1]

res=corr.test(ba,cal,method = "spearman")
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
df$factor_1="Fecal Calprotectin"
df$stars=ifelse(abs(df$correlation)<0.4,"",
                ifelse(df$p_value < 0.01, "**",
                       ifelse(df$p_value < 0.05, "*", "")
                )
)

df$factor_2=factor(df$factor_2,levels = levels(fc_data$Bile_Acid))

p1=ggplot(df, aes(factor_1,factor_2)) +
  geom_point(aes(fill=correlation,size=abs(correlation)),
             color="black",shape=21) +
  geom_text(aes(label=stars), color="black", size=2) + 
  scale_fill_gradient2(low='#4671e0', high='#e04146',mid = 'white', limit=c(-1,1))+
  labs(x=NULL,y=NULL) + theme_bw(base_size = 8)+
  theme(panel.border=element_rect(fill = NULL,linewidth = 0.5),
        panel.background = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        legend.key.size = unit(0.3,"cm"),
        legend.position = "right" )
p1

library(dplyr)

pvalue_df2 <- pvalue_df %>%
  mutate(Bile_Acid = as.character(imp$Bile_Acid), .before = 1) %>%
  rename(Fecal_Calprotectin_pvalue = Fecal_Calprotectin)

cor_df2 <- cor_df %>%
  mutate(Bile_Acid = as.character(imp$Bile_Acid), .before = 1) %>%
  rename(Fecal_Calprotectin_cor = Fecal_Calprotectin)

fc_data2 <- fc_data %>%
  mutate(Bile_Acid = as.character(Bile_Acid))

imp2 <- imp %>%
  mutate(Bile_Acid = as.character(Bile_Acid))

out <- imp2 %>%
  left_join(fc_data2, by = "Bile_Acid") %>%
  left_join(pvalue_df2, by = "Bile_Acid") %>%
  left_join(cor_df2, by = "Bile_Acid")
out
write.csv(out,"./radar_data.csv")
library(dplyr)
library(fmsb)

radar_vars <- c(
  "MeanDecreaseAccuracy",
  "MeanDecreaseGini",
  "FoldChange",
  "Fecal_Calprotectin_pvalue"
)

radar_df <- out %>%
  select(Bile_Acid, all_of(radar_vars)) %>%
  mutate(
    Fecal_Calprotectin_pvalue = -log10(Fecal_Calprotectin_pvalue)
  ) %>%
  mutate(
    across(
      all_of(radar_vars),
      ~ (. - min(.)) / (max(.) - min(.))
    )
  )

radar_df=radar_df %>%
  left_join(cor_df2, by = "Bile_Acid")
radar_df$Fecal_Calprotectin_cor=-radar_df$Fecal_Calprotectin_cor
plot_data <- rbind(
  rep(1, 5),
  rep(0, 5),
  radar_df %>% select(-Bile_Acid) %>% as.data.frame()
)

colnames(plot_data) <- c(
  "MeanDecreaseAccuracy",
  "MeanDecreaseGini",
  "FoldChange",
  "-log10(pvalue)",
  "cor(-r)"
)

pdf("./radar_screening.pdf",width=8/1.2,height = 5/1.2)
radarchart(
  plot_data,
  axistype = 4,
  pcol = c("steelblue", "darkorange", "darkgreen"),
  pfcol = c(
    rgb(0.27, 0.51, 0.71, 0.3),
    rgb(1.00, 0.55, 0.00, 0.3),
    rgb(0.00, 0.39, 0.00, 0.3)
  ),
  plwd = 2,
  cglcol = "grey70",
  cglty = 1,
  axislabcol = "grey30",
  vlcex = 0.8,
  title = ""
)

legend(
  "topright",
  legend = radar_df$Bile_Acid,
  col = c("steelblue", "darkorange", "darkgreen"),
  lwd = 3,
  bty = "n"
)
dev.off()