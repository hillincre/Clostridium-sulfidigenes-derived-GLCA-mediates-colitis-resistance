library(ggplot2)
library(dplyr)
library(patchwork)

setwd("path")
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

summ_wide <- summ %>%
  pivot_wider(
    names_from  = group,
    values_from = c(mean, sem),
    names_glue  = "{group}_{.value}"
  ) %>%
  select(Bile_Acid, PBS_mean, PBS_sem, CS_mean, CS_sem)

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

all=left_join(or_df,pval_df,by = "Bile_Acid")
all=left_join(all,summ_wide,by = "Bile_Acid")
write.csv(all,"./BA_screen_data.csv")

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
  scale_fill_manual(values = c("PBS" = col_pbs, "CS" = col_cs)) +
  labs(x = NULL, y = "Abundance (mean ± SEM)", title = "") +
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
  plot_layout(ncol = 3, widths = c(1.4, 1, 1))&theme(plot.margin = margin(0,0,0,0))
p_comb
ggsave("bile_acid_3panel.pdf", p_comb, width = 12/2.53, height = 5/2.53)