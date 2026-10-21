setwd("path")
library(ggplot2)
library(ggridges)
library(dplyr)
library(tidyr)
library(readxl)
library(patchwork)

df <- read_excel("target_BA_abun.csv")

df_long <- df %>%
  pivot_longer(cols = starts_with("sample"),
               names_to = "sample",
               values_to = "value") %>%
  mutate(value = as.numeric(value))

group_colors <- c("Control"        = "#B6CAE7",
                  "C.sulfidigenes" = "#E7CACC",
                  "Vector"         = "#7ABBCF",
                  "BSH"            = "#D7B3D5")

make_panel <- function(ba, g1, g2) {
  sub <- df_long %>%
    filter(BA == ba, group %in% c(g1, g2)) %>%
    mutate(group = factor(group, levels = c(g1, g2)))
  
  wt <- t.test(value ~ group,alternative ="less", data = sub)
  p  <- wt$p.value
  p_lab <- ifelse(p < 0.0001, "p < 0.0001", paste0("p = ", signif(p, 3)))
  
  ggplot(sub, aes(x = value, y = group, fill = group)) +
    geom_density_ridges(alpha = 0.7, scale = 2.5,
                        quantile_lines = TRUE,
                        quantiles = 0.5, ) +
    scale_fill_manual(values = group_colors) +
    scale_y_discrete(expand = expansion(mult = c(0.05, 0.35))) +
    labs(x = paste0(ba, " abundance"), y = NULL,
         title = "") +
    
    coord_cartesian(clip = "off") +
    annotate("text", x = max(sub$value), y = 1.5,
             label = p_lab, hjust = 1, vjust = 0, size = 2) +
    theme_minimal(base_size = 6) +
    theme(
      legend.position = "right",legend.title = element_blank(),
      legend.key.size = unit(0.3,"cm"),
      axis.text.y = element_blank(),
      panel.border = element_rect(fill = NULL),
      plot.title = element_text(hjust = 0.5, size = 6))
}

p1 <- make_panel("LCA",  "Control", "C.sulfidigenes")
p2 <- make_panel("LCA",  "Vector",  "BSH")
p3 <- make_panel("GLCA", "Control", "C.sulfidigenes")
p4 <- make_panel("GLCA", "Vector",  "BSH")

combined <- (p1 | p3 | p2 | p4) &theme(legend.position = "none")+plot_layout(guides = "collect")
combined

ggsave("ridge_plot.pdf", combined, width = 10/2.53, height = 2.9/2.53, dpi = 300)