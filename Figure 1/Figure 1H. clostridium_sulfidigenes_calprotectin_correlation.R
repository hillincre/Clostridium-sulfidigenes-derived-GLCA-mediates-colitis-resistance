library(ggplot2)
library(ggpubr)
library(dplyr)
library(data.table)

setwd("path")

# ---------- 1. ----------
# metadata: Sample + Fecal.Calprotectin
meta <- read.delim("metadata.tsv", stringsAsFactors = FALSE, check.names = FALSE)

target <- "d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium_J;s__Clostridium_J sulfidigenes"

counts <- fread("species.counts.tsv",
                select = c("Sample", target),
                header = TRUE, check.names = FALSE)
counts <- as.data.frame(counts)

# ---------- 2. ----------
dat <- merge(meta[, c("Sample", "Fecal.Calprotectin")], counts, by = "Sample")

names(dat)[names(dat) == target] <- "Abundance"
dat$Fecal.Calprotectin <- suppressWarnings(as.numeric(dat$Fecal.Calprotectin))
dat$Abundance          <- suppressWarnings(as.numeric(dat$Abundance))

dat <- dat[!is.na(dat$Fecal.Calprotectin), ]

# ---------- 3. ----------
p <- ggscatter(dat,
               x = "Abundance", y = "Fecal.Calprotectin",
               color = "black",
               add = "reg.line",
               conf.int = TRUE,
               add.params = list(color = "red", fill = "grey"),
               cor.coef = TRUE,
               cor.method = "spearman",
               size = 1.5,
               cor.coef.size = 3) +
  labs(y = "Fecal Calprotectin",
       x = "Clostridium_J sulfidigenes abundance") +
  theme_bw(base_size = 9) +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", size = 0.3),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5))

p

ggsave("clostridium_sulfidigenes_calprotectin_correlation.pdf",
       p, width = 4, height = 3.5)