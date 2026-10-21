suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))

base_dir <- "path"
setwd(base_dir)

read_known <- function(path, grp) {
  x <- fread(path, sep = "\t", encoding = "UTF-8", fill = TRUE)
  x <- x[!is.na(`Motif Name`) & `Motif Name` != ""]
  x[, group := grp]
  x[, target_pct := as.numeric(sub("%", "", `% of Target Sequences with Motif`))]
  x[, bg_pct     := as.numeric(sub("%", "", `% of Background Sequences with Motif`))]
  x[, fold       := target_pct / bg_pct]
  x[, neglog10p  := -as.numeric(`Log P-value`) / log(10)]
  x
}

glca <- read_known("GLCA_group_encodeidr.optimal_peak.narrowPeak_homer_out/knownResults.txt", "GLCA")
pbs  <- read_known("PBS_group_encodeidr.optimal_peak.narrowPeak_homer_out/knownResults.txt", "PBS")

tf_map <- c(
  "RORgt (RORC)"      = "RORg|Rorc",
  "AHR"               = "Ahr",
  "TCF7 (TCF-1)"      = "Tcf7\\(",
  "GATA3"             = "Gata3",
  "NFIL3 (E4BP4)"     = "NFIL3",
  "RUNX1"             = "^RUNX1\\(",
  "RBPJ (Notch)"      = "Rbpj",
  "IRF4"              = "IRF4",
  "BATF"              = "BATF",
  "PU.1 (SPI1)"       = "PU\\.1",
  "AP-1 (JUN/FOS)"    = "^Jun-AP1|^AP-1",
  "NF-kB (p65/RELA)"  = "p65",
  "STAT3"             = "Stat3",
  "STAT5"             = "Stat5",
  "STAT4"             = "Stat4",
  "STAT1"             = "Stat1"
)

extract_tf <- function(x) {
  out <- lapply(names(tf_map), function(tf) {
    hit <- grep(tf_map[[tf]], x[["Motif Name"]], ignore.case = TRUE, perl = TRUE)
    if (length(hit) == 0) return(NULL)
    sel <- x[hit][which.max(neglog10p)]
    sel[, TF := tf]
    sel[, .(group, TF, `Motif Name`, `P-value`, target_pct, bg_pct, fold, neglog10p)]
  })
  rbindlist(out, use.names = TRUE, fill = TRUE)
}

res <- rbind(extract_tf(glca), extract_tf(pbs))
res[, TF := factor(TF, levels = rev(names(tf_map)))]
res[, group := factor(group, levels = c("PBS", "GLCA"))]

p <- ggplot(res, aes(x = group, y = TF)) +
  geom_point(aes(size = fold, color = neglog10p), shape = 16) +
  scale_size_continuous(
    name = "Fold enrichment\n(Target% / Background%)",
    range = c(2.5, 10),
    breaks = c(1.2, 1.5, 2.0, 2.5)
  ) +
  scale_color_gradient(
    name = expression(-log[10](italic(P))),
    low  = "#fee5d9", high = "#b30000",
    breaks = c(10, 50, 100, 200, 400)
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 8) +
  theme(
    axis.text.y = element_text(size = 8, hjust = 1),
    axis.text.x = element_text(size = 9),
    axis.ticks.length = unit(1.5, "mm"),
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 7),
    legend.key.size = unit(3, "mm"),
    legend.position = "right",
    plot.margin = margin(2, 2, 2, 2)
  )
p
ggsave("ILC3_TF_bubble.pdf", p, width = 10/2.53, height = 10/2.53)