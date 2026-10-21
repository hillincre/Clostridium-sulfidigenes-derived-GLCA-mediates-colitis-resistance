
suppressMessages({
  library(data.table)
  library(ggplot2)
  library(gridExtra)
  library(grid)
})

data_dir <- "TSS_enrichment"

prof     <- fread(file.path(data_dir, "profiles.csv"))
mat_pbs  <- fread(file.path(data_dir, "matrix_PBS.csv"))
mat_glca <- fread(file.path(data_dir, "matrix_GLCA.csv"))

m_pbs  <- as.matrix(mat_pbs[, -1, with = FALSE])
m_glca <- as.matrix(mat_glca[, -1, with = FALSE])
pos    <- prof$pos

core <- which(pos >= -100 & pos <= 100)
m_pbs  <- m_pbs[order(rowMeans(m_pbs[, core]), decreasing = FALSE), , drop = FALSE]
m_glca <- m_glca[order(rowMeans(m_glca[, core]), decreasing = FALSE), , drop = FALSE]

C_PBS  <- "#b4c6e2"
C_GLCA <- "#e4c8ca"
pal    <- c("#ffffff", "#fcbba1", "#fb6a4a", "#cb181d")
zmax   <- quantile(c(as.vector(m_pbs), as.vector(m_glca)), 0.99, na.rm = TRUE)
heat_rows <- 2000

df_pbs  <- data.frame(pos = pos, signal = prof$PBS_mean,  se = prof$PBS_se)
df_glca <- data.frame(pos = pos, signal = prof$GLCA_mean, se = prof$GLCA_se)
ymax <- max(max(df_pbs$signal + df_pbs$se),
            max(df_glca$signal + df_glca$se)) * 1.05

mk_hist <- function(df, col, title) {
  ggplot(df, aes(x = pos, y = signal)) +
    geom_line(colour = col, linewidth = 1)+
    geom_errorbar(aes(ymin = pmax(0, signal - se), ymax = signal + se),
                  width = 0, colour = "grey30", linewidth = 0.25) +
    geom_vline(xintercept = 0, colour = "grey40",
               linetype = "dashed", linewidth = 0.4) +
    scale_x_continuous(breaks = c(-2000, -1000, 0, 1000, 2000),
                       limits = c(-2000, 2000), expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    coord_cartesian(ylim = c(0, ymax)) +
    labs(title = title, x = "Distance from TSS (bp)",
         y = "Fold-enrichment signal") +
    theme_classic(base_size = 6) +
    theme(plot.title = element_text(hjust = 0.5,  size = 6),
          axis.text  = element_text(colour = "grey20"))
}

aggregate_rows <- function(m, target = 2000) {
  nr <- nrow(m)
  if (nr <= target) return(m)
  g   <- ceiling(seq_len(nr) * target / nr)
  s   <- rowsum(m, g)
  cnt <- tabulate(g, nbins = target)
  s / cnt
}

make_heat_df <- function(m) {
  nr <- nrow(m); nb <- ncol(m)
  data.frame(
    y = rep(seq_len(nr), times = nb),
    x = rep(pos, each = nr),
    z = as.vector(pmin(m, zmax))
  )
}

mk_heat <- function(m, title) {
  m  <- aggregate_rows(m, heat_rows)
  df <- make_heat_df(m)
  ggplot(df, aes(x = x, y = y, fill = z)) +
    geom_raster() +
    scale_fill_gradientn(
      colours = pal, limits = c(0, zmax),
      breaks = c(0, zmax / 2, zmax),
      labels = c("0", sprintf("%.1f", zmax / 2), sprintf("%.1f", zmax)),
      guide = guide_colorbar(barwidth = 0.6, barheight = 5,
                             title = "fc signal")) +
    scale_x_continuous(breaks = c(-2000, -1000, 0, 1000, 2000),
                       expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(title = title, x = "", y = "Genes (ranked)") +
    theme_classic(base_size = 6) +
    theme(plot.title   = element_text(hjust = 0.5,  size = 6),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          legend.key.size = unit(0.3,"cm"),
          axis.line = element_blank(),
          axis.text.x  = element_text(colour = "grey20",size = 6))
}

p1 <- mk_hist(df_pbs,  C_PBS,  "PBS")
p2 <- mk_hist(df_glca, C_GLCA, "GLCA")
h1 <- mk_heat(m_pbs,  "")
h2 <- mk_heat(m_glca, "")

out_pdf <- file.path(data_dir, "TSS_enrichment_2x2_2.pdf")

draw_all <- function() {
  grid.arrange(
    ggplotGrob(p1), ggplotGrob(p2), ggplotGrob(h1), ggplotGrob(h2),
    layout_matrix = matrix(1:4, nrow = 2, byrow = TRUE),
    heights = c(1, 3)
  )
}

pdf(out_pdf, width = 10/2.53, height = 10/2.53)
draw_all()
dev.off()

