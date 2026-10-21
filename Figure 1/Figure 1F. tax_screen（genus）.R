library(ggplot2)
library(dplyr)
library(tidyr)
setwd("path")

##contribution
con=read.csv("./Contribute_Genus_KEGG_Enzyme.csv",row.names = 2)
df=con[-1,-1]
df=as.data.frame(t(df))
df=df[order(rownames(df),decreasing = F),]

df$group=c(rep("DSS",7),rep("Control",8))
df$group=factor(df$group,levels = c("DSS","Control"))

df2=pivot_longer(df,cols = c(1:50),names_to="Tax",values_to = "Contribution")

df_summ=df2 %>%
  group_by(group,Tax) %>%
  summarise(mean=mean(Contribution),
            SD=sd(Contribution),
            .groups = "drop"
  )


abun=read.csv("./tax_g_rpkm.csv",row.names = 1)
abun=abun[,colnames(abun)!="Total"]
rownames(abun) <- sub("^.*(g__)", "\\1", rownames(abun))
abun=abun[,order(colnames(abun),decreasing = F)]
abun=abun[,!colnames(abun) %in% c("A216","A218","A231","A233","A234")]
group=df[,51,drop=F]
group=group[rownames(group) %in% colnames(abun),,drop=F]

library(dplyr)

group1 <- "Control"
group2 <- "DSS"

safe_wilcox <- function(g1, g2) {
  tryCatch({
    wt <- wilcox.test(g1, g2, conf.int = TRUE)
    list(p = wt$p.value, lower = wt$conf.int[1], upper = wt$conf.int[2])
  }, error = function(e) {
    list(p = NA, lower = NA, upper = NA)
  })
}

result <- lapply(rownames(abun), function(tax){
  val <- as.numeric(abun[tax, ])
  names(val) <- colnames(abun)
  
  g1 <- val[group$group == group1]
  g2 <- val[group$group == group2]
  
  if (length(g1) == 0 | length(g2) == 0) return(NULL)
  
  mean1 <- mean(g1, na.rm = TRUE)
  sd1 <- sd(g1, na.rm = TRUE)
  mean2 <- mean(g2, na.rm = TRUE)
  sd2 <- sd(g2, na.rm = TRUE)
  fold <-  (mean1+0.0001) / (mean2+0.0001)
  
  wt <- safe_wilcox(g1, g2)
  pval <- wt$p
  lower.ci <- wt$lower
  upper.ci <- wt$upper
  
  n1 <- length(g1)
  n2 <- length(g2)
  z <- qnorm(pval / 2)
  r <- ifelse(is.na(pval), NA, abs(z) / sqrt(n1 + n2))
  
  data.frame(
    Name = tax,
    Control.Mean = mean1,
    Control.Sd = sd1,
    DSS.Mean = mean2,
    DSS.Sd = sd2,
    Pvalue = pval,
    Fold.Change = fold,
    Lower.ci = lower.ci,
    Upper.ci = upper.ci,
    Effectsize = r
  )
})

res_df <- do.call(rbind, result[!sapply(result, is.null)])
res_df$Corrected.pvalue <- p.adjust(res_df$Pvalue, method = "fdr")

res_df <- res_df %>%
  select(Name, Control.Mean, Control.Sd, DSS.Mean, DSS.Sd, Pvalue, Fold.Change, Corrected.pvalue, Lower.ci, Upper.ci, Effectsize) %>%
  mutate(across(where(is.numeric), round, 5))



########plot########
data=left_join(df_summ,res_df,by = join_by("Tax"=="Name"))
data=data[order(data$mean,decreasing = F),]
data$label=ifelse(data$Pvalue>0.05,"",
                  ifelse(data$Pvalue>0.01,"*",
                         ifelse(data$Pvalue>0.001,"**","***")))
data=data[!is.na(data$Pvalue),]

library(ggplot2)
library(dplyr)
library(ggtext)
library(patchwork)

data_clean <- data %>%
  filter(!is.na(Tax), !is.na(mean), !is.na(Fold.Change),group=="DSS") %>%
  arrange(mean) %>%
  arrange(Fold.Change) %>%
  mutate(Tax = factor(Tax, levels = unique(Tax)))
data_clean=data_clean[-grep("sp.",data_clean$Tax),]
write.csv(data_clean,"tax_screnn_data.csv")

data_clean=tail(data_clean,10)
p1 <- ggplot(data_clean, aes(x = mean, y = Tax)) +
  geom_col(fill = "#3B82F6", width = 0.7) +
  labs(x = "Contribution to BSH", y = "") +
  theme_bw(base_size = 6) +
  theme(
    axis.text.y = element_blank(),
    panel.grid = element_blank()
  )
p1

p2 <- ggplot(data_clean, aes(x = Fold.Change, y = Tax)) +
  geom_col(fill = "#10B981", width = 0.7) +
  labs(x = "log2(FC)", y = "") +
  theme_bw(base_size = 6) +
  theme(
    axis.text.y = element_text(size = 6),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank()
  )


p3 <- ggplot(data_clean, aes(x = 0, y = Tax)) +
  geom_text(aes(label = label), size = 3, hjust = 0,vjust=0.75, color = "#1E293B") +
  labs(x = "Significance", y = "") +
  xlim(-0.1, 1) +
  theme_bw(base_size = 6) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank()
  )

combined_plot <-  p2 + p3 + p1 +
  plot_layout(widths = c(4, 1.5, 5)) 

print(combined_plot)

# ggsave("./genus_diff&con.pdf",plot = combined_plot,width = 12/2.53,height = 15/2.53)






#########clostridium abundance##########
abun_clos=read.csv("./tax_s_rpkm.csv")
abun_clos=abun_clos[grep("g__Clostridium",abun_clos$X.Taxonomy),]
abun_clos=abun_clos[-grep("sp.",abun_clos$X.Taxonomy),]
abun_clos=abun_clos[,colnames(abun_clos)!="Total"]
abun_clos=pivot_longer(abun_clos,cols = c(2:16),names_to = "sample",values_to = "abundance")
group=data.frame(sample=rownames(df),group=df$group)
abun_clos=left_join(abun_clos,group,by = "sample")

abun_clos$X.Taxonomy=gsub("d__Bacteria;k__Bacillati;p__Bacillota;c__Clostridia;o__Eubacteriales;f__Clostridiaceae;g__Clostridium;",
                          "",abun_clos$X.Taxonomy)
head(abun_clos)



library(tidyverse)

# ======================

diff_result <- abun_clos %>%
  
  group_by(X.Taxonomy) %>%
  
  filter(length(unique(group)) == 2) %>%
  
  summarise(
    
    mean_Control = mean(abundance[group == "Control"], na.rm = T),
    mean_mettl3ff  = mean(abundance[group == "DSS"], na.rm = T),
    fold_change   = (mean_Control+0.0001) / (mean_mettl3ff+0.0001),
    log2FC        = log2(fold_change),
    p_value       = wilcox.test(abundance ~ group, exact = FALSE)$p.value,
    .groups       = "drop"
  ) %>%
  arrange(p_value) %>%
  select(
    tax       = X.Taxonomy,
    mean_Control = mean_Control,
    mettl3ff= mean_mettl3ff,
    FC = fold_change,
    log2FC     = log2FC,
    pvalue        = p_value
  )

# write.csv(diff_result, "clostridium_diff.csv", row.names = FALSE)
# diff_result=read.csv("./两组差异检验统计表_DiffTwoGroup_NR_Species_Clostridium.csv")
# diff_result=diff_result[-grep("sp.",diff_result$Name),]
library(tidyverse)
library(ggrepel)

diff_plot <- diff_result %>%
  mutate(
    color_group = case_when(
      pvalue < 0.05 & log2FC > 0 ~ "Up",
      pvalue < 0.05 & log2FC < 0 ~ "Down",
      TRUE ~ "NS"
    ),
    
    up_rank = ifelse(color_group == "Up", rank(-log2FC), NA),
    down_rank = ifelse(color_group == "Down", rank(log2FC), NA),
    
    show_label = case_when(
      color_group == "Up"  ~ as.character(tax),
      color_group == "Down"  ~ as.character(tax),
      TRUE ~ NA_character_
    )
  ) %>%
  arrange(log2FC)


diff_plot=diff_plot[diff_plot$log2FC>0,]
p5=ggplot(diff_plot, aes(x = log2FC, y = fct_reorder(tax, log2FC,.desc = F))) +
  geom_col(aes(fill = color_group), width = 0.7) +
  geom_text(
    aes(label = show_label, hjust = ifelse(log2FC < 0, 1, 0)), 
    size = 8 / .pt,
    color = "#ad1002",
    na.rm = TRUE
  ) +
  scale_fill_manual(values = c("Up" = "#f46d43", "Down" = "#66c2a5", "NS" = "gray80")) +
  labs(x = "log2(Fold Change)", y = "Rank of log2(Fold Change) of species within Clostridium")+
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "inside",
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.key.size = unit(0.3,"cm"),
    legend.position.inside = c(0.9,0.1),
    plot.title = element_text(hjust = 0.5)
  )

combined_plot <-  p2 + p3 + p1 + p5+
  plot_layout(widths = c(4, 1.5, 5,10)) 
combined_plot
ggsave("./clostridium_species_screen.pdf",width = 16/2.53,height = 6/2.53)
