setwd("path")
library(drc)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
data <- read.csv("./CETSA.csv")

data_long <- data %>%
  pivot_longer(cols = -Temp, names_to = "Sample", values_to = "Signal") %>%
  mutate(Group = case_when(
    str_detect(Sample, "C") ~ "DMSO",
    str_detect(Sample, "G") ~ "GLCA"
  ))

model <- drm(
  Signal ~ Temp, 
  data = data_long, 
  fct = LL.4(names = c("HillSlope", "Bottom", "Top", "Tm")),
  curveid = Group,
  robust = "median"
)

summary(model)

newdate <- expand.grid(
  Temp = seq(37, 67, 0.5),
  Group = c("DMSO", "GLCA")
)
newdate$pred <- predict(model, newdata = newdate)

p1=ggplot(data_long, aes(x = Temp, y = Signal, color = Group, fill = Group)) +
  geom_point(size = 1, alpha = 0.8) +
  geom_line(data = newdate, aes(y = pred), linewidth = 0.75) +
  scale_color_manual(values = c("grey","#4393C4")) +
  ylim(0, 130) +
  labs(
    x = "Temperature (°C)",
    y = "Relative intensity\n GPRC5A(%)",
    title = ""
  ) +
  theme_bw(base_size = 6) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),legend.key.size = unit(0.3,"cm"),
        legend.position = "inside",legend.background = element_blank(),
        legend.position.inside = c(0.15,0.3),
        panel.grid = element_blank())
p1
ggsave("./CETSA_CURVE.pdf",plot = p1,width = 5.9/2.53,height = 4/2.53)