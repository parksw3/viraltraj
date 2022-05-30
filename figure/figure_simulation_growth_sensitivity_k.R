library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family="Times"))
library(egg)
library(shellpipes)

summarize_growth_sensitivity_k <- rdsRead()

summarize_growth_sensitivity_k2 <- summarize_growth_sensitivity_k %>%
  group_by(sim, k, type) %>%
  summarize(
    mean=mean(ct),
    median=median(ct),
    skewness=mean((ct-mean)^3)/mean((ct-mean)^2)^(3/2)
  )

g1 <- ggplot(summarize_growth_sensitivity_k2) +
  geom_point(aes(k, mean, col=type, shape=type)) +
  geom_smooth(aes(k, mean, col=type, fill=type)) +
  scale_x_log10("Overdispersion in incidence, k") +
  scale_y_continuous("Mean") +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

g2 <- ggplot(summarize_growth_sensitivity_k2) +
  geom_point(aes(k, median, col=type, shape=type)) +
  geom_smooth(aes(k, median, col=type, fill=type)) +
  scale_x_log10("Overdispersion in incidence, k") +
  scale_y_continuous("Median") +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

g3 <- ggplot(summarize_growth_sensitivity_k2) +
  geom_point(aes(k, skewness, col=type, shape=type)) +
  geom_smooth(aes(k, skewness, col=type, fill=type)) +
  scale_x_log10("Overdispersion in incidence, k") +
  scale_y_continuous("Skewness") +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_blank()
  )

gfinal <- ggarrange(g1, g2, g3, nrow=1, draw=FALSE)

saveGG(gfinal, width=8, height=3)
