library(dplyr)
library(tidyr)
library(ggplot2); theme_set(theme_bw(base_family="Times"))
library(egg)
library(shellpipes)

summarize_growth_sensitivity_corr <- rdsRead()

summarize_growth_sensitivity_corr2 <- summarize_growth_sensitivity_corr %>%
  group_by(sim, type, group,
           log_v_ini_a_corr, 
           log_v_max_a_corr, log_v_ini_tau_onset_corr, log_v_max_tau_onset_corr, a_tau_onset_corr, 
           log_v_max_beta_max_corr, log_v_ini_log_v_max_corr) %>%
  summarize(
    mean_ct=mean(ct),
    lwr_ct=t.test(ct)[[4]][1],
    upr_ct=t.test(ct)[[4]][2],
    mean_a=mean(a),
    lwr_a=t.test(a)[[4]][1],
    upr_a=t.test(a)[[4]][2],
    mean_log_v_max=mean(log_v_max),
    lwr_log_v_max=t.test(log_v_max)[[4]][1],
    upr_log_v_max=t.test(log_v_max)[[4]][2]
  )

summarize_growth_sensitivity_corr3 <- summarize_growth_sensitivity_corr2 %>%
  gather(key, value, -sim, -type, -mean_ct, -lwr_ct, -upr_ct, 
         -mean_a, -lwr_a, -upr_a,
         -mean_log_v_max, -lwr_log_v_max, -upr_log_v_max, -group) %>%
  filter(
    group==key
  )

g1 <- ggplot(summarize_growth_sensitivity_corr3) +
  geom_point(aes(value, mean_ct, col=type, shape=type)) +
  geom_errorbar(aes(value, ymin=lwr_ct, ymax=upr_ct, col=type), width=0) +
  geom_smooth(aes(value, mean_ct, col=type, fill=type)) +
  scale_x_continuous("Correlation coefficient") +
  scale_y_continuous("Mean ct") +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  facet_wrap(~group, scale="free_x") +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    legend.position=c(0.5, 0.2),
    legend.title=element_blank()
  )

g2 <- ggplot(summarize_growth_sensitivity_corr3) +
  geom_point(aes(value, mean_a, col=type, shape=type)) +
  geom_errorbar(aes(value, ymin=lwr_a, ymax=upr_a, col=type), width=0) +
  geom_smooth(aes(value, mean_a, col=type, fill=type)) +
  scale_x_continuous("Correlation coefficient") +
  scale_y_continuous("Mean viral growth rate") +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  facet_wrap(~group, scale="free_x") +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    legend.position=c(0.5, 0.2),
    legend.title=element_blank()
  )

g3 <- ggplot(summarize_growth_sensitivity_corr3) +
  geom_point(aes(value, mean_log_v_max, col=type, shape=type)) +
  geom_errorbar(aes(value, ymin=lwr_log_v_max, ymax=upr_log_v_max, col=type), width=0) +
  geom_smooth(aes(value, mean_log_v_max, col=type, fill=type)) +
  scale_x_continuous("Correlation coefficient") +
  scale_y_continuous("Mean log10(peak viral load)") +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  facet_wrap(~group, scale="free_x") +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    legend.position=c(0.5, 0.2),
    legend.title=element_blank()
  )

gfinal <- ggarrange(g1, g2, g3, nrow=2, draw=FALSE)

saveGG(gfinal, width=16, height=12)
