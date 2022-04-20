library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family="Times"))
library(gridExtra)
library(egg)
source("../R/viraltraj.R")
source("../R/util.R")
load("../simulation/simulation_base_example.rda")

t_ct <- seq(40, 120, by=1)

daily_data <- incidence_to_daily(simulation_base_example)

ct_data_all <- lapply(t_ct, function(t) {
  get_ct_all(simulation_base_example, t, return_param = TRUE, type="all")
}) %>%
  bind_rows %>%
  mutate(
    type="All infections"
  )

ct_data_onset <- lapply(t_ct, function(t) {
  get_ct_all(simulation_base_example, t, return_param = TRUE, type="onset")
}) %>%
  bind_rows %>%
  mutate(
    type="Symptom onset"
  )

ct_data_report <- lapply(t_ct, function(t) {
  get_ct_all(simulation_base_example, t, return_param = TRUE, type="report")
}) %>%
  bind_rows %>%
  mutate(
    type="Case report"
  )

ct_data_comb <- bind_rows(
  ct_data_all,
  ct_data_onset,
  ct_data_report
)

ct_data_summary <- ct_data_comb  %>%
  group_by(type, t_sample) %>%
  summarize(
    mean_ct=mean(ct),
    mean_tinf=mean(t_sample - t_infected),
    mean_inc=mean(tau_onset),
    mean_tau_max=mean(tau_max),
    mean_a=mean(a),
    mean_log_v_max=mean(log_v_max)
  )

g1 <- ggplot(daily_data) +
  geom_point(aes(time, incidence)) +
  geom_line(aes(time, incidence)) +
  scale_x_continuous("Time (days)", limits=c(0, 140), expand=c(0, 0)) +
  scale_y_continuous("Daily new infections", expand=c(0, 0), limits=c(0, 1800)) +
  theme(
    panel.grid = element_blank()
  )

g2 <- ggplot(ct_data_summary) +
  geom_point(aes(t_sample, mean_tinf, col=type)) +
  geom_smooth(aes(t_sample, mean_tinf, col=type, fill=type)) +
  scale_x_continuous("Time (days)", limits=c(0, 140), expand=c(0, 0)) +
  scale_y_continuous("Mean time since infection (days)") +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.1, 0.8),
    legend.background = element_rect(fill=NA),
    panel.grid = element_blank()
  )

g3 <- ggplot(ct_data_summary) +
  geom_point(aes(t_sample, -mean_ct, col=type)) +
  geom_smooth(aes(t_sample, -mean_ct, col=type, fill=type)) +
  scale_x_continuous("Time (days)", limits=c(0, 140), expand=c(0, 0)) +
  scale_y_continuous("Mean Ct values (days)",
                     breaks=c(-21, -24, -27, -30, -33),
                     labels=c(21, 24, 27, 30, 33)) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme(
    legend.position = "none",
    panel.grid = element_blank()
  )

gcomb <- ggarrange(g1, g2, g3, nrow=3, draw=FALSE)

ggsave("figure_simulation_base_example.pdf", gcomb, width=8, height=8)

g4 <- ggplot(ct_data_summary) +
  geom_point(aes(t_sample, mean_a, col=type)) +
  geom_smooth(aes(t_sample, mean_a, col=type, fill=type)) +
  scale_x_continuous("Time (days)", limits=c(0, 140), expand=c(0, 0)) +
  scale_y_continuous("Mean viral load growth rate (1/day)",) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme(
    legend.position = "none",
    panel.grid = element_blank()
  )

g5 <- ggplot(ct_data_summary) +
  geom_point(aes(t_sample, mean_log_v_max, col=type)) +
  geom_smooth(aes(t_sample, mean_log_v_max, col=type, fill=type)) +
  scale_x_continuous("Time (days)", limits=c(0, 140), expand=c(0, 0)) +
  scale_y_continuous("Mean maximum log10(viral load)") +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme(
    legend.position = "none",
    panel.grid = element_blank()
  )

g6 <- ggplot(ct_data_summary) +
  geom_point(aes(t_sample, mean_tau_max, col=type)) +
  geom_smooth(aes(t_sample, mean_tau_max, col=type, fill=type)) +
  scale_x_continuous("Time (days)", limits=c(0, 140), expand=c(0, 0)) +
  scale_y_continuous("Mean time to peak viral load (days)") +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme(
    legend.position = "none",
    panel.grid = element_blank()
  )

g7 <- ggplot(ct_data_summary) +
  geom_point(aes(t_sample, mean_inc, col=type)) +
  geom_smooth(aes(t_sample, mean_inc, col=type, fill=type)) +
  scale_x_continuous("Time (days)", limits=c(0, 140), expand=c(0, 0)) +
  scale_y_continuous("Mean incubation period (days)") +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme(
    legend.position = "none",
    panel.grid = element_blank()
  )

gcomb2 <- ggarrange(g4, g5, g6, g7, nrow=2, draw=FALSE)

ggsave("figure_simulation_base_example_correlation.pdf", gcomb2, width=8, height=8)

ct_data_report2 <- ct_data_report %>%
  mutate(
    t_sample2=floor(t_sample/10) * 10
  )

g8 <- ggplot(ct_data_report2) +
  geom_smooth(aes(tau_report, ct, group=t_sample2, col=factor(t_sample2), fill=factor(t_sample2)), alpha=0.1) +
  scale_x_continuous("Time since symptom onset") +
  scale_y_reverse("Ct values") +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  facet_wrap(~t_sample2, nrow=1) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

ggsave("figure_simulation_base_example_phase.pdf", g8, width=8, height=4)
