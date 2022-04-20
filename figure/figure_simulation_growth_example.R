library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family="Times"))
library(gridExtra)
library(egg)
library(viridis)
source("../R/viraltraj.R")
source("../R/simgrowth.R")
source("../R/util.R")

seed <- 101
tmax <- 40
r <- 0.1
I0 <- 10

sim1 <- simgrowth(k=1, seed=seed, tmax=tmax, r=r, I0=I0)
sim2 <- simgrowth(k=10, seed=seed, tmax=tmax, r=r, I0=I0)
sim3 <- simgrowth(k=100, seed=seed, tmax=tmax, r=r, I0=I0)

daily_data1 <- incidence_to_daily(sim1)

daily_data2 <- incidence_to_daily(sim2)

daily_data3 <- incidence_to_daily(sim3)

g1 <- ggplot(daily_data1) +
  stat_function(fun= function(x) I0 * exp(r * x), lwd=2, col="gray") +
  geom_point(aes(time, incidence), shape=1) +
  scale_x_continuous("Time (days)") +
  scale_y_continuous("Daily incidence", limits=c(0, 900)) +
  theme(
    panel.grid = element_blank()
  )

g2 <- ggplot(daily_data2) +
  stat_function(fun= function(x) I0 * exp(r * x), lwd=2, col="gray") +
  geom_point(aes(time, incidence), shape=1) +
  scale_x_continuous("Time (days)") +
  scale_y_continuous("Daily incidence", limits=c(0, 900)) +
  theme(
    panel.grid = element_blank()
  )

g3 <- ggplot(daily_data3) +
  stat_function(fun= function(x) I0 * exp(r * x), lwd=2, col="gray") +
  geom_point(aes(time, incidence), shape=1) +
  scale_x_continuous("Time (days)") +
  scale_y_continuous("Daily incidence", limits=c(0, 900)) +
  theme(
    panel.grid = element_blank()
  )

ct_data_all1 <- get_ct_all(sim1, tmax, return_param = TRUE, type="all") %>% mutate(type="All")
ct_data_all2 <- get_ct_all(sim2, tmax, return_param = TRUE, type="all") %>% mutate(type="All")
ct_data_all3 <- get_ct_all(sim3, tmax, return_param = TRUE, type="all") %>% mutate(type="All")

ct_data_onset1 <- get_ct_all(sim1, tmax, return_param = TRUE, type="onset") %>% mutate(type="Symptom onset")
ct_data_onset2 <- get_ct_all(sim2, tmax, return_param = TRUE, type="onset") %>% mutate(type="Symptom onset")
ct_data_onset3 <- get_ct_all(sim3, tmax, return_param = TRUE, type="onset") %>% mutate(type="Symptom onset")

ct_data_report1 <- get_ct_all(sim1, tmax, return_param = TRUE, type="report") %>% mutate(type="Case report")
ct_data_report2 <- get_ct_all(sim2, tmax, return_param = TRUE, type="report") %>% mutate(type="Case report")
ct_data_report3 <- get_ct_all(sim3, tmax, return_param = TRUE, type="report") %>% mutate(type="Case report")

ct_data1 <- bind_rows(
  ct_data_all1,
  ct_data_onset1,
  ct_data_report1
)

ct_data2 <- bind_rows(
  ct_data_all2,
  ct_data_onset2,
  ct_data_report2
)

ct_data3 <- bind_rows(
  ct_data_all3,
  ct_data_onset3,
  ct_data_report3
)

g4 <- ggplot(ct_data1) +
  geom_bar(aes(t_sample - t_infected, y = (..prop..), fill=type, group=type), position=position_dodge(1, preserve = "single")) +
  scale_x_continuous("Time since infection (days)", limits=c(0, 19), expand=c(0, 0))  +
  scale_y_continuous("Frequency", limits=c(0, 0.35), expand=c(0, 0)) +
  scale_fill_viridis_d() +
  theme(
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.position = c(0.75, 0.8),
    legend.background = element_rect(fill=NA)
  )

g5 <- ggplot(ct_data2) +
  geom_bar(aes(t_sample - t_infected, y = (..prop..), fill=type, group=type), position=position_dodge(1, preserve = "single")) +
  scale_x_continuous("Time since infection (days)", limits=c(0, 19), expand=c(0, 0))  +
  scale_y_continuous("Frequency", limits=c(0, 0.35), expand=c(0, 0)) +
  scale_fill_viridis_d() +
  theme(
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  )

g6 <- ggplot(ct_data3) +
  geom_bar(aes(t_sample - t_infected, y = (..prop..), fill=type, group=type), position=position_dodge(1, preserve = "single")) +
  scale_x_continuous("Time since infection (days)", limits=c(0, 19), expand=c(0, 0))  +
  scale_y_continuous("Frequency", limits=c(0, 0.35), expand=c(0, 0)) +
  scale_fill_viridis_d() +
  theme(
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  )

g7 <- ggplot(ct_data1) +
  geom_density(aes(ct, group=type, fill=type, col=type, lty=type), alpha=0.3) +
  scale_x_continuous("Ct value", limits=c(14, 40), expand=c(0, 0))  +
  scale_y_continuous("Density", limits=c(0, 0.095), expand=c(0, 0)) +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  theme(
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  )

g8 <- ggplot(ct_data2) +
  geom_density(aes(ct, group=type, fill=type, col=type, lty=type), alpha=0.3) +
  scale_x_continuous("Ct value", limits=c(14, 40), expand=c(0, 0))  +
  scale_y_continuous("Density", limits=c(0, 0.095), expand=c(0, 0)) +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  theme(
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  )

g9 <- ggplot(ct_data3) +
  geom_density(aes(ct, group=type, fill=type, col=type, lty=type), alpha=0.3) +
  scale_x_continuous("Ct value", limits=c(14, 40), expand=c(0, 0))  +
  scale_y_continuous("Density", limits=c(0, 0.095), expand=c(0, 0)) +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  theme(
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  )

g10 <- ggplot(filter(ct_data1, type != "Symptom onset")) +
  geom_smooth(aes(t_sample - tau_onset - t_infected, ct, group=type, col=type, fill=type)) +
  scale_x_continuous("Time since symptom onset (days)", limits=c(-18, 18), expand=c(0, 0)) +
  scale_y_reverse("Ct value", limits=c(45, 14), expand=c(0, -2)) +
  scale_fill_manual(values=viridis(3)[1:2]) +
  scale_color_manual(values=viridis(3)[1:2]) +
  theme(
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  )
  
g11 <- ggplot(filter(ct_data2, type != "Symptom onset")) +
  geom_smooth(aes(t_sample - tau_onset - t_infected, ct, group=type, col=type, fill=type)) +
  scale_x_continuous("Time since symptom onset (days)", limits=c(-18, 18), expand=c(0, 0)) +
  scale_y_reverse("Ct value", limits=c(45, 14), expand=c(0, -2)) +
  scale_fill_manual(values=viridis(3)[1:2]) +
  scale_color_manual(values=viridis(3)[1:2]) +
  theme(
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  )

g12 <- ggplot(filter(ct_data3, type != "Symptom onset")) +
  geom_smooth(aes(t_sample - tau_onset - t_infected, ct, group=type, col=type, fill=type)) +
  scale_x_continuous("Time since symptom onset (days)", limits=c(-18, 18), expand=c(0, 0)) +
  scale_y_reverse("Ct value", limits=c(45, 14), expand=c(0, -2)) +
  scale_fill_manual(values=viridis(3)[1:2]) +
  scale_color_manual(values=viridis(3)[1:2]) +
  theme(
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  )

gfinal <- ggarrange(g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11, g12, nrow=3, byrow=FALSE, draw=FALSE)

ggsave("figure_simulation_growth_example.pdf", gfinal, width=12, height=8)
