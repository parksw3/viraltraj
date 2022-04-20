library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family="Times"))
library(gridExtra)
library(ggpubr)
library(viridis)
source("../R/viraltraj.R")
source("../R/simgrowth.R")
source("../R/util.R")

tau_onset_logmean <- log(4)
tau_onset_logsd <- 0.5
tau_onset_mean <- exp(tau_onset_logmean + tau_onset_logsd^2/2)

seed <- 101
tmax <- 40
rvec <- seq(-0.1, 0.2, by=0.01)

simall <- lapply(rvec, function(r) {
  print(r)
  
  if (r==0) {
    I0 <- 20000/(tmax+1)
  } else {
    I0  <- 20000/(exp(r * tmax)/r - 1/r) 
  }
  
  
  out <- simgrowth(k=100, seed=seed, tmax=tmax, r=r, I0=I0)
  
  out$r <- r
  
  out
})

ctall <- lapply(simall, function(x) {
  d1 <- get_ct_all(x, tmax, return_param = TRUE, type="all") %>% mutate(type="All")
  d2 <- get_ct_all(x, tmax, return_param = TRUE, type="onset") %>% mutate(type="Symptom onset")
  d3 <- get_ct_all(x, tmax, return_param = TRUE, type="report") %>% mutate(type="Case report")
  
  dall <- bind_rows(d1, d2, d3)
  
  dall$r <- x$r
  
  dall
}) %>%
  bind_rows

ctall2 <- ctall %>%
  group_by(r, type) %>%
  summarize(
    ct_mean=mean(ct),
    ct_lwr=t.test(ct)[[4]][1],
    ct_upr=t.test(ct)[[4]][2],
    a_mean=mean(a),
    a_lwr=t.test(a)[[4]][1],
    a_upr=t.test(a)[[4]][2],
    log_v_max_mean=mean(log_v_max),
    log_v_max_lwr=t.test(log_v_max)[[4]][1],
    log_v_max_upr=t.test(log_v_max)[[4]][2],
    tau_onset_mean=mean(tau_onset),
    tau_onset_lwr=t.test(tau_onset)[[4]][1],
    tau_onset_upr=t.test(tau_onset)[[4]][2]
  )

g1 <- ggplot(ctall2) +
  geom_errorbar(aes(r, ymin=ct_lwr, ymax=ct_upr, col=type), width=0, position = position_dodge(0.01)) +
  geom_point(aes(r, ct_mean, col=type, shape=type), position = position_dodge(0.01)) +
  scale_x_continuous("Epidemic growth rate (1/day)") +
  scale_y_continuous("Mean Ct value") +
  scale_color_viridis_d() +
  theme(
    panel.grid = element_blank(),
    legend.title = element_blank()
  )

g2 <- ggplot(ctall2) +
  geom_errorbar(aes(r, ymin=a_lwr, ymax=a_upr, col=type), width=0, position = position_dodge(0.01)) +
  geom_point(aes(r, a_mean, col=type, shape=type), position = position_dodge(0.01)) +
  scale_x_continuous("Epidemic growth rate (1/day)") +
  scale_y_continuous("Viral load growth rate (1/day)") +
  scale_color_viridis_d() +
  theme(
    panel.grid = element_blank(),
    legend.title = element_blank()
  )

g3 <- ggplot(ctall2) +
  geom_errorbar(aes(r, ymin=log_v_max_lwr, ymax=log_v_max_upr, col=type), width=0, position = position_dodge(0.01)) +
  geom_point(aes(r, log_v_max_mean, col=type, shape=type), position = position_dodge(0.01)) +
  scale_x_continuous("Epidemic growth rate (1/day)") +
  scale_y_continuous("log10(peak viral load)") +
  scale_color_viridis_d() +
  theme(
    panel.grid = element_blank(),
    legend.title = element_blank()
  )

g4 <- ggplot(ctall2) +
  geom_hline(yintercept=tau_onset_mean, lty=2) +
  geom_errorbar(aes(r, ymin=tau_onset_lwr, ymax=tau_onset_upr, col=type), width=0, position = position_dodge(0.01)) +
  geom_point(aes(r, tau_onset_mean, col=type, shape=type), position = position_dodge(0.01)) +
  scale_x_continuous("Epidemic growth rate (1/day)") +
  scale_y_continuous("Mean incubation period (days)") +
  scale_color_viridis_d() +
  theme(
    panel.grid = element_blank(),
    legend.title = element_blank()
  )

gfinal <- ggpubr::ggarrange(g1, g2, g3, g4, nrow=2, ncol=2, common.legend = TRUE)

ggsave("figure_simulation_growth_r.pdf", gfinal, width=8, height=6)
ggsave("figure_simulation_growth_r.png", gfinal, width=8, height=6)

g5 <- ggplot(filter(ctall, type != "Symptom onset")) +
  geom_smooth(aes(t_sample - tau_onset - t_infected, ct, group=r, col=r, fill=r), se=FALSE, alpha=0.5) +
  scale_x_continuous("Time since symptom onset (days)", limits=c(-18, 18), expand=c(0, 0)) +
  scale_y_reverse("Ct value", limits=c(45, 14), expand=c(0, -2)) +
  scale_fill_viridis_c() +
  scale_color_viridis_c() +
  facet_wrap(~type) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )

ggsave("figure_simulation_growth_VL.pdf", g5, width=10, height=5)
ggsave("figure_simulation_growth_VL.png", g5, width=10, height=5)
