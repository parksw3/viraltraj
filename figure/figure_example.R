library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(gridExtra)
library(egg)
source("../R/viraltraj.R")

nsim <- 1000
taudiff <- 0.1
taumax <- 20
tauvec <- seq(0, taumax, by=taudiff)

param <- singanayagam_prior(nsim, seed=101)

trajdata <- apply(param, 1, function(x) {
  vt <- do.call(singanayagam_vt, c(list(tau=tauvec), as.list(x)))
  ct <- singanayagam_vt_to_ct(vt)
  beta <- singanayagam_vt_to_beta(vt, theta=x[["theta"]], beta_max=x[["beta_max"]], symptom=x[["symptom"]])
  
  data.frame(
    time=tauvec,
    vt=vt,
    ct=ct,
    beta=beta
  )
}, simplify = FALSE) %>%
  bind_rows(.id="sim")

trajdata2 <- trajdata %>%
  filter(
    sim %in% sample(1:nsim, 100)
  ) %>%
  mutate(
    ct=ifelse(ct > 40, 40, ct)
  )

g1 <- ggplot(trajdata2) +
  geom_line(aes(time, vt, group=sim), alpha=0.1) +
  scale_x_continuous("Time since infection (days)") +
  scale_y_log10("Viral load (viral copies per ml)", limits=c(1e-5, NA)) +
  theme(
    panel.grid = element_blank()
  )

g2 <- ggplot(trajdata2) +
  geom_line(aes(time, -ct, group=sim), alpha=0.1) +
  scale_x_continuous("Time since infection (days)") +
  scale_y_continuous("Ct value",
                     breaks=-c(40, 35, 30, 25, 20),
                     labels=c(40, 35, 30, 25, 20)) +
  theme(
    panel.grid = element_blank()
  )

g3 <- ggplot(param) +
  geom_point(aes(log_v_ini, log_v_max)) +
  scale_x_continuous("log10(initial viral load)") +
  scale_y_continuous("log10(peak viral load)") +
  theme(
    panel.grid = element_blank()
  )

g4 <- ggplot(param) +
  geom_point(aes(log_v_ini, a)) +
  scale_x_continuous("log10(initial viral load)") +
  scale_y_continuous("Viral load growth rate (1/days)") +
  theme(
    panel.grid = element_blank()
  )

g5 <- ggplot(param) +
  geom_point(aes(a, log_v_max)) +
  scale_x_continuous("Viral load growth rate (1/days)") +
  scale_y_continuous("log10(peak viral load)") +
  theme(
    panel.grid = element_blank()
  )

g6 <- ggplot(param) +
  geom_point(aes(log_v_max, tau_max)) +
  scale_x_continuous("log10(peak viral load)") +
  scale_y_continuous("Timing of peak viral load (days)") +
  theme(
    panel.grid = element_blank()
  )

gfinal <- ggarrange(g1, g2, g3, g4, g5, g6, nrow=3, draw=FALSE)

ggsave("figure_example.pdf", gfinal, width=10, height=8)
