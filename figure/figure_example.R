library(dplyr)
library(ggplot2); theme_set(theme_bw())
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

paramdata <- apply(param, 1, function(x) {
  vt <- do.call(singanayagam_vt, c(list(tau=tauvec), as.list(x)))
  ct <- singanayagam_vt_to_ct(vt)
  beta <- singanayagam_vt_to_beta(vt, theta=x[["theta"]], beta_max=x[["beta_max"]], symptom=x[["symptom"]])
  
  nu <- sum(beta) * taudiff
  duration_ct <- sum(ct < 40) * 0.1
  
  c(nu=nu, duration_ct=duration_ct, x)
}, simplify = FALSE) %>%
  bind_rows

g1 <- ggplot(trajdata2) +
  geom_line(aes(time, vt, group=sim), alpha=0.1) +
  scale_x_continuous("Time since infection (days)") +
  scale_y_log10("Viral load") +
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

g3 <- ggplot(trajdata2) +
  geom_line(aes(time, beta, group=sim), alpha=0.1) +
  scale_x_continuous("Time since infection (days)") +
  scale_y_continuous("Infection kernal (1/day)") +
  theme(
    panel.grid = element_blank()
  )

g4 <- ggplot(paramdata) +
  geom_point(aes(log_v_max, nu)) +
  scale_x_continuous("log10(peak viral load)") +
  scale_y_continuous("Individual reproduction number") +
  theme(
    panel.grid = element_blank()
  )

g5 <- ggplot(paramdata) +
  geom_point(aes(duration_ct, nu)) +
  scale_x_continuous("Duration of viral detection, ct < 40, (days)") +
  scale_y_continuous("Individual reproduction number") +
  theme(
    panel.grid = element_blank()
  )

g6 <- ggplot(paramdata) +
  geom_point(aes(a, tau_onset)) +
  scale_x_continuous("Viral load growth rate (1/days)") +
  scale_y_continuous("Symptom onset (days)") +
  theme(
    panel.grid = element_blank()
  )

g7 <- ggplot(paramdata) +
  geom_point(aes(a, tau_onset)) +
  scale_x_continuous("Viral load growth rate (1/days)") +
  scale_y_continuous("Timing of symptom onset (days)") +
  theme(
    panel.grid = element_blank()
  )

g8 <- ggplot(paramdata) +
  geom_point(aes(tau_max, tau_onset)) +
  scale_x_continuous("Timing of peak viral load (days)") +
  scale_y_continuous("Timing of symptom onset (days)") +
  theme(
    panel.grid = element_blank()
  )

gcomb1 <- ggarrange(g1, g2, g3, nrow=1, draw=FALSE)
gcomb2 <- ggarrange(g4, g5, g6, g7, nrow=2, draw=FALSE)

gfinal <- arrangeGrob(gcomb1, gcomb2, nrow=2, heights=c(1, 2))

ggsave("figure_example.pdf", gfinal, width=10, height=8)
