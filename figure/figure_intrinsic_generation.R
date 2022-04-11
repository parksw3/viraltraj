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

trajdata3 <- trajdata %>%
  group_by(time) %>%
  summarize(
    sum_beta=sum(beta)
  ) %>%
  mutate(
    intrinsic=sum_beta/sum(sum_beta)/taudiff
  )

trajdata_summ <- trajdata3 %>%
  summarize(
    mean=sum(time*intrinsic)/sum(intrinsic),
    var=sum((time-mean)^2*intrinsic)/sum(intrinsic),
    cv2=var/mean^2
  )

g1 <- ggplot(trajdata2) +
  geom_line(aes(time, beta, group=sim), alpha=0.1) +
  scale_x_continuous("Time since infection (days)") +
  scale_y_continuous("Infection kernal (1/day)") +
  theme(
    panel.grid = element_blank()
  )

g2 <- ggplot(trajdata3) +
  geom_line(aes(time, intrinsic)) +
  stat_function(geom="line", fun=function(x) dgamma(x, shape=1/trajdata_summ$cv2, rate=1/trajdata_summ$cv2/1/trajdata_summ$mean), lty=2, n=100) +
  scale_x_continuous("Intrinsic generation interval (days)") +
  scale_y_continuous("Density (1/day)") +
  theme(
    panel.grid = element_blank()
  )

gcomb <- ggarrange(g1, g2, nrow=1, draw=FALSE)

ggsave("figure_intrinsic_generation.pdf", gcomb, width=10, height=5)

