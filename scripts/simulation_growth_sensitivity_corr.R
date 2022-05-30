library(dplyr)
library(shellpipes)

loadEnvironments()

corrmin <- 0.1
corrmax <- 0.9
nsim <- 21

default <- data.frame(
  log_v_ini_log_v_max_corr=0.5,
  log_v_ini_a_corr=0.5,
  log_v_max_a_corr=0.5,
  log_v_ini_tau_onset_corr=-0.5,
  log_v_max_tau_onset_corr=-0.5,
  a_tau_onset_corr=-0.5,
  log_v_max_beta_max_corr=0.5
)

replace <- data.frame(
  log_v_ini_log_v_max_corr=seq(corrmin, corrmax, length.out=nsim),
  log_v_ini_a_corr=seq(corrmin, corrmax, length.out=nsim),
  log_v_max_a_corr=seq(corrmin, corrmax, length.out=nsim),
  log_v_ini_tau_onset_corr=seq(-corrmax, -corrmin, length.out=nsim),
  log_v_max_tau_onset_corr=seq(-corrmax, -corrmin, length.out=nsim),
  a_tau_onset_corr=seq(-corrmax, -corrmin, length.out=nsim),
  log_v_max_beta_max_corr=seq(corrmin, corrmax, length.out=nsim)
)

paramdata_corr <- lapply(1:ncol(replace), function(x) {
  out <- merge(default[,-x], replace[,x])
  
  names(out)[ncol(out)] <- names(default)[x]
  
  out$group <- names(default)[x]
  
  out
}) %>%
  bind_rows

simulation_growth_sensitivity_corr <- vector('list', nrow(paramdata_corr))

set.seed(seed)
for (i in 1:nrow(paramdata_corr)) {
  priorfun <- function(n) {
    do.call("singanayagam_prior", c(n=n, paramdata_corr[i,-ncol(paramdata_corr)]))
  }
  
  simulation_growth_sensitivity_corr[[i]] <- simgrowth(k=k, tmax=tmax, r=r, I0=I0, priorfun=priorfun)    
} 

saveVars(simulation_growth_sensitivity_corr, paramdata_corr)
