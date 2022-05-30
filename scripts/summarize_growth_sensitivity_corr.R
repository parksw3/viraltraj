library(dplyr)
library(shellpipes)
rpcall("summarize_growth_sensitivity_corr.Rout summarize_growth_sensitivity_corr.R simulation_growth_sensitivity_corr.rda baseparam.rda R/viraltraj.rda R/util.rda")

loadEnvironments()

summarize_growth_sensitivity_corr <- lapply(simulation_growth_sensitivity_corr, function(x) {
  x$param <- x$param[,-c(11,12)]
  
  d1 <- get_ct_all(x, t_sample=tmax, return_param = TRUE, type="all") %>% 
    mutate(type="All")
  
  d2 <- get_ct_all(x, t_sample=tmax, return_param = TRUE, type="onset") %>% 
    mutate(type="Symptom onset")
  
  d3 <- get_ct_all(x, t_sample=tmax, return_param = TRUE, type="report") %>% 
    mutate(type="Case report")
  
  dall <- bind_rows(
    d1, d2, d3
  )
}) %>%
  bind_rows(.id="sim")

paramdata_corr$sim <- 1:nrow(paramdata_corr)

summarize_growth_sensitivity_corr <- merge(summarize_growth_sensitivity_corr, paramdata_corr)

rdsSave(summarize_growth_sensitivity_corr)
