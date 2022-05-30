library(dplyr)
library(shellpipes)
rpcall("summarize_growth_sensitivity_k.Rout summarize_growth_sensitivity_k.R simulation_growth_sensitivity_k.rds baseparam.rda R/viraltraj.rda R/util.rda")

loadEnvironments()

simulation_growth_sensitivity_k <- rdsRead()

summarize_growth_sensitivity_k <- lapply(simulation_growth_sensitivity_k, function(x) {
  d1 <- lapply(x, function(y) {
    k <- y$param$k[1]
    r <- y$param$r[1]
    y$param <- y$param[,-c(11,12)]
    
    gca <- get_ct_all(y, t_sample=tmax, return_param = TRUE, type="all")
    
    gca$k <- k
    gca$r <- r
    gca
  }) %>%
    bind_rows(.id="sim") %>% 
    mutate(type="All")
  
  d2 <- lapply(x, function(y) {
    k <- y$param$k[1]
    r <- y$param$r[1]
    y$param <- y$param[,-c(11,12)]
    
    gca <- get_ct_all(y, t_sample=tmax, return_param = TRUE, type="onset")
    
    gca$k <- k
    gca$r <- r
    gca
  }) %>%
    bind_rows(.id="sim") %>% 
    mutate(type="Symptom onset")
  
  d3 <- lapply(x, function(y) {
    k <- y$param$k[1]
    r <- y$param$r[1]
    y$param <- y$param[,-c(11,12)]
    
    gca <- get_ct_all(y, t_sample=tmax, return_param = TRUE, type="report")
    
    gca$k <- k
    gca$r <- r
    gca
  }) %>%
    bind_rows(.id="sim") %>% 
    mutate(type="Case report")
  
  dall <- bind_rows(
    d1, d2, d3
  )
}) %>%
  bind_rows

rdsSave(summarize_growth_sensitivity_k)
