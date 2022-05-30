library(shellpipes)

loadEnvironments()

kvec <- exp(seq(log(1), log(100), length.out=11))
nsim <- 10

simulation_growth_sensitivity_k <- vector('list', length(kvec))

set.seed(seed)
for (i in 1:length(kvec)) {
  simlist <- vector('list', nsim)
  
  for (j in 1:nsim) {
    simlist[[j]] <- simgrowth(k=kvec[i], tmax=tmax, r=r, I0=I0)    
  }
  
  simulation_growth_sensitivity_k[[i]] <- simlist
} 

rdsSave(simulation_growth_sensitivity_k)
