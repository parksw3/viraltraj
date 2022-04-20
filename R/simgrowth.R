simgrowth <- function(r=0.1,
                      I0=10,
                      priorfun=singanayagam_prior,
                      k=10,
                      seed = NULL,
                      tmax=40){
  if (!is.null(seed)) set.seed(seed)
  
  incidence <- rnbinom(tmax+1, mu=I0 * exp(r * 0:tmax), size=k)
  
  t_infected <- rep(0:tmax, incidence)
  
  infected_by <- rep(NA, length(t_infected))
    
  param <- priorfun(length(t_infected))
  
  return(
    list(
      data=data.frame(
        time=t_infected,
        infected=1:length(t_infected)
      ),
      t_infected=t_infected,
      infected_by=infected_by,
      param=param
    )
  )
}
