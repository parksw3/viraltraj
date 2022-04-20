##' @param x a vector of nodes
##' @param size number of nodes to pick at random
sample2 <- function(x, size) {
  if(length(x)==1) {
    rep(x, size)
  } else{
    sample(x, size, replace=TRUE)
  }
}

## creating a full graph takes up a lot of memory
## so here's an optimized version that doesn't rely on igraph but does the same thing
##' @param n number of individuals
##' @param I0 number of initial infections
simepi <- function(n=50000,
                   I0=10,
                   priorfun=singanayagam_prior,
                   taumax=20,
                   taudiff=0.1,
                   seed = NULL,
                   imax,
                   tmax=300){
  if (!is.null(seed)) set.seed(seed)
  
  V <- 1:n
  
  initial_infected <- 1:I0
  
  if (missing(imax)) imax <- n
  
  param <- priorfun(n)
  
  queue_v <- queue_t <- queue_infector <- rep(NA, I0)
  
  queue_v[1:I0] <- initial_infected 
  queue_t[1:I0] <- 0
  
  t_infected <- rep(NA, n)
  t_infected[initial_infected] <- 0
  
  t_gillespie <- NULL
  c_infected <- 0
  
  done <- rep(FALSE, n)
  infected_by <- rep(NA, n)
  
  tauvec <- seq(0, taumax, by=taudiff)
  taulength <- length(tauvec)
  
  stop <- FALSE
  
  while (!stop) {
    j.index <- which.min(queue_t)
    j <- queue_v[j.index]
    
    c_infected <- c_infected +1
    
    infected_by[j] <- queue_infector[j.index]
    t_infected[j] <- queue_t[j.index]
    
    t <- queue_t[j.index]; t_gillespie <- c(t_gillespie, t)
    
    pp <- param[j,]
    
    beta <- do.call(singanayagam_beta, c(param[j,], list(tau=tauvec)))
    
    contactset <- V[V != j]
    
    ncontact <- rbinom(taulength, 1, 1-exp(-beta*taudiff))
    
    if (sum(ncontact) > 0) {
      queue_v <- c(queue_v, sample2(contactset, sum(ncontact)))
      queue_infector <- c(queue_infector, rep(j, sum(ncontact)))
      
      queue_t <- c(queue_t, t + tauvec[which(ncontact != 0)])
    }  
    
    done[j] <- TRUE
    
    filter2 <- !done[queue_v]
    queue_v <- queue_v[filter2]
    queue_infector <- queue_infector[filter2]
    queue_t <- queue_t[filter2]
    
    stop <- (c_infected == length(V) || all(done[queue_v]) || c_infected == imax || t > tmax)
  }
  
  return(
    list(
      data=data.frame(
        time=t_gillespie[(I0):c_infected],
        infected=(I0):c_infected
      ),
      t_infected=t_infected,
      infected_by=infected_by,
      param=param
    )
  )
}
