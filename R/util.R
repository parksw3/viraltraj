##' @param n number of days in a window
incidence_to_daily <- function(simulation,
                               n=1) {
  data <- simulation$data
  
  tt <- table(floor(data$time/n))
  
  out <- data.frame(
    time=0:as.numeric(tail(names(tt), 1)),
    incidence=0
  )
  
  out$incidence[as.numeric(names(tt))+1] <- unname(tt)
  
  out$time <- out$time * n
  
  out
}

get_ct_individual <- function(param, t_infected, t_sample) {
  tau <- t_sample-t_infected
  
  do.call("singanayagam_ct", c(param, tau=tau))
}

get_ct_internal <- function(ww, 
                            param, 
                            t_infected, 
                            t_sample,
                            cutoff,
                            return_param) {
  ct <- sapply(ww, function(x) get_ct_individual(param[x,], t_infected[x], t_sample))
  
  if (return_param) {
    param_filter <- param[ww,]
    
    param_filter$ct <- ct
    
    param_filter$t_sample <- t_sample
    param_filter$t_infected <- t_infected[ww]
    
    param_filter[ct < cutoff,]
  } else {
    ct[ct < cutoff]
  }
}

get_ct_all <- function(simulation, 
                       t_sample,
                       cutoff=40,
                       taumax=20,
                       return_param=FALSE,
                       type=c("all", "onset", "report")) {
  type <- match.arg(type)
  
  param <- simulation$param
  t_infected <- simulation$t_infected
  
  if (type=="all") {
    ww <- which(t_sample - taumax < t_infected & t_infected <= t_sample)
  } else if (type=="onset") {
    t_report <- t_infected + param$tau_onset
    
    ww <- which(t_sample-1 < t_report & t_report <= t_sample)
  } else if (type=="report") {
    t_report <- t_infected + param$tau_onset+param$tau_report
    
    ww <- which(t_sample-1 < t_report & t_report <= t_sample)
  }
  
  if (length(ww) > 0) {
    get_ct_internal(ww, param, t_infected, t_sample, cutoff, return_param)
  } else {
    NULL
  }
}
