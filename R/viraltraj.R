singanayagam_tau_max <- function(a,
                                 b,
                                 log_v_max,
                                 log_v_ini,
                                   ...) {
  uniroot(function(x) {
    10^(log_v_max-log_v_ini) * (a+b) - b * exp(a * x) - a * exp(- b* x)
  }, c(0, 50))[[1]]
}

singanayagam_vt <- function(tau,
                            a,
                            b,
                            log_v_max,
                            tau_max,
                            ...) {
  10^log_v_max * (a+b)/(b * exp(-a * (tau - tau_max)) + a * exp(b*(tau - tau_max)))
}

singanayagam_vt_to_ct <- function(x,
                                  kappa=3.2,
                                  alpha=1.4,
                               ...) {
  40 - (log(x) - kappa) * alpha
}

singanayagam_ct <- function(tau,
                               a,
                               b,
                               log_v_max,
                               tau_max,
                               kappa=3.2,
                               alpha=1.4,
                               ...) {
  singanayagam_vt_to_ct(singanayagam_vt(tau, a, b, log_v_max, tau_max), kappa, alpha)
}

singanayagam_vt_to_beta <- function(x, 
                                    theta,
                                    beta_max,
                                    tau_onset,
                                    symptom,
                                    ...) {
  beta <- plogis(log(x)-theta)*beta_max
  
  if (symptom) {
    
  }
  
  beta
}

singanayagam_beta <- function(tau,
                                 a,
                                 b,
                                 log_v_max,
                                 tau_max,
                                 theta,
                                 beta_max,
                                 tau_onset,
                                 symptom,
                                 ...) {
  singanayagam_vt_to_beta(singanayagam_vt(tau, a, b, log_v_max, tau_max), theta, beta_max, tau_onset, symptom)
}

singanayagam_prior <- function(n=1000, 
                               a_logmean=log(6),
                               a_logsd=0.35,
                               b_logmean=log(2),
                               b_logsd=0.2,
                               log_v_max_mean=8,
                               log_v_max_sd=0.5,
                               log_v_ini_mean=-4,
                               log_v_ini_sd=0.5,
                               theta_logmean=log(13),
                               theta_logsd=0.1,
                               beta_max_logmean=log(0.5),
                               beta_max_logsd=0.2,
                               tau_onset_logmean=log(4),
                               tau_onset_logsd=0.5,
                               tau_report_logmean=log(3),
                               tau_report_logsd=0.5,
                               p_s=0.5,
                               log_v_ini_log_v_max_corr=0.5,
                               log_v_ini_a_corr=0.5,
                               log_v_max_a_corr=0.5,
                               log_v_ini_tau_onset_corr=-0.5,
                               log_v_max_tau_onset_corr=-0.5,
                               a_tau_onset_corr=-0.5,
                               log_v_max_beta_max_corr=0.5,
                               seed) {
  if (!missing(seed)) set.seed(seed)
  
  ## higher initial viral load correlates with higher peak viral load
  ## higher initial viral load correlates with faster viral growth
  ## higher peak viral load correlates with faster viral growth
  ## higher initial viral load correlates with earlier symptom onset (negative correlation)
  ## higher peak viral load correlates with earlier symptom onset (negative correlation)
  ## faster viral growth correlates with earlier symptom onset (negative correlation)
  ## higher peak viral load correlations with higher transmission rate
  
  sigma <- matrix(
    c(
      log_v_ini_sd^2, log_v_ini_log_v_max_corr * log_v_ini_sd * log_v_max_sd, log_v_ini_a_corr * log_v_ini_sd * a_logsd, log_v_ini_tau_onset_corr * log_v_ini_sd * tau_onset_logsd, 0, 
      log_v_ini_log_v_max_corr * log_v_ini_sd * log_v_max_sd, log_v_max_sd^2, log_v_max_a_corr * log_v_max_sd * a_logsd, log_v_max_tau_onset_corr * log_v_max_sd * tau_onset_logsd, log_v_max_beta_max_corr * log_v_max_sd * beta_max_logsd,
      log_v_ini_a_corr * log_v_ini_sd * a_logsd, log_v_max_a_corr * log_v_max_sd * a_logsd, a_logsd^2, a_tau_onset_corr * a_logsd * tau_onset_logsd, 0,
      log_v_ini_tau_onset_corr * log_v_ini_sd * tau_onset_logsd, log_v_max_tau_onset_corr * log_v_max_sd * tau_onset_logsd, a_tau_onset_corr * a_logsd * tau_onset_logsd, tau_onset_logsd^2, 0,
      0, log_v_max_beta_max_corr * log_v_max_sd * beta_max_logsd, 0, 0, beta_max_logsd^2
    ),
    5, 5
  )
  
  parset <- mvtnorm::rmvnorm(n, mean=c(log_v_ini_mean, log_v_max_mean, a_logmean, tau_onset_logmean, beta_max_logmean), sigma=sigma)
  
  out <- data.frame(
    a=exp(parset[,3]),
    b=rlnorm(n, b_logmean, b_logsd),
    log_v_ini=parset[,1],
    log_v_max=parset[,2],
    theta=rlnorm(n, theta_logmean, theta_logsd),
    beta_max=exp(parset[,5]),
    tau_onset=exp(parset[,4]),
    tau_report=rlnorm(n, tau_report_logmean, tau_report_logsd),
    symptom=rbinom(n, 1, p_s)
  )
  
  out$tau_max <- apply(out, 1, function(x) do.call("singanayagam_tau_max", as.list(x)))
  
  out
}
