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
                               tau_max_logmean=log(6),
                               tau_max_logsd=0.2,
                               theta_logmean=log(13),
                               theta_logsd=0.1,
                               beta_max_logmean=log(0.5),
                               beta_max_logsd=0.2,
                               beta_max_log_v_max_corr=0.5,
                               tau_onset_logmean=log(4),
                               tau_onset_logsd=0.5,
                               tau_max_tau_onset_corr=0.5,
                               a_tau_onset_corr=-0.5,
                               tau_report_logmean=log(3),
                               tau_report_logsd=0.5,
                               p_s=0.5,
                               seed) {
  if (!missing(seed)) set.seed(seed)
  
  sigma1 <- matrix(
    c(tau_max_logsd^2, tau_max_tau_onset_corr * tau_max_logsd * tau_onset_logsd, 0,
      tau_max_tau_onset_corr * tau_max_logsd * tau_onset_logsd, tau_onset_logsd^2, a_tau_onset_corr * a_logsd * tau_onset_logsd,
      0, a_tau_onset_corr * a_logsd * tau_onset_logsd, a_logsd^2),
    3, 3
  )
  
  tauset <- exp(mvtnorm::rmvnorm(n, mean=c(tau_max_logmean, tau_onset_logmean, a_logmean), sigma=sigma1))
  
  sigma2 <- matrix(
    c(beta_max_logsd^2, beta_max_log_v_max_corr * beta_max_logsd * log_v_max_sd,
      beta_max_log_v_max_corr * beta_max_logsd * log_v_max_sd, log_v_max_sd^2),
    2, 2
  )
  
  betaset <- mvtnorm::rmvnorm(n, mean=c(beta_max_logmean, log_v_max_mean), sigma=sigma2)
  
  data.frame(
    a=tauset[,3],
    b=rlnorm(n, b_logmean, b_logsd),
    log_v_max=betaset[,2],
    tau_max=tauset[,1],
    theta=rlnorm(n, theta_logmean, theta_logsd),
    beta_max=exp(betaset[,1]),
    tau_onset=tauset[,2],
    tau_report=rlnorm(n, tau_report_logmean, tau_report_logsd),
    symptom=rbinom(n, 1, p_s)
  )
}
