#' Update nodes sociability parameters given the rest, using HMC as
#' described in Caron and Fox, Appendix F.1 a
#'
#' @param w: vector of positive scalars, nodes sociability (previous sample)
#' @param w_rem: positive scalar; total mass of discarded nodes (previous sample)
#' @param sigma: real in (-Inf, 1) (previous sample or fixed)
#' @param tau: positive scalar (previous sample or fixed)
#' @param N_counts: vector of positive integers, nodes degrees (previous sample, or observed)
#' @param L: HMC parameter, num. of leapfrog steps at each iteration
#' @param epsilon:  HMC parameter, leapfrog step-size
#' 
#' @references 
#' Adapted to R from original Matlab code by F. Caron:
#' Copyright (C) Francois Caron, University of Oxford
#' see: http://www.stats.ox.ac.uk/~caron/code/bnpgraph/index.html
#' Caron, F., & Fox, E. B. (2017). Sparse graphs using exchangeable random measures. Journal of the Royal Statistical Society. Series B, Statistical Methodology, 79(5), 12
#' @keywords internal
update_w <- function(w, w_rem, sigma, tau, N_counts, L, epsilon) {
  
  # Sample proposal ----------
  
  p <- rnorm(length(w))
  
  p_prop <- p - (epsilon / 2) * grad_loglik_log_w(N_counts, w, sigma, tau, w_rem)
  log_w_prop <- log(w)
  
  for (l in 1:L) {
    
    log_w_prop <- log_w_prop + epsilon * p_prop
    
    if (l < L){
      p_prop <- p_prop - epsilon * grad_loglik_log_w(N_counts, 
                                                     exp(log_w_prop),
                                                     sigma,
                                                     tau,
                                                     w_rem)
    }
  }
  
  w_prop <- exp(log_w_prop)
  
  p_prop <- p_prop - (epsilon / 2) * grad_loglik_log_w(N_counts, 
                                                       w_prop,
                                                       sigma,
                                                       tau,
                                                       w_rem)
  
  # Compute acceptance rate ----------
  
  sum_w_all <- sum(w) + w_rem
  sum_w_prop_all <- sum(w_prop) + w_rem
  log_w <- log(w)
  log_posterior_ratio <- sum( (N_counts - sigma) * (log_w_prop - log_w) ) +
    sum_w_all^2 - sum_w_prop_all^2 +  
    tau * (sum(w) - sum(w_prop))
  
  log_K_p <- -0.5 * sum(p_prop^2 - p^2)
  
  log_acceptance <- log_posterior_ratio + log_K_p
  
  # Sample ----------
  
  if (is.na(log_acceptance)) {
    log_acceptance <- -Inf
  }
  
  if (log(runif(1)) < log_acceptance) {
    w <- w_prop
  } 
  
  rate <- min(1, exp(log_acceptance))
  
  w_updated <- list(w = w,
                    rate = rate,
                    log_posterior_ratio = log_posterior_ratio)
  
  return(w_updated)
}


#' Compute gradient of negative log-likelihood of log w ----------------------------
#' @param N_counts: vector of positive integers, nodes degrees (previous sample, or observed)
#' @param w: vector of positive scalars, nodes sociability (previous sample)
#' @param sigma: real in (-Inf, 1) (previous sample or fixed)
#' @param tau: positive scalar (previous sample or fixed)
#' @param w_rem: positive scalar; total mass of discarded nodes (previous sample)
#' 
#' @references 
#' Adapted to R from original Matlab code by F. Caron:
#' Copyright (C) Francois Caron, University of Oxford
#' see: http://www.stats.ox.ac.uk/~caron/code/bnpgraph/index.html
#' Caron, F., & Fox, E. B. (2017). Sparse graphs using exchangeable random measures. Journal of the Royal Statistical Society. Series B, Statistical Methodology, 79(5), 12
#' @keywords internal
#' 
# Gradient of negative log-likelihood of log w ----------------------------

grad_loglik_log_w <- function(N_counts, w, sigma, tau, w_rem) {
  
  - (N_counts - sigma) + w * (tau + 2 * sum(w) + 2 * w_rem)
  
}
