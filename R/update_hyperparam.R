#' Update total mass of removed nodes and hyperparameters given the rest,
#' using MH as described in Caron and Fox, Appendix F.2
#'
#' @param w: vector of positive scalars, nodes sociability (previous sample)
#' @param w_rem: positive scalar; total mass of discarded nodes (previous sample)
#' @param alpha: positive scalar; truncation of GGP (previous sample or fixed)
#' @param sigma: real in (-Inf, 1) (previous sample or fixed)
#' @param tau: positive scalar (previous sample or fixed)
#' @param ggp_param: list including information on alpha, sigma and tau
#' 
#' @references 
#' Adapted to R from original Matlab code by F. Caron:
#' Copyright (C) Francois Caron, University of Oxford
#' see: http://www.stats.ox.ac.uk/~caron/code/bnpgraph/index.html
#' Caron, F., & Fox, E. B. (2017). Sparse graphs using exchangeable random measures. Journal of the Royal Statistical Society. Series B, Statistical Methodology, 79(5), 12
#' @keywords internal
update_ggp_hyperparam <- function(w, w_rem, alpha, sigma, tau,
                                  ggp_param, rw_alpha = FALSE) {
  
  # Settings
  rw_std <- c(0.2, 0.2)
  log_w <- log(w)
  log_alpha <- log(alpha)
  num_nodes <- length(w)
  hyper_alpha <- c(0, 0)
  hyper_sigma <- c(0, 0)
  hyper_tau <- c(0, 0)
  
  # Sample proposals ----------
  
  if (ggp_param$estimate_sigma) {
    sigma_prop <- 1 - exp( log(1 - sigma) + rw_std[1] * rnorm(1))
  } else{
    sigma_prop <- sigma
  }
  
  if (ggp_param$estimate_tau) {
    tau_prop <- exp( log(tau) + rw_std[2] * rnorm(1))
  } else {
    tau_prop <- tau
  }
  
  sum_w <- sum(w)
  sum_all_w <- sum_w + w_rem
    
    if (sigma_prop > -1) { # Taken from Caron's code, not mentioned in the paper
      # This appears to be done for numerical stability,
      # as alpha can take very large values for sigma < -1
      # so it is preferable to adopt a different strategy
      # for the two cases
      
      if (ggp_param$estimate_alpha) {
        if (rw_alpha == FALSE) { # Gamma proposal, as in Appendix F.2
          alpha_prop <- rgamma(1, 
                               num_nodes, 
                               psi_GGP(2 * sum_all_w, 1, sigma_prop, tau_prop))
        } else { # Random walk proposal
          alpha_prop <- alpha * exp(0.02 * rnorm(1))
        }
      } else {
        alpha_prop <- alpha
      }    
      log_alpha_prop <- log(alpha_prop)
      w_prop_rem <- rGGPsum(alpha_prop, sigma_prop, tau_prop + 2*sum_w + 2*w_rem)
      
    } else { # more stable numerically as alpha can take very large values 
      # in that case, we sample alpha2=alpha*tau^sigma
      if (ggp_param$estimate_alpha) {
        if(rw_alpha == FALSE) { # gamma proposal
          alpha2_prop <- rgamma(1, 
                                num_nodes,
                                ( psi_GGP((2 * sum_w + 2 * w_rem) / tau_prop,
                                          1, 
                                          sigma_prop, 1) ))
          log_alpha_prop <- log(alpha2_prop) - sigma_prop*log(tau_prop)
        } else { # Random walk
          log_alpha_prop <- log_alpha + .02 * rnorm(1)
        }
        alpha_prop <- exp(log_alpha_prop)
        rate_K <- exp( log_alpha_prop - log(-sigma_prop) + 
                         sigma_prop*log(tau_prop + 2*sum_w + 2*w_rem ) )
        num_clust <- rpois(1, rate_K)      
        w_prop_rem = rgamma(1, 
                            -sigma_prop* num_clust, 
                            1 / (tau_prop + 2*sum_w + 2*w_rem))
      } else {
        alpha_prop <- alpha
        log_alpha_prop <- log_alpha
        w_prop_rem <- rGGPsum(alpha_prop, sigma_prop, tau_prop + 2*sum_w + 2*w_rem)
        
      }
    }
    
    # Compute acceptance rate ----------
    
    sum_w_prop <- sum(w)
    sum_all_w_prop <- sum_w_prop + w_prop_rem
    
    log_r_1 <- - sum_all_w_prop^2 + sum_all_w^2 +
      - (tau_prop - tau + 2 * w_rem - 2 * w_prop_rem) * sum_w +
      (sigma - sigma_prop) * sum(log_w) 
    
    log_r_2 <- num_nodes * (lgamma(1 - sigma) - lgamma(1 - sigma_prop))
    
    log_acceptance <- log_r_1 + log_r_2
    
    # Adjust for priors
    
    if (ggp_param$estimate_alpha) {
      
      if (rw_alpha == FALSE) {
        log_acceptance <- log_acceptance +
          + num_nodes * 
          (log(psi_GGP((2 * sum_w_prop + 2*w_prop_rem) / tau, 1, sigma, 1) ) 
           + sigma * log(tau)
           - log(psi_GGP((2*sum_w + 2*w_rem) / tau_prop, 1, sigma_prop, 1) ) - 
             sigma_prop * log(tau_prop))
        
      } else {
        log_acceptance <- log_acceptance +
          - exp(log_alpha_prop + sigma_prop * log(tau_prop)) *  
          psi_GGP((2 * sum_w + 2 * w_rem) / tau_prop, 1, sigma_prop, 1) +
          exp(log_alpha + sigma*log(tau)) * 
          psi_GGP((2*sum_w_prop + 2*w_prop_rem) / tau, 1, sigma, 1) +
          num_nodes * (log_alpha_prop - log_alpha)
      
      }
      if (hyper_alpha[1] > 0) {
        log_acceptance <- log_acceptance + 
          hyper_alpha[1] * ( log_alpha_prop - log_alpha)
      } 
      
      if (hyper_alpha[2] > 0) {
        log_acceptance <- log_acceptance - hyper_alpha(2) * (alpha_prop - alpha)
      }
      
    } else {
      log_acceptance <- log_acceptance +
        - psi_GGP(2*sum_w + 2*w_rem, alpha_prop, sigma_prop, tau_prop) +
        + psi_GGP(2*sum_w_prop + 2*w_prop_rem, alpha, sigma, tau)
      
    }
    
    if (ggp_param$estimate_tau) {
      log_acceptance <- log_acceptance +
        hyper_tau[1] * ( log(tau_prop) - log(tau)) - 
        hyper_tau[2] * (tau_prop - tau)
      
    }
    
    if (ggp_param$estimate_sigma) {
      
      log_acceptance <- log_acceptance +
        hyper_sigma[1] * ( log(1 - sigma_prop) - log(1 - sigma)) -
        hyper_sigma[2] * (1 - sigma_prop - 1 + sigma)
      
    }
    
    if (is.na(log_acceptance)) {
      stop("Something is not right, log acceptance of MH step for hyperparams is NA")
    }
    
    # Accept/Reject step ------------------------------------------------------
    
    if (log(runif(1)) < log_acceptance) {
      
      w_rem <- w_prop_rem
      alpha <- alpha_prop
      sigma <- sigma_prop
      tau <- tau_prop
      
    }
    
    rate_MH <- min(1, exp(log_acceptance))
    
    ggp_param_updated <- list(w_rem = w_rem,
                              alpha = alpha,
                              sigma = sigma,
                              tau = tau,
                              rate_MH = rate_MH)
    
    return(ggp_param_updated)

}


#' Laplace exponent of GGP: Psi function
#' 
#' psi(tt) = -log ( E[exp(-tt * sum_i w_i)] )
#'        = alpha/sigma * ( (tt+tau)^sigma - tau^sigma))
#' where
#' (w_i)_{i=1,2,..} are the points of a Poisson process on R_+ of mean measure 
#' rho(dw) = alpha/Gamma(1-sigma) * w^{-1-sigma} * exp(-tau*w)dw
#'
#' @param tt: value at which to evaluate the Laplace exponent
#' @param alpha: positive scalar; truncation of GGP (previous sample or fixed)
#' @param sigma: real in (-Inf, 1) (previous sample or fixed)
#' @param tau: positive scalar (previous sample or fixed)
#' 
#' @references 
#' Adapted to R from original Matlab code by F. Caron:
#' Copyright (C) Francois Caron, University of Oxford
#' see: http://www.stats.ox.ac.uk/~caron/code/bnpgraph/index.html
#' Caron, F., & Fox, E. B. (2017). Sparse graphs using exchangeable random measures. Journal of the Royal Statistical Society. Series B, Statistical Methodology, 79(5), 12
#' @keywords internal
psi_GGP <- function(tt, alpha, sigma, tau) {

  if (sigma == 0) {  # gamma process
    alpha * log( 1 + (tt / tau) )
    
  } else {
    (alpha / sigma) * ((tt + tau)^sigma - tau^sigma)
    
  }
}

#' Sample from distribution of GGP total mass
#' 
#' It generates a realization of the random variable S with Laplace transform
#' E[e^-(tt*S)] = exp(-alpha/sigma * [(tt+tau)^sigma - tau^sigma])
#'
#' @param tt: value at which to evaluate the Laplace exponent
#' @param alpha: positive scalar; truncation of GGP (previous sample or fixed)
#' @param sigma: real in (-Inf, 1) (previous sample or fixed)
#' @param tau: positive scalar (previous sample or fixed)
#' 
#' @references 
#' Adapted to R from original Matlab code by F. Caron:
#' Copyright (C) Francois Caron, University of Oxford
#' see: http://www.stats.ox.ac.uk/~caron/code/bnpgraph/index.html
#' Caron, F., & Fox, E. B. (2017). Sparse graphs using exchangeable random measures. Journal of the Royal Statistical Society. Series B, Statistical Methodology, 79(5), 12
#' @keywords internal
rGGPsum <- function(alpha, sigma, tau) {
  
  if (sigma < -10^-8) {
    
    # Compound Poisson case (dense graph)
    # S is distributed from a Poisson mixture of gamma variables
    D <- rpois(1, - alpha / sigma / tau^(-sigma))
    S <- rgamma(1, -sigma * D, tau)
    
  } else if (sigma < 10^-8) {
    
    # Gamma process case
    # S is Gamma distributed
    S <- rgamma(alpha, tau)
    
  } else if (sigma == 0.5 & tau == 0) { 
    
    # Inverse Gaussian process case
    # S is inverse Gaussian distributed
    lambda <- 2 * alpha^2
    mu <- alpha / sqrt(tau)
    S <- statmod::rinvgauss(1, mu, lambda)
    
  } else {
    
    # General case
    # S is distributed from an exponentially tilted stable distribution
    S <- r_et_stable(alpha/sigma, sigma, tau)
  }  
  
  return(S)

}


#' Sample from from the exponentially tilted stable distribution
#' distributed fron an exponentially tilted stable distribution 
#' with Laplace transform (in z) exp(-V0 * ((z + tau)^alpha - tau^alpha))
#' Uses the algorithm proposed in (Devroye, 2009) with corrections pointed 
#' out in (Hofert, 2011)
#' @param V0: positive scalar
#' @param alpha: positive scalar; truncation of GGP (previous sample or fixed)
#' @param tau: positive scalar (previous sample or fixed)
#' @param n: integer
#' 
#' @value
#' S: vector of length n
#' 
#' @references 
#' Adapted to R from original Matlab code by F. Caron:
#' Copyright (C) Francois Caron, University of Oxford
#' see: http://www.stats.ox.ac.uk/~caron/code/bnpgraph/index.html
#' Caron, F., & Fox, E. B. (2017). Sparse graphs using exchangeable random measures. Journal of the Royal Statistical Society. Series B, Statistical Methodology, 79(5), 12
#' Luc Devroye. Random variate generation for exponentially and polynomially tilted stable distributions. ACM Transactions on Modeling and Computer Simulation, vol. 19(4), 2009.
#' Marius Hofert. Sampling exponentially tilted stable distributions. ACM Transactions on Modeling and Computer Simulation, vol. 22(1), 2011
#'
#' @keywords internal
r_et_stable <- function(V0, alpha, tau, n = 1) {
  
  lambda_alpha <- tau^alpha * V0
  
  # Sample from an exponentially tilted distribution of parameters
  # sigma, lambda, as in (Devroye, 2009)
  
  gam <- lambda_alpha * alpha * (1 - alpha)
  
  xi <- (1 / pi) * ((2 + sqrt(pi / 2)) * sqrt(2 * gam) + 1) # Correction in Hofert
  psi = (1 / pi) * exp(-gam * pi^2/8) * (2 + sqrt(pi/2)) * sqrt(gam * pi)
  w1 = xi * sqrt(pi / 2 / gam)
  w2 = 2 * psi * sqrt(pi)
  w3 = xi * pi
  b = (1 - alpha) / alpha
  
  samples <- rep(0, n)
  
  for (i in c(1:n)) {
    
    while(TRUE) { # Generate U with density g*/G*
      
      while(TRUE) { # Generate U with density proportional to g**
        
        U <- gen_U(w1, w2, w3, gam)
        W <- runif(1)
        zeta <- sqrt(ratio_B(U, alpha))
        z <-  1 / (1 - (1 + alpha * zeta / sqrt(gam))^(-1 / alpha))
        
        rho <- pi * exp(-lambda_alpha * (1 - zeta^(-2))) *
          ( xi * exp(-gam * U^2 / 2) * (U >= 0) * (gam >= 1) + 
              psi / sqrt(pi - U) * (U > 0) * (U < pi) + 
              xi *(U >= 0) * (U <= pi) * (gam < 1)) /
          ( (1 + sqrt(pi / 2)) * sqrt(gam) / zeta + z)
        
        if ( (U < pi) & (W * rho <= 1) ) {
          break
        }
      }
      
      # Generate X with density proportional to g(x, U)
      
      a <- zolotarev(U, alpha)
      m <- (b / a)^alpha * lambda_alpha
      delta <- sqrt(m * alpha / a)
      a1 <- delta * sqrt(pi / 2)
      a2 <- a1 + delta # correction in Hofert
      a3 <- z / a
      s <- a1 + delta + a3 # correction in Hofert
      V_p <- runif(1)
      N_p <- rnorm(1)
      E_p <- -log(runif(1))
      
      if (V_p < (a1 / s)) {
        X <- m - delta * abs(N_p)
      } else if (V_p < (a2 / s)) {
        X <- delta * runif(1) + m
      } else {
        X <- m + delta + a3 * E_p
      }
      
      E <- -log(runif(1))
      cond <- ( a * (X - m)) + 
        exp(1 / alpha * log(lambda_alpha) - b * log(m)) *
        ((m / X)^b - 1) - ( (N_p^2 / 2) * (X < m) ) -
        ( E_p * (X > m + delta) )
      
      if ( (X>= 0) & (cond <= E)) {
        break
      }
    }
    samples[i] <-  exp( 1 / alpha* log(V0) -b * log(X)) # more stable than V0^(1/alpha) * X^(-b);
  }
  
  return(samples)
}


#' Utility function for r_et_stable (generates U)
#' @param w1: scalar
#' @param w2: scalar
#' @param w3: scalar
#' @param gam: scalar
#' @value U: scalar
#' @keywords internal
gen_U <- function(w1, w2, w3, gam) {
  
  V <- runif(1)
  W_p <- runif(1)
  
  if (gam > 1) {
    
    if (V < w1 / (w1 + w2)) {
      U <- abs(rnorm(1)) / sqrt(gam)
    } else {
      U <- pi * (1 - W_p^2)
    }
    
  } else {
    
    if (V < w3 / (w3 + w2)) {
      U <- pi * W_p
    } else {
      U <- pi * (1 - W_p^2)
    }
  }
  
}

#' Utility function for r_et_stable (computes sin(x)/x)
#' @param x: scalar
#' @keywords internal
sinc <- function(x) {
  sin(x) / x
}

#' Utility function for r_et_stable
#' @param x: scalar
#' @param sigma: real in (-Inf, 1) (previous sample or fixed)
#' @keywords internal
ratio_B <- function(x, sigma) {
  sinc(x) / 
    (sinc(sigma * x))^sigma / 
    (sinc((1 - sigma) * x))^(1 - sigma)
}

#' Utility function for r_et_stable (evaluates Zolotarev function)
#' @param u: scalar
#' @param sigma: real in (-Inf, 1) (previous sample or fixed)
#' @keywords internal
zolotarev <- function(u, sigma) {
  # Zolotarev function, cf (Devroye, 2009)
  ( (sin(sigma * u))^sigma * 
      (sin((1 - sigma) * u))^(1 - sigma) / sin(u) )^(1 / (1-sigma) )
}


#' Update total mass of removed nodes and hyperparameters given the rest,
#' using MH as described in Caron and Fox, Appendix F.2
#'
#' @param alpha: positive scalar; truncation of GGP (previous sample or fixed)
#' @param sigma: real in (-Inf, 1) (previous sample or fixed)
#' @param tau: positive scalar (previous sample or fixed)
#' @param beta_vec: vector of overall frequencies of each K_max communities (previous sample or fixed)
#' @param num_N_tot_samples: integer: number of samples to base update on
#' 
#' @keywords internal
update_N_tot_to_N_ratio <- function(alpha,
                                    sigma,
                                    tau,
                                    beta_vec,
                                    num_N_samples) {
  N_tot_to_N_ratio <- 0
  K_max <- length(beta_vec)
  
  for (s in 1:num_N_samples){
    tggp_s <- rand_tggp_network(alpha = alpha, 
                                sigma = round(sigma, 1), # This is because of the weird behavior or rand_ggp when 0.1 < sigma < 0.15
                                tau = tau, 
                                K = K_max,
                                beta_vec = beta_vec,
                                only_simulate_N = TRUE)
    N <- max(tggp_s$N_unthinned, 1)
    N_tot <- N + tggp_s$N_thinned
    N_tot_to_N_ratio <- N_tot_to_N_ratio + (N_tot / N)
  }
  
  # Take mean
  N_tot_to_N_ratio <- N_tot_to_N_ratio / num_N_samples
  
  return(N_tot_to_N_ratio)
}

#' Update total mass of removed nodes and hyperparameters given the rest,
#' using MH as described in Caron and Fox, Appendix F.2
#'
#' @param N_tot: int; total number of nodes (unthinned and thinned)
#' @param N_tot_old: int; total number of nodes (unthinned and thinned) before latest update
#' @param N: int; total number of observed (i.e. unthinned) nodes
#' @param w_nobs: vector of scalars: nodes sociabilities at which likelihood is evaluated (unobserved nodes)
#' @param Pi_nobs: matrix of probability vectors: nodes community memberships at which likelihood is evaluated (unobserved nodes)
#' @param beta_vec: vector of community frequencies
#' 
#' @keywords internal
trim_enlarge_w_Pi_nobs <- function(N_tot, N_tot_old, N, w_nobs, Pi_nobs, beta_vec) {
  
  N_nobs <- N_tot - N
  w_nobs_old <- w_nobs
  Pi_nobs_old <- Pi_nobs
  
  if (N_tot < N_tot_old) {
    # Retain only the largest N_nobs w for thinned nodes
    w_nobs <- sort(w_nobs_old)[1:N_nobs]
    # Retain Pi of corresponding nodes
    Pi_nobs <- Pi_nobs_old[order(w_nobs_old),][(1:N_nobs),]
    
  } else {
    # Number of extra thinned nodes
    N_nobs_extra <- N_tot - N_tot_old
    # Resample among extra w_nobs among existing w_nobs
    w_nobs_extra <- sample(w_nobs_old, N_nobs_extra, replace = TRUE)
    w_nobs <- c(w_nobs_old, w_nobs_extra)
    # Resample extra community memberships from prior
    Pi_nobs_extra <- MCMCpack::rdirichlet(N_nobs_extra, beta_vec)
    Pi_nobs <- rbind(Pi_nobs_old, Pi_nobs_extra)
  }
  
  return(list(w_nobs = w_nobs, Pi_nobs = Pi_nobs))
}

