#' (Approximate) sample from the Generalized Gamma Process 
#'
#' Samples points of the Generalized Gamma Process (GGP) with Lévy measure
#' \deqn{[alpha / Gamma(1-sigma)] * w^(-1-sigma)  * exp(-tau * w)}
#' When sigma < 0, exact sampling is possible.
#' For sigma >= 0, approximate samples are obtained using the adaptive 
#' thinning strategy described in Favaro and Teh (2013).  
#'
#' @param alpha positive real number. Affects the number of nodes and edges (higher alpha leads to larger networks).
#' @param sigma real number. in (-Inf, 1). sigma < 0 corresponds to a finite-activity GGP (dense case).
#' @param tau positive real number. Affects the tails of the degree distribution (higher tau leads to faster decay).
#' 
#' @return A real-valued vector with an approximate sample from the GGP. 
#' 
#' @examples
#' w <- simulate_ggp(20, 0.5, 1)
#'
#' @source 
#' R Code adapted from original Matlab code by F. Caron:
#' Copyright (C) Francois Caron, University of Oxford.
#' 
#' See: http://www.stats.ox.ac.uk/~caron/code/bnpgraph/index.html
#' 
#' @references 
#' Caron, F., & Fox, E. B. (2017). Sparse graphs using exchangeable random measures. Journal of the Royal Statistical Society. Series B, Statistical Methodology, 79(5), 12
#' 
#' S. Favaro and Y.W. Teh. MCMC for normalized random measure mixture models. Statistical Science, vol.28(3), pp.335-359, 2013.
#' 
#' @export
rand_ggp <- function(alpha, sigma, tau) {
  
  if ((alpha <= 0) | (sigma >= 1) | (tau <= 0)) {
    stop("Inputs must be as follows: alpha > 0, sigma < 1 and tau > 0")
  }
  

  if (sigma < 0) {
    
    # Finite activity GGP ----------------------------------------------------
    
    rate <- exp( log(alpha) - log(-sigma) + sigma * log(tau) )
    K <- rpois(1, rate)
    w <- rgamma(K, -sigma, tau)
    w <- w[(w > 0)] # Taken from original code, not sure why it's needed
    
  } else {
    
    # Infinite activity GGP --------------------------------------------------
    
    # Set truncation and expected number of jumps
    if (sigma > 0.1) {
      
      Njumps <- 20000  
      trunc <- exp( 1 / sigma * 
                  (log(alpha) - log(sigma) - lgamma(1 - sigma) - log(Njumps))) 
      
    } else {
      
      trunc <- 1e-10
      
      if (sigma > 0) {
        
        Njumps <- floor(alpha / sigma / gamma(1 - sigma) * trunc^(-sigma))
        
      } else {
        
        Njumps <- floor(-alpha * log(trunc))
        
      }
      
    } 
    
    # Use adaptive thinning strategy from Favaro and Teh (2013)
    
    w <- array(0, c(1, ceiling(Njumps + 3 * sqrt(Njumps))))
    k <- 1
    trunc_old <- trunc
    iter <- 0
    keep_adapting <- TRUE

    while (keep_adapting) {
      
      # Sample exponential random variable of unit rate
      e <- -log(runif(1))
      
      if (e > compute_W(trunc_old, Inf, alpha, sigma, tau)) {
        
        w <- w[1 : (k-1)]
        return(w)
        
      } else{
        
        trunc_new <- compute_inv_W(trunc_old, e, alpha, sigma, tau)
        
      }
      
      if ( tau == 0 || (log(runif(1)) < ((-1 - sigma) * log(trunc_new / trunc_old))) ) {
        
        # if (tau > 0), adaptive thinning - otherwise accept always
        w[k] <- trunc_new
        k <- k + 1
        
      }
      
      trunc_old <- trunc_new
      iter <- iter + 1
      
      if (iter > 10^8){ 
        
        # If too many computations, should lower the threshold trunc to
        # trunc_new = round(trunc/10, 10). This has not been implemented yet
        stop('Too many computations were required. Try rerun or change parameters.')

      }
    }
  }
  
}


# SUBFUNCTIONS ------------------------------------------------------------


compute_W <- function(trunc_old, x, alpha, sigma, tau) {
  
  if (tau > 0){
    
    logout <- log(alpha) + log(1 - exp(-tau * (x - trunc_old))) +
      (-1 - sigma) * log(trunc_old) + (-trunc_old * tau) - log(tau) - lgamma(1 - sigma)
    
  } else {
    
    logout <- log(alpha) - lgamma(1 - sigma) - log(sigma) + 
      log(trunc_old^(-sigma) - x^(-sigma))
    
  }
  
  return(exp(logout))
  
}


compute_inv_W <- function(trunc_old, x, alpha, sigma, tau) {
  
  if (tau > 0){
    
    out <- trunc_old - (1 / tau) * 
      log( 1 - gamma(1 - sigma) * x * tau / 
             (alpha * trunc_old^(-1 - sigma) * exp(-trunc_old * tau))
      )
    
  } else {
    
    logout <- -1 / sigma *
      log( trunc_old^(-sigma)  - sigma * gamma(1 - sigma) / alpha * x)
    
    out <- exp(logout)
    
  }
  
  return(out)
  
}

