
#' Computes log-likelihood of data given specified w under TGGP or GGP model
#' @param rows_obs_edges: vector of integers; row indices of 1s in lower triangular binary adjacency matrix of undirected network
#' @param cols_obs_edges: vector of integers; column indices of 1s in lower triangular binary adjacency matrix of undirected network
#' @param rows_nobs_edges: vector of integers; row indices of 0s in lower triangular binary adjacency matrix of undirected network
#' @param cols_nobs_edges: vector of integers; column indices of 0s in lower triangular binary adjacency matrix of undirected network
#' @param w: vector of scalars: nodes sociabilities at which likelihood is evaluated 
#' @param Pi: matrix of probability vectors: nodes community memberships at which likelihood is evaluated (if model is TGGP)
#' @param model: string; either "tggp" or "ggp"
#' @export
compute_log_likelihood_Z <- function(rows_obs_edges, 
                                     cols_obs_edges,
                                     rows_nobs_edges,
                                     cols_nobs_edges,
                                     w,
                                     Pi = NULL,
                                     model = "tggp") {
  
  log_likelihood_Z <- 0
  w_i_w_j_obs <- w[rows_obs_edges] * w[cols_obs_edges]
  w_i_w_j_nobs <- w[rows_nobs_edges] * w[cols_nobs_edges]
  
  if (model == "tggp") {
    
    sum_pi_i_pi_j_obs <- rowSums(Pi[rows_obs_edges, ] * Pi[cols_obs_edges, ])
    sum_pi_i_pi_j_nobs <- rowSums(Pi[rows_nobs_edges, ] * Pi[cols_nobs_edges, ])
    
    # Log-likelihood of edges
    log_likelihood_Z <- log_likelihood_Z +
      sum(
        log(
          1 - 
            exp(- 2 * w_i_w_j_obs) -
            exp(- sum_pi_i_pi_j_obs) +
            exp(- 2 * w_i_w_j_obs - sum_pi_i_pi_j_obs + 
                  2 * w_i_w_j_obs * sum_pi_i_pi_j_obs)
          
        ))

    # Log-likelihood of non-edges
    log_likelihood_Z <- log_likelihood_Z +
      sum(
        log(
          exp(- 2 * w_i_w_j_nobs) +
            exp(- sum_pi_i_pi_j_nobs) -
            exp(- 2 * w_i_w_j_nobs - sum_pi_i_pi_j_nobs + 
                  2 * w_i_w_j_nobs * sum_pi_i_pi_j_nobs)
          
        ))    
    
  } else if (model == "ggp") {
    
    # Log-likelihood of edges
    log_likelihood_Z <- log_likelihood_Z +
      sum(log(1 - exp(- 2 * w_i_w_j_obs)))
    
    # Log-likelihood of non-edges
    log_likelihood_Z <- log_likelihood_Z +
      sum(- 2 * w_i_w_j_nobs)
    
  }

  return(log_likelihood_Z)
  
}


#' Computes log-posterior of log w given current GGP parameters values
#' @param N_counts vector with total number of edges for each node
#' @param w: vector of positive scalars
#' @param w_rem: positive scalar; total mass of discarded nodes (previous sample)
#' @param sigma: real in (-Inf, 1) (previous sample or fixed)
#' @param tau: positive scalar (previous sample or fixed)
#' @export
compute_log_posterior_w <- function(N_counts, w, w_rem, sigma, tau){
  
  log_w <- log(w)
  w_sum <- sum(w)
  
  return(sum( (N_counts - sigma) * log_w - (w_sum + w_rem)^2 - tau * w_sum ))
  
}

#' Computes log-posterior of log w given current GGP parameters values
#' @param M_counts matrix with number of edges in every community (TGGP) or vector with total number of edges (GGP), for each node
#' @param Pi matrix of nodes community memberships
#' @param beta_vec vector of overall community frequencies
#' @param alpha_0 hyperparameter of prior distribution of Pi
#' @export
compute_log_posterior_Pi <- function(M_counts, Pi, beta_vec, alpha_0) {
  
  beta_updated <-  M_counts + alpha_0 * beta_vec[col(M_counts)]
  log_posterior_Pi <- compute_log_posterior_Pi_cpp(Pi, beta_updated)
  
  return(log_posterior_Pi)
}

#' Computes Shannon-entropy for each row of a NxK matrix of N probability vectors over K classes
#' @param Pi  matrix of nodes community memberships
#' @export
compute_entropy_Pi <- function(Pi) {
  # Shannon entropy is defined as  -∑ p*log(p)
  # If p is very close to 0, then log(p) could be unstable
  # For numerical stability, rewrite
  # -∑ p*log(p) = log{ exp[ -∑ p*log(p) ] }
  #             = log { prod [ exp (- log(p^p)) ]} 
  
  compute_entropy_Pi_cpp(Pi)
  
}


