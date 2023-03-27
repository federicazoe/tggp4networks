#' Fit Generalized Gamma Process to a simple undirected network via Gibbs Sampling 
#' following the method described in Caron and Fox (2017)
#'
#' ADD DESCRIPTION HERE
#' 
#' @param Z TsparseMatrix; lower-triangular binary adjacency matrix of undirected network.
#' @param mcmc_settings list providing desired mcmc settings (e.g. number of iterations). See Details.
#' @param ggp_param list (optional) can include fields alpha, sigma and tau
#' @param general_settings optional list specifying general settings, such as whether to print hyperparameter values over iterations. See Details.
#' 
#' @details 
#' 
#' @value
#' \subsection{A list with the following components: } {
#'    \describe{
#'      \item{mcmc_samples} {MCMC samples of model parameters (nodes' sociabilities) and hyperparameters (GGP parameters)).}
#'      \item{mcmc_acceptance} {Hamiltonian Monte Carlo and Metropolis Hastings acceptance rates.}
#'      \item{mcmc_settings} {List with selected MCMC settings (e.g. number of iterations).}
#'      \item{general_settings} {List with selected logistics settings (e.g. path for saving).}
#'      \item{initial_values} {List with lists of initial hyperparameters' and parameters' values.}
#'    }
#' }
#' 
#' @examples
#' \dontrun{
#' # Fit TGGP model on a simulated network with 3 communities of similar size
#' tggp_network <- sample_tggp_network(alpha = 50, sigma = 0.5, tau = 1, K = 3)
#' Z <- tggp_network$Z
#' mcmc_settings <- list(n_iter = 10)
#' general_settings <- list(path_for_saving = "experiments/results"),
#' ggp_param <- list(K = 5)
#' Z_fit <- fit_tggp_mcmc(Z, mcmc_settings, general_settings, ggp_param)
#' }
#' 
#' @references 
#' F. Ricci, M. Guindani, and E. Sudderth (2022). Thinned Random Measures for Sparse Graphs with Overlapping Communities. Neural Information Processing Systems 
#' 
#' Caron, F., & Fox, E. B. (2017). Sparse graphs using exchangeable random measures. Journal of the Royal Statistical Society. Series B, Statistical Methodology, 79(5), 12
#' 
#' @export
fit_ggp_mcmc <- function(Z,
                         mcmc_settings, 
                         ggp_param = list(), 
                         general_settings = list()) {
  
  # Extract number of nodes, indices of nodes with(out) an edge ------------
  N <- Z@Dim[1]
  rows_obs_edges <- Z@i + 1
  cols_obs_edges <- Z@j + 1
  all_nodes_pairs <- as.matrix(expand.grid(1:N, 1:N))
  all_nodes_pairs <- all_nodes_pairs[all_nodes_pairs[, 1] > all_nodes_pairs[, 2],]
  nodes_pairs_edges_collapsed <- paste(rows_obs_edges, cols_obs_edges)
  all_nodes_pairs_collapsed <- paste(all_nodes_pairs[, 1], all_nodes_pairs[, 2])
  non_edges <- !(all_nodes_pairs_collapsed %in% nodes_pairs_edges_collapsed)
  rows_nobs_edges <- all_nodes_pairs[non_edges, 1]
  cols_nobs_edges <- all_nodes_pairs[non_edges, 2]

  # Initialize hyperparameters alpha, log_alpha, sigma, tau ----------------
  ggp_param <- initialize_model_param(ggp_param, N, model = "ggp")
  alpha <- ggp_param$alpha
  log_alpha <- ggp_param$log_alpha
  sigma <- ggp_param$sigma
  tau <- ggp_param$tau
  w_rem <- ggp_param$w_rem

  # Initialize w -----------------------------------------------------------
  w <- rgamma(N, 1, 1)
  log_w <- log(w)

  # General and MCMC settings ----------------------------------------------
  mcmc_settings <- set_mcmc_settings(mcmc_settings, N)
  general_settings <- set_general_settings(general_settings, mcmc_settings$n_iter)

  # Define HMC parameters ---------------------------------------------------
  epsilon <- mcmc_settings$epsilon
  L <- mcmc_settings$L
  
  # Store initialized hyperparameters and parameters -----------------------
  initial_values <- list(
    hyperparams = list(
      alpha = alpha,
      sigma = sigma,
      tau = tau,
      w_rem = w_rem
    ),
    params = list(
      w = w
    )
  )

  # Prepare storage --------------------------------------------------------
  # Samples
  mcmc_samples <- list(
    w = matrix(NA, nrow = mcmc_settings$n_samples, ncol = N), 
    w_rem = rep(NA, mcmc_settings$n_samples), 
    alpha = rep(NA, mcmc_settings$n_samples),
    sigma = rep(NA, mcmc_settings$n_samples), 
    tau = rep(NA, mcmc_settings$n_samples),
    log_likelihood_Z = rep(NA, mcmc_settings$n_iter),
    log_posterior_w = rep(NA, mcmc_settings$n_iter)
  )  
  # Acceptance
  mcmc_acceptance <- list(
    rate_HMC = rep(NA, mcmc_settings$n_iter), 
    rate_MH_hyper = rep(NA, mcmc_settings$n_iter),
    epsilon = rep(NA, mcmc_settings$n_iter)
  )
  # Diagnostics
  mcmc_diagnostics <- list(
    log_likelihood_Z = rep(NA, mcmc_settings$n_iter),
    log_posterior_w = rep(NA, mcmc_settings$n_iter)
  )
  
  # Run MCMC ---------------------------------------------------------------
  for (i in 1:mcmc_settings$n_iter) {
    
    # Communicate the state of the MCMC, if wanted
    if (general_settings$verbose == TRUE  &&  
        (i %% general_settings$speak_every == 0)) {
      print("----------------------------------------")
      print(paste("Iter: ", i))
      print(paste("alpha: ", round(alpha, 2)))
      print(paste("sigma: ", round(sigma, 2)))
      print(paste("tau: ", round(tau, 2)))
      if (i > 1) {
        print(paste("rate_HMC: ", round(mcmc_acceptance$rate_HMC[i-1], 2)))
        print(paste("rate_MC_hyper: ", round(mcmc_acceptance$rate_MH_hyper[i-1], 2))) 
      }
      print("----------------------------------------")
    }
    
    # Update number of edges in which each node participates (N_counts) ------
    N_counts <- update_N_counts_ggp(rows_obs_edges, 
                                    cols_obs_edges,
                                    w)
    
    # Update nodes sociabilities w using HMC
    w_updated <- update_w(w, w_rem, sigma, tau, N_counts, L, epsilon)
    w <- w_updated$w
    mcmc_acceptance$rate_HMC[i] <- w_updated$rate
    mcmc_acceptance$log_posterior_ratio_w[i] <- w_updated$log_posterior_ratio
    mcmc_acceptance$epsilon[i] <- epsilon
    
    # Adapt epsilon
    if (i < mcmc_settings$n_adapt) {
      if (is.null(mcmc_settings$eps_avg_window) || 
          (i < mcmc_settings$eps_avg_window) ) {
        
        epsilon <- exp(log(epsilon) + 
                         mcmc_settings$eps_adjust_rate * 
                         (mean(mcmc_acceptance$rate_HMC[1:i]) - 0.6))
        
      } else {
        start_from <- (i - mcmc_settings$eps_avg_window + 1)
        epsilon <- exp(log(epsilon) + 
                         mcmc_settings$eps_adjust_rate * 
                         (mean(mcmc_acceptance$rate_HMC[start_from:i]) - 0.6))
        
      }
    }
    
    if (ggp_param$estimate_alpha | 
        ggp_param$estimate_sigma | 
        ggp_param$estimate_tau) {
      
      # Update total sociability of removed nodes and hyperparams using MH
      # (alternate between gamma and random walk for alpha)
      if (i %% 2 == 0) {
        rw_alpha <- FALSE
      } else {
        rw_alpha <- TRUE
      }    
      
      ggp_param_updated <- update_ggp_hyperparam(w, w_rem, alpha, sigma, tau,
                                                 ggp_param, rw_alpha = rw_alpha)
      w_rem <- ggp_param_updated$w_rem
      alpha <- ggp_param_updated$alpha
      sigma <- ggp_param_updated$sigma
      tau <- ggp_param_updated$tau
      mcmc_acceptance$rate_MH_hyper[i] <- ggp_param_updated$rate_MH
      
    }
    
    # Store log-likelihood of observed data and log-posterior of Pi and w
    # at current iteration
    mcmc_diagnostics$log_likelihood_Z[i] <- compute_log_likelihood_Z(
      rows_obs_edges, 
      cols_obs_edges,
      rows_nobs_edges,
      cols_nobs_edges,
      w,
      model = "ggp"
    )
    mcmc_diagnostics$log_posterior_w[i] <- compute_log_posterior_w(
      N_counts,
      w,
      w_rem,
      sigma,
      tau
    )
    
    # Store samples
    if (i %in% mcmc_settings$iters_to_save) {
      
      ii <- which(mcmc_settings$iters_to_save == i)
      mcmc_samples$w[ii, ] <- w
      mcmc_samples$w_rem[ii] <- w_rem
      mcmc_samples$alpha[ii] <- alpha
      mcmc_samples$sigma[ii] <- sigma
      mcmc_samples$tau[ii] <- tau
      
    }
    
  }
  
  
  return(list(mcmc_samples = mcmc_samples,
              mcmc_acceptance = mcmc_acceptance,
              mcmc_diagnostics = mcmc_diagnostics,
              initial_values = initial_values))
  
}
