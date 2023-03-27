#' Fit Thinned Generalized Gamma Process to a simple undirected network via Gibbs Sampling
#'
#' ADD DESCRIPTION HERE
#' 
#' @param Z TsparseMatrix; lower-triangular binary adjacency matrix of undirected network.
#' @param mcmc_settings list specifying desired number of iterations and (optionally) other mcmc settings. See Details.
#' @param tggp_param list specifying desired maximum number of communities and (optionally) other hyperparameters (only if user wishes to fix them or provide initial values). See Details.
#' @param general_settings optional list specifying general settings, such as whether to print hyperparameter values over iterations. See Details.
#' 
#' @details 
#' 
#' @return
#' \subsection{A list with the following components: } {
#'    \describe{
#'      \item{mcmc_samples} {MCMC samples of model parameters (nodes' sociabilities and community memberships) and hyperparameters (relative size of communities, GGP parameters)).}
#'      \item{mcmc_acceptance} {Hamiltonian Monte Carlo and Metropolis Hastings acceptance rates.}
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
#' model_hyper <- list(K = 5)
#' Z_fit <- fit_tggp_mcmc(Z, mcmc_settings, general_settings, model_hyper)
#' }
#' 
#' @details 
#' ADD DETAILS ON ARGUMENTS' NAMED FIELDS
#' 
#' @references 
#' Ricci, F. Z., Guindani, M., & Sudderth, E. (2022). Thinned Random Measures for Sparse Graphs with Overlapping Communities. Advances in neural information processing systems, 35
#' 
#' Caron, F., & Fox, E. B. (2017). Sparse graphs using exchangeable random measures. Journal of the Royal Statistical Society. Series B, Statistical Methodology, 79(5), 12
#' 
#' Teh, Y., Jordan, M., Beal, M., & Blei, D. (2004). Sharing clusters among related groups: Hierarchical Dirichlet processes. Advances in neural information processing systems, 17.
#' 
#' @export
fit_tggp_mcmc <- function(Z,
                          mcmc_settings, 
                          tggp_param, 
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

  # General and MCMC settings ----------------------------------------------
  mcmc_settings <- set_mcmc_settings(mcmc_settings, N)
  general_settings <- set_general_settings(general_settings, mcmc_settings$n_iter)
  
  # Initialize all hyperparameters and parameters --------------------------
  initialization_lists <- initialize_tggp(Z,
                                          mcmc_settings$n_iter_initial,
                                          tggp_param,
                                          general_settings$verbose)
  initial_values <- initialization_lists$initial_values
  tggp_param <- initialization_lists$tggp_param
  
  K <- tggp_param$K_max
  alpha <- initial_values$alpha
  sigma <- initial_values$sigma
  tau <- initial_values$tau
  w_rem <- initial_values$w_rem
  beta_vec <- initial_values$beta
  w <- initial_values$w
  log_w <- initial_values$log_w
  Pi <- initial_values$Pi
  N_tot <- initial_values$N_tot
  w_nobs <- initial_values$w_nobs
  Pi_nobs <- initial_values$Pi_nobs
  gamma_0 <- tggp_param$gamma_0
  alpha_0 <- tggp_param$alpha_0
  
  # Initial values diagnostics ---------------------------------------------
  initial_values$log_likelihood_Z <- compute_log_likelihood_Z(
    rows_obs_edges, 
    cols_obs_edges,
    rows_nobs_edges,
    cols_nobs_edges,
    w,
    Pi
  )

  # Define HMC parameters --------------------------------------------------
  epsilon <- mcmc_settings$epsilon
  L <- mcmc_settings$L
  
  # Prepare storage --------------------------------------------------------
  # Samples
  mcmc_samples <- list(
    w = matrix(NA, nrow = mcmc_settings$n_samples, ncol = N), 
    Pi = array(NA, c(mcmc_settings$n_samples, N, K)),
    w_rem = rep(NA, mcmc_settings$n_samples), 
    alpha = rep(NA, mcmc_settings$n_samples),
    sigma = rep(NA, mcmc_settings$n_samples), 
    tau = rep(NA, mcmc_settings$n_samples),
    beta = matrix(NA, 
                  nrow = mcmc_settings$n_samples, 
                  ncol = K)
  )   
  # Acceptance
  mcmc_acceptance <- list(
    rate_HMC = rep(NA, mcmc_settings$n_iter), 
    rate_MH_hyper = rep(NA, mcmc_settings$n_iter),
    epsilon = rep(NA, mcmc_settings$n_iter),
    log_posterior_ratio_w = rep(NA, mcmc_settings$n_iter)
  )
  # Diagnostics
  mcmc_diagnostics <- list(
    log_likelihood_Z = rep(NA, mcmc_settings$n_iter),
    log_posterior_w = rep(NA, mcmc_settings$n_iter),
    log_posterior_Pi = rep(NA, mcmc_settings$n_iter)
  )
  
  # Run MCMC ---------------------------------------------------------------
  for (i in 1:mcmc_settings$n_iter) {
    
    # Communicate the state of the MCMC, if wanted
    if (general_settings$verbose == TRUE && 
        (i %% general_settings$speak_every == 0)) {
      cat("----------------------------------------\n")
      cat(paste("Iteration:", i, "\n"))
      cat(paste("alpha:", round(alpha, 2), "\n"))
      cat(paste("sigma:", round(sigma, 2), "\n"))
      cat(paste("tau:", round(tau, 2), "\n"))
      cat(paste0("beta: ", "(", paste(round(beta_vec, 2), 
                                      collapse = ", "), ")", "\n"))
      if (i > 1) {
        cat(paste("rate_HMC:", round(mcmc_acceptance$rate_HMC[i-1], 2), "\n"))
        cat(paste("rate_MC_hyper:", round(mcmc_acceptance$rate_MH_hyper[i-1], 2), "\n")) 
      }
      cat("----------------------------------------")
    }
    
    # Update number of edges and their community assignments
    M_counts <- update_M_counts(rows_obs_edges,
                                cols_obs_edges,
                                w,
                                Pi,
                                w_nobs,
                                Pi_nobs)
    
    # Update community memberships
    beta_Pi <- update_beta_Pi(M_counts, beta_vec, gamma_0, alpha_0)
    beta_vec <- beta_Pi$beta_vec
    Pi <- beta_Pi$Pi_obs_nobs[1:N,]
    Pi_nobs <- beta_Pi$Pi_obs_nobs[(N+1):N_tot, ]
    
    # Update sociabilities
    N_counts <- rowSums(M_counts)
    w_updated <- update_w(c(w, w_nobs), w_rem, sigma, tau, N_counts, L, epsilon)
    w <- w_updated$w[1:N]
    w_nobs <- w_updated$w[(N+1):N_tot]
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
    
    # Update GGP hyperparameters
    if (tggp_param$estimate_alpha | 
        tggp_param$estimate_sigma | 
        tggp_param$estimate_tau) {
      
      # Update total sociability of removed nodes and hyperparams using MH
      # (alternate between gamma and random walk for alpha)
      if (i %% 2 == 0) {
        rw_alpha <- FALSE
      } else {
        rw_alpha <- TRUE
      }    
      
      ggp_param_updated <- update_ggp_hyperparam(c(w, w_nobs), 
                                                 w_rem, 
                                                 alpha, 
                                                 sigma, 
                                                 tau,
                                                 tggp_param, 
                                                 rw_alpha = rw_alpha)
      w_rem <- ggp_param_updated$w_rem
      alpha <- ggp_param_updated$alpha
      sigma <- ggp_param_updated$sigma
      tau <- ggp_param_updated$tau
      mcmc_acceptance$rate_MH_hyper[i] <- ggp_param_updated$rate_MH
      
    }
    
    # Compute diagnostics
    mcmc_diagnostics$log_likelihood_Z[i] <- compute_log_likelihood_Z(
      rows_obs_edges, 
      cols_obs_edges,
      rows_nobs_edges,
      cols_nobs_edges,
      w,
      Pi
    )
    mcmc_diagnostics$log_posterior_w[i] <- compute_log_posterior_w(N_counts, 
                                                                   c(w, w_nobs), 
                                                                   w_rem, 
                                                                   sigma, 
                                                                   tau)
    mcmc_diagnostics$log_posterior_Pi[i] <-  compute_log_posterior_Pi(M_counts, 
                                                                      Pi, 
                                                                      beta_vec, 
                                                                      alpha_0)
    
    # Update total number of nodes (unthinned and thinned)
    if ((i < tggp_param$stop_resample_N_tot_after) & 
        (i %% tggp_param$resample_N_tot_every == 0)){
      
      N_tot_old <- N_tot
      N_tot_to_N_ratio <-  update_N_tot_to_N_ratio(alpha,
                                                   sigma,
                                                   tau,
                                                   beta_vec,
                                                   tggp_param$num_N_tot_samples)
      N_tot <- round(N * N_tot_to_N_ratio)  
      
      if (N_tot_old != N_tot) {
        # Trim or enlarge w_not and Pi_not
        w_Pi_nobs <- trim_enlarge_w_Pi_nobs(N_tot, N_tot_old, N, w_nobs, Pi_nobs, beta_vec)
        w_nobs <- w_Pi_nobs$w_nobs
        Pi_nobs <- w_Pi_nobs$Pi_nobs 
      }
      
    }
    
    # Store samples
    if (i %in% mcmc_settings$iters_to_save) {
      
      ii <- which(mcmc_settings$iters_to_save == i)
      mcmc_samples$w[ii, ] <- w
      mcmc_samples$Pi[ii, , ] <- Pi
      mcmc_samples$w_rem[ii] <- w_rem
      mcmc_samples$alpha[ii] <- alpha
      mcmc_samples$sigma[ii] <- sigma
      mcmc_samples$tau[ii] <- tau
      mcmc_samples$beta[ii, ] <- beta_vec
      
    }
  }
  
  return(list(
    mcmc_samples = mcmc_samples,
    mcmc_acceptance = mcmc_acceptance,
    mcmc_diagnostics = mcmc_diagnostics
  ))
}




