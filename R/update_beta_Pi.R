
update_beta_Pi <- function(M_counts, beta_vec, gamma_0, alpha_0) {

  # Declare dimensions
  K <- dim(M_counts)[2]

  # Update auxiliary variables
  auxiliary_tables <- update_auxiliary_tables(M_counts, beta_vec, alpha_0)
  
  # Update vector of relative community frequencies
  beta_vec <- c(MCMCpack::rdirichlet(1, auxiliary_tables + rep(gamma_0/K, K)))
  
  # Update matrix of community memberships
  updated_alpha_beta <- beta_vec * alpha_0
  concentrations <- updated_alpha_beta + M_counts
  Pi_obs_nobs <- rand_dirichlet_multiconc(concentrations)
  
  return(
    list(
      beta_vec = beta_vec,
      Pi_obs_nobs = Pi_obs_nobs
    )
  )
  
}