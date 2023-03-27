#' Initializes all hyperparameters and parameters for fit_tggp_mcmc
#'
#' @param Z sparse symmetric binary adjacency matrix of undirected network.
#' @param n_iter_initial int; desired number of GGP mcmc runs
#' 
#' @keywords internal
initialize_tggp <- function(Z, 
                            n_iter_initial,
                            tggp_param,
                            verbose) {
  
  initial_values <- list()
  N <-  Z@Dim[1]
  
  if (verbose) {
    cat("Running GGP mcmc for initialization...\n")
    cat("... \n")
  }
  
  # Run GGP mcmc -----------------------------------------------------------
  ggp_initialized_param <- initialize_with_ggp(Z, 
                                               n_iter_initial)
  
  # Initialize GGP parameters ----------------------------------------------
  initial_values$alpha <- ifelse(is.null(tggp_param$alpha),
                                 ggp_initialized_param$alpha,
                                 tggp_param$alpha)
  initial_values$log_alpha <- log(initial_values$alpha)
  initial_values$sigma <- ifelse(is.null(tggp_param$sigma),
                                 ggp_initialized_param$sigma,
                                 tggp_param$sigma)
  initial_values$tau <- ifelse(is.null(tggp_param$tau),
                               ggp_initialized_param$tau,
                               tggp_param$tau)
  initial_values$w_rem <- ifelse(is.null(tggp_param$w_rem),
                                 ggp_initialized_param$w_rem,
                                 tggp_param$w_rem)
  
  # Initialize hyperparameteres for community memberships ------------------
  tggp_param <- initialize_model_param(tggp_param, N, model = "tggp")
  initial_values$beta <- tggp_param$beta
  
  # Initialize nodes sociabilities  ----------------------------------------
  if (!is.null(tggp_param$w) && length(tggp_param$w) ==  N) {
    initial_values$w <- tggp_param$w
  } else {
    initial_values$w <- ggp_initialized_param$w
  }
  
  # Initialize nodes community memberships ---------------------------------
  if (!is.null(tggp_param$Pi) && 
      (dim(tggp_param$Pi) ==  N) && 
      (dim(tggp_param$Pi)[2] == tggp_param$K_max)) {
    initial_values$Pi <- tggp_param$Pi
  } else {
    initial_values$Pi <- MCMCpack::rdirichlet(N, initial_values$beta)
  }
  
  # Initialize total number of nodes (unthinned and thinned) ---------------
  if (!is.null(tggp_param$N_tot) && tggp_param$N_tot > N) {
    initial_values$N_tot <- tggp_param$N_tot
  } else {
    beta_N_tot_init <-  MCMCpack::rdirichlet(
      1, 
      rep(tggp_param$gamma_0 / tggp_param$K_max, 
          tggp_param$K_max))
    N_tot_to_N_ratio <-  update_N_tot_to_N_ratio(
      initial_values$alpha,
      initial_values$sigma,
      initial_values$tau,
      initial_values$beta,
      tggp_param$num_N_tot_samples
    )
    initial_values$N_tot <- round(N * N_tot_to_N_ratio)  
  }
  
  # Initialize parameters of thinned nodes --------------------------------
  N_nobs <- initial_values$N_tot - N
  initial_values$w_nobs <- sample(initial_values$w,
                                  N_nobs, 
                                  replace = TRUE,
                                  prob = 1/(1 + initial_values$w))
  initial_values$Pi_nobs <- MCMCpack::rdirichlet(N_nobs, 
                                                 initial_values$beta)
  
  if (verbose) {
    cat("...Done with initialization \n")
    cat("Initialization values:\n")
    cat(paste("alpha: ", round(initial_values$alpha, 2), "\n"))
    cat(paste("sigma: ", round(initial_values$sigma, 2), "\n"))
    cat(paste("tau: ", round(initial_values$tau, 2), "\n"))
    cat(paste0("beta: ", "(", paste(round(initial_values$beta, 2), 
                                    collapse = ", "), ")", "\n"))
  }
  
  return(list(
    initial_values = initial_values,
    tggp_param = tggp_param
  ))
  
  
}

#' Initializes nodes sociabilities and GGP hyperparameters running GGP mcmc
#'
#' @param Z sparse symmetric binary adjacency matrix of undirected network.
#' @param n_iter_initial int; desired number of GGP mcmc runs
#' 
#' @keywords internal
initialize_with_ggp <- function(Z, 
                                n_iter_initial) {
  
  mcmc_settings_initial <- list(n_iter = n_iter_initial,
                                warmup = round(n_iter_initial/2))
  general_settings_initial <- list(verbose = FALSE)
  mcmc_samples <- fit_ggp_mcmc(Z,
                               mcmc_settings = mcmc_settings_initial, 
                               general_settings = general_settings_initial)$mcmc_samples
  
  ggp_initialized_param <- list()
  ggp_initialized_param$alpha <- mean(mcmc_samples$alpha)
  ggp_initialized_param$sigma <- mean(mcmc_samples$sigma)
  ggp_initialized_param$tau <- mean(mcmc_samples$tau)
  ggp_initialized_param$w_rem <- mean(mcmc_samples$w_rem)
  ggp_initialized_param$w <- apply(mcmc_samples$w, 2, mean)
  
  return(ggp_initialized_param)
  
}


#' Initializes model-related stuff (and for fit_ggp_mcmc, initialize hyperparams)
#'
#' @param model_param list; if ``model = tggp`` then the list needs to include a field K specifying the upper bound to the number of communities
#' @param N int; number of nodes in observed network
#' @param model string; defaults to "tggp". If given a different value, only hyperparameters from GGP model are initialized
#' 
#' @keywords internal
#'
initialize_model_param <- function(model_param, N, model = "tggp") {
  
  if (is.null(model_param$alpha)) {
    model_param$alpha <- 100 * runif(1)
  }
  
  model_param$log_alpha <- log(model_param$alpha)
  
  if (is.null(model_param$estimate_alpha)){
    model_param$estimate_alpha <- TRUE
  }
  
  if (is.null(model_param$sigma)) {
    model_param$sigma <- 2 * runif(1) - 1
  } 
  if (is.null(model_param$estimate_sigma)){
    model_param$estimate_sigma <- TRUE
  } 
  
  if (is.null(model_param$tau)) {
    model_param$tau <- 10 * runif(1)
  }
  if (is.null(model_param$estimate_tau)){
    model_param$estimate_tau <- TRUE
  }  
  
  if (is.null(model_param$w_rem)) {
    model_param$w_rem <- rgamma(1, 1, 1)
  }

  if (model == "tggp") {
    
    if (is.null(model_param$gamma_0)) {
      model_param$gamma_0 <- 0.2 * model_param$K_max
    } 
    if (is.null(model_param$alpha_0)) {
      model_param$alpha_0 <- 0.5
    }    
    
    if (is.null(model_param$beta)) {
      model_param$beta <- rep(1/model_param$K_max, model_param$K_max)
    }
    
    if (is.null(model_param$num_N_tot_samples)) {
      model_param$num_N_tot_samples <- 10
    }   
    
    if (is.null(model_param$resample_N_tot_every)) {
      model_param$resample_N_tot_every <- 100
    }
    
    if (is.null(model_param$stop_resample_N_tot_after)) {
      model_param$stop_resample_N_tot_after <- Inf
    }  
    
  }
  
  return(model_param)
}


#' Sets all MCMC settings
#'
#' @param mcmc_settings list; field "n_iter" must indicate desired number of iteration; other fields are optional.
#' @param N_tot int; number of nodes in observed network
#' @param model string; defaults to "tggp". If given a different value, only hyperparameters from GGP model are initialized
#' 
#' @keywords internal
#'
set_mcmc_settings <- function(mcmc_settings, N) {
  if (is.null(mcmc_settings$L)) {
    mcmc_settings$L <- 5
  }
  if (is.null(mcmc_settings$epsilon)) {
    mcmc_settings$epsilon <- 0.1
  }  
  if (is.null(mcmc_settings$eps_adjust_rate)){
    mcmc_settings$eps_adjust_rate <- 0.01
  }
  if (is.null(mcmc_settings$eps_avg_window)) {
    mcmc_settings$eps_avg_window <- 50
  }
  if (is.null(mcmc_settings$n_adapt)) {
    mcmc_settings$n_adapt <- mcmc_settings$n_iter / 4
  }
  
  if (is.null(mcmc_settings$warmup)) {
    mcmc_settings$warmup <- 0
  }
  if (is.null(mcmc_settings$n_adapt)) {
    mcmc_settings$n_adapt = floor(mcmc_settings$n_iter / 4)
  }
  
  if (is.null(mcmc_settings$thin)) {
    mcmc_settings$thin <- 1
  }
  
  if (is.null(mcmc_settings$n_iter_initial)) {
    mcmc_settings$n_iter_initial <- 10000
  }  
  
  save_iter <- which((mcmc_settings$warmup+1):mcmc_settings$n_iter %% mcmc_settings$thin == 0)
  mcmc_settings$iters_to_save <- c((mcmc_settings$warmup+1):mcmc_settings$n_iter)[save_iter]
  mcmc_settings$n_samples <- length(mcmc_settings$iters_to_save)
  mcmc_settings$epsilon <- mcmc_settings$epsilon / N^(1/4)

  return(mcmc_settings)
}

#' Sets all general settings
#'
#' @param general_settings list; field "path_for_saving" must indicate desired number of iteration; other fields are optional.
#' @param n_iter int; number of mcmc iterations
#' 
#' @keywords internal
#'
set_general_settings <- function(general_settings, n_iter) {
  if (is.null(general_settings$verbose)) {
    general_settings$verbose <- TRUE
  }  
  if (general_settings$verbose == TRUE) {
    if (is.null(general_settings$speak_every)) {
      general_settings$speak_every <- 50
    }
  }
  if (is.null(general_settings$save_every)){
    general_settings$save_every <- n_iter
  } else {
    if (general_settings$save_every > n_iter | 
        !is.numeric(general_settings$save_every) | general_settings$save_every <= 0) {
      stop("If specified, parameter save_every must be a positive number smaller than n_iter")
    }
  }
  return(general_settings)
}