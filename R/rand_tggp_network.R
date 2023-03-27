#' Sample of a network from the Thinned Generalized gamma_0 Process 
#'
#' ADD DESCRIPTION HERE
#' 
#' @param alpha positive real number. Affects the number of nodes and edges (higher alpha leads to larger networks).
#' @param sigma real number. in (-Inf, 1). sigma < 0 corresponds to a finite-activity GGP (dense case).
#' @param tau positive real number. Affects the tails of the degree distribution (higher tau leads to faster decay).
#' @param K positive integer. Number of communities that nodes can have memberships in.
#' @param beta_vec positive real-valued vector. beta_vec[k] affects the overall size of community k. Defaults to a vector with all elements equal to 1/K. If K is specified, length of beta_vec must be equal to K. 
#' @param gamma_0 real number. If gamma_0 is specified, beta_vec is sampled from a K-dimensional Dirichlet with concentration parameters gamma_0/K. If both gamma_0 and beta_vec are provided, provided gamma_0 is ignored.
#' @param zeta real number. Node community memberships are sampled from a K-dimensional Dirichlet with concentration parameter (zeta*beta_vec)/K. Defaults to 0.5.
#' @param return_thinned logical. Defaults to FALSE. If TRUE, thinned edges, as well as sociabilities and community memberships of nodes with no unthinned edge, are also returned. 
#' @param only_simulate_N logical. Defaults to FALSE. If TRUE, return only number of unthinned nodes and of thinned nodes.
#' 
#' @return
#' \describe{
#' A list with at least the following components: 
#' \item{edges}{a matrix storing source node, target node, and community 
#'   assignment of source and target nodes for each unthinned edge.}
#' \item{Z}{sparse symmetric binary adjacency matrix of undirected network.}
#' \item{nodes}{a matrix storing sociability and community memberships of 
#'   all nodes with at least one unthinned edge.}
#' If ``return_thinned=TRUE``, the following additional components are also 
#' included in the returned list:
#' \item{edges_thinned}{a matrix storing source node, target node, and community 
#'   assignments of source and target nodes for each thinned edge.}
#' \item{Z_thinned}{sparse symmetric binary adjacency matrix of thinned edges.}
#' \item{nodes_thinned}{a matrix storing sociability and community memberships 
#'   of nodes with no unthinned edge.}
#' }
#' 
#' @examples
#' # Sample a TGGP network with 50 communities of similar size
#' tggp_1 <- rand_tggp_network(200, 0.5, 1, K = 50)
#'
#' 
#' @references 
#' F. Ricci, M. Guindani, and E. Sudderth (2022). Thinned Random Measures for Sparse Graphs with Overlapping Communities. Neural Information Processing Systems 
#' 
#' Caron, F., & Fox, E. B. (2017). Sparse graphs using exchangeable random measures. Journal of the Royal Statistical Society. Series B, Statistical Methodology, 79(5), 12
#' 
#' @export

rand_tggp_network <- function(alpha, 
                              sigma, 
                              tau, 
                              K = 2, 
                              beta_vec = NULL, 
                              gamma_0 = NULL, 
                              zeta = NULL,
                              return_thinned = FALSE,
                              only_simulate_N = FALSE) {
  
  # Prepare list to return --------------------------------------------------
  tggp_network <- list()
  
  # Sample potential nodes --------------------------------------------------
  w_pot <- rand_ggp(alpha, sigma, tau)
  N_pot <- length(w_pot)
  sum_w_pot <- sum(w_pot)
  
  # Sample potential edges sources (i) and targets (j) ----------------------
  D_star <- rpois(1, sum_w_pot^2)
  nodes_i_j <- sample.int(N_pot, (2 * D_star), replace = TRUE, prob = w_pot)
  
  # Remove self-edges -------------------------------------------------------
  nodes_i <- nodes_i_j[1:D_star]
  nodes_j <- nodes_i_j[(D_star + 1) : (2 * D_star)]
  nonself_edges <- (nodes_i != nodes_j)
  nodes_i <- nodes_i[nonself_edges]
  nodes_j <- nodes_j[nonself_edges]
  
  # Rearrange pairs so that nodes_i[e] > nodes_j[e] -------------------------
  to_rearrange <- nodes_i < nodes_j
  nodes_i_pre_rearrange <- nodes_i 
  nodes_j_pre_rearrange <- nodes_j
  nodes_i[to_rearrange] <- nodes_j_pre_rearrange[to_rearrange]
  nodes_j[to_rearrange] <- nodes_i_pre_rearrange[to_rearrange]
  nodes_i_j <- c(nodes_i, nodes_j)
  N_e <- length(nodes_i)
  
  if (K >= 2) {
    
    # Sample relative community sizes -----------------------------------------
    if (!is.null(beta_vec)) { # beta_vec was provided
      
      # Check that beta_vec is admissible
      if (length(beta_vec) != K | !(all(beta_vec > 0))) {
        stop("beta_vec must be a positive real-valued vector of length K")
      } 
    } else if (!is.null(gamma_0)) { # beta_vec was not provided, but gamma_0 was provided
      gamma_0_aux_rv <- rgamma_0(K, shape = gamma_0/K, scale = 1)
      beta_vec <- gamma_0_aux_rv / sum(gamma_0_aux_rv)
    } else { # neither beta_vec nor gamma_0 were provided
      beta_vec <- rep(1/K, K)
    }
    
    # Sample potential nodes' memberships -------------------------------------
    if (is.null(zeta)) {
      zeta <- 0.5
    }
    shapes <- beta_vec * zeta
    K <- length(shapes)
    Pi <- rand_dirichlet_singleconc(shapes, N_pot)
    
    # Assign nodes to communities for each of their potential edge ------------
    nodes_i_j_split <- split(seq_along(nodes_i_j), nodes_i_j)
    edges_presplit_order <- unlist(nodes_i_j_split, use.names = FALSE)
    N_i_j_sorted <- lengths(nodes_i_j_split, use.names = FALSE) 
    nodes_i_j_idx <- as.integer(names(nodes_i_j_split))
    kappas_proposed <- rand_communities(K, nodes_i_j_idx, 
                                        N_i_j_sorted, edges_presplit_order, 
                                        Pi[nodes_i_j_idx,])
    kappas_proposed_i <- kappas_proposed[1:N_e] 
    kappas_proposed_j <- kappas_proposed[(N_e + 1) : (2 * N_e)]
    
    unthinned_edges <- (kappas_proposed_i == kappas_proposed_j)
    thinned_edges <- (kappas_proposed_i != kappas_proposed_j)
    
    if (only_simulate_N) {
      nodes_i_j_unth <- c(nodes_i[unthinned_edges],
                          nodes_j[unthinned_edges])
      nodes_i_j_th <- c(nodes_i[thinned_edges],
                        nodes_j[thinned_edges])
      N_unthinned <- length(unique(nodes_i_j_unth))
      N_thinned <- length(setdiff(nodes_i_j_th, 
                                  nodes_i_j_unth))
      return(list(
        N_unthinned = N_unthinned,
        N_thinned = N_thinned
      ))
    }
    
    # Separate unthinned edges ------------------------------------------------
    N_unthinned_edges <- sum(unthinned_edges)
    nodes_i_unth <- nodes_i[unthinned_edges]
    nodes_j_unth <- nodes_j[unthinned_edges]
    kappa_unth <- kappas_proposed_i[unthinned_edges]

    # Reindex unthinned nodes -------------------------------------------------
    nodes_i_j_unth_reindex <- dplyr::dense_rank(
      c(nodes_i_unth, nodes_j_unth)
    ) 
    N <- length(unique(nodes_i_j_unth_reindex))
    nodes_i_unth_reindex <- nodes_i_j_unth_reindex[1:N_unthinned_edges]
    nodes_j_unth_reindex <- nodes_i_j_unth_reindex[(N_unthinned_edges+1):(N_unthinned_edges*2)]
    
    # Prepare unthinned edges and nodes summaries -----------------------------
    tggp_network$Z <- Matrix::sparseMatrix(
      i = nodes_i_unth_reindex - 1,
      j = nodes_j_unth_reindex - 1,
      x = 1,
      dims = c(N, N), 
      use.last.ij = TRUE,
      index1 = FALSE,
      repr = "T"
    )
    tggp_network$edges <- matrix(c(tggp_network$Z@i + 1, 
                                   tggp_network$Z@j + 1,
                                   kappa_unth), 
                                 ncol = 3, 
                                 byrow = FALSE)
    colnames(tggp_network$edges) <- c("source", "target", "community")

    unth_nodes <- sort(unique(c(nodes_i[unthinned_edges], 
                                nodes_j[unthinned_edges])))
    Pi_unth <- Pi[unth_nodes, ]
    w_unth <- w_pot[unth_nodes]
    tggp_network$nodes <- cbind(1:N, w_unth, Pi_unth)
    colnames(tggp_network$nodes) <- c("v", "w", paste0("pi_", 1:K))
    
    # Prepare thinned edges and nodes summaries -------------------------------
    if (return_thinned) {
      
      # Reindex thinned edges  -----------------------------------------
      nodes_i_th <- nodes_i[thinned_edges]
      nodes_j_th <- nodes_j[thinned_edges]  
      kappa_i_th <- kappas_proposed_i[thinned_edges]
      kappa_j_th <- kappas_proposed_j[thinned_edges]
      thinned_nodes <- sort(setdiff(c(nodes_i_th, nodes_j_th), 
                                    c(nodes_i_unth, nodes_j_unth)))
      N_th_nodes <- length(thinned_nodes)
      thinned_nodes_reindx <- (N + 1):(N + N_th_nodes) 
      nodes_i_j_th_reindex <- c(nodes_i_th, nodes_j_th)
      nodes_i_j_th_reindex[nodes_i_j_th_reindex %in% thinned_nodes] <- 
        dplyr::dense_rank(
          c(nodes_i_th[nodes_i_th %in% thinned_nodes],
            nodes_j_th[nodes_j_th %in% thinned_nodes])
        ) + N
      nodes_i_j_th_reindex[nodes_i_j_th_reindex %in% unth_nodes] <- 
        match(nodes_i_j_th_reindex[nodes_i_j_th_reindex %in% unth_nodes],
              unth_nodes)
      # Prepare thinned edges and nodes summaries ---------------------
      tggp_network$edges_thinned <- matrix(c(nodes_i_j_th_reindex,
                                             kappa_i_th,
                                             kappa_j_th), 
                                           ncol = 4, byrow = FALSE)
      colnames(tggp_network$edges_thinned) <- c("source", 
                                                "target",
                                                "community_source",
                                                "community_target")
      tggp_network$Z_thinned <- new(
        "dgTMatrix",
        i = as.integer(c(tggp_network$edges_thinned[, 1] - 1,
                         tggp_network$edges_thinned[, 2] - 1)),
        j = as.integer(c(tggp_network$edges_thinned[, 2] - 1,
                         tggp_network$edges_thinned[, 1] - 1)),
        x = rep(1, length(nodes_i_th) * 2),
        Dim = c(N + N_th_nodes, N + N_th_nodes))    
      Pi_th <- Pi[thinned_nodes, ]
      w_th <- w_pot[thinned_nodes]
      tggp_network$nodes_thinned <- cbind(thinned_nodes_reindx,
                                          w_th,
                                          Pi_th)
      colnames(tggp_network$nodes_thinned) <- c("v", "w", paste0("pi_", 1:K))
      
    }
  } else {
    # Prepare edges and nodes summaries --------------------------------
    nodes_i_j_reindex <- dplyr::dense_rank(
      nodes_i_j
    ) 
    N_edges <- length(nodes_i_j_reindex)/2
    nodes_i_unth_reindex <- nodes_i_j_reindex[1:N_edges]
    nodes_j_unth_reindex <- nodes_i_j_reindex[(N_edges+1):(N_edges*2)]
    N <- length(unique(nodes_i_j_reindex))
    tggp_network$edges <- matrix(c(nodes_i_j_reindex, 
                                   rep(1, N_edges)), 
                                 ncol = 3, byrow = FALSE)
    colnames(tggp_network$edges) <- c("source", "target", "community")
    tggp_network$Z <- Matrix::sparseMatrix(
      i = nodes_i_unth_reindex - 1,
      j = nodes_j_unth_reindex - 1,
      x = 1,
      dims = c(N, N), 
      use.last.ij = TRUE,
      index1 = FALSE,
      repr = "T"
    )
    nodes <- sort(unique(c(nodes_i, 
                           nodes_j)))
    if (only_simulate_N) {
      stop("N_tot = N if K = 1")
    }
    w <- w_pot[nodes]
    tggp_network$nodes <- cbind(c(1:N), w, rep(1, N))
    colnames(tggp_network$nodes) <- c("v", "w", "pi_1")
    
  }
  
  
  return(tggp_network)
  
}
