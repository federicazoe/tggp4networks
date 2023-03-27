#' Samples number of edges that each node participates in and their community assignments, given current nodes'sociabilities and community memberships
#' @param rows_obs_edges: vector of integers; row indices of 1s in lower triangular binary adjacency matrix of undirected network
#' @param cols_obs_edges: vector of integers; column indices of 1s in lower triangular binary adjacency matrix of undirected network
#' @param w: vector of scalars: nodes sociabilities at which likelihood is evaluated 
#' @param Pi: matrix of probability vectors: nodes community memberships at which likelihood is evaluated
#' @param w_nobs: vector of scalars: nodes sociabilities at which likelihood is evaluated (unobserved nodes)
#' @param Pi_nobs: matrix of probability vectors: nodes community memberships at which likelihood is evaluated (unobserved nodes)
#' @keywords internal
update_M_counts <- function(rows_obs_edges, 
                            cols_obs_edges,
                            w,
                            Pi,
                            w_nobs,
                            Pi_nobs) {
  
  # Declare dimensions
  N <- length(w)
  N_nobs <- length(w_nobs)
  N_tot <- N + N_nobs
  K <- dim(Pi)[2]
  
  # Prepare object to store samples
  M_counts <- matrix(0L, nrow = N_tot, ncol = K)
  
  # Sample edges and community assignments corresponding to observed edges
  M_counts <- update_within_community_edges(rows_obs_edges, 
                                            cols_obs_edges,
                                            w, 
                                            Pi, 
                                            M_counts)
  
  # Sample edges and community assignments that have been thinned (observed nodes)
  M_counts <- update_across_community_edges(w, 
                                            Pi, 
                                            M_counts)  

  # Sample edges and community assignments that have been thinned (unobserved nodes)
  if (N_nobs >= 1) {
    M_counts <- update_across_community_edges_unobs_nodes(w, 
                                                          Pi, 
                                                          w_nobs, 
                                                          Pi_nobs, 
                                                          M_counts)
  }
  
  
  
}


#' Samples within-community edges in the latent multigraph
#' @param rows_obs_edges: vector of integers; row indices of 1s in lower triangular binary adjacency matrix of undirected network
#' @param cols_obs_edges: vector of integers; column indices of 1s in lower triangular binary adjacency matrix of undirected network
#' @param w: vector of scalars: nodes sociabilities at which likelihood is evaluated 
#' @param Pi: matrix of probability vectors: nodes community memberships at which likelihood is evaluated
#' @param M_counts: matrix indicating number of edges each node has in each community
#' @keywords internal
update_within_community_edges <- function(rows_obs_edges, 
                                          cols_obs_edges,
                                          w, 
                                          Pi, 
                                          M_counts) {
  # Declare dimensions
  N <- length(w)
  K <- dim(Pi)[2]
  
  # Compute summary statistics
  w_i <- w[rows_obs_edges]
  w_j <- w[cols_obs_edges]
  two_w_i_w_j <- 2 * w_i * w_j
  Pi_ij <- Pi[rows_obs_edges, ] * Pi[cols_obs_edges, ]
  Pi_ij_sum <- rowSums(Pi_ij)

  # Sample total number of edges in the latent multigraph between all pairs of 
  # nodes with an observed edge in the observed simple graph
  n_ij <- rztpois_vectorized(two_w_i_w_j * Pi_ij_sum)
  
  # Sample community assignments for each edge (using Rcpp function)
  M_counts <- sample_book_keep_concordant(n_ij, 
                                          Pi_ij,
                                          M_counts,
                                          rows_obs_edges,
                                          cols_obs_edges)
  
  return (M_counts)

}

#' Samples across-community edges in the latent multigraph between observed nodes
#' @param w: vector of scalars: nodes sociabilities at which likelihood is evaluated 
#' @param Pi: matrix of probability vectors: nodes community memberships at which likelihood is evaluated
#' @param M_counts: matrix indicating number of edges each node has in each community
#' @keywords internal
update_across_community_edges <- function(w, Pi, M_counts) {
  
  # Declare dimensions
  N <- length(w)
  K <- dim(Pi)[2]
  
  # Propose edges only based on nodes sociabilities
  N_edges <- rpois(1, sum(w)^2)
  nodes_idx <- sample.int(N, 
                          N_edges*2,
                          replace = TRUE, 
                          prob = w)  
  rows_prop_edges <- nodes_idx[1:N_edges]
  cols_prop_edges <- nodes_idx[(N_edges + 1) : (2 * N_edges)]

  # SELF-EDGES -------------------------------------------------------------
  # All proposed self-edges are kept, regardless of their community assignment
  # (they would be unobserved even if they were within-community)
  which_self_edges <- (rows_prop_edges == cols_prop_edges)
  if (sum(which_self_edges) > 0) {
    
    # Extract what nodes have self edges, how many self edges they have
    self_edges <- rows_prop_edges[which_self_edges]
    n_ii <- tabulate(self_edges, N)
    ii_node_indices <- which(n_ii > 0)
    n_ii <- n_ii[ii_node_indices]
    p_ii <- Pi[ii_node_indices,]
    
    if (length(n_ii) > 1) {
      
      # Sample community assignments of all self edges
      K_ii_counts <- sample_multiple_n_and_p(n_ii, p_ii) 
      M_counts[ii_node_indices, ] <- M_counts[ii_node_indices, ] + 
        K_ii_counts
      
    } else {
 
      # Sample community assignment of the only self edge
      k <- sample.int(K, 2, prob = p_ii)
      M_counts[ii_node_indices, k[1]] <- M_counts[ii_node_indices, k[1]]+ 1
      M_counts[ii_node_indices, k[2]] <- M_counts[ii_node_indices, k[2]]+ 1
      
    }
    
  }

  # NON SELF-EDGES ---------------------------------------------------------
  # Proposed non-self edges are independently assigned to communities 
  # and are retained only if they are across-community (otherwise they
  # would have been observed)
  which_non_self_edges <- (rows_prop_edges != cols_prop_edges)
  if (sum(which_non_self_edges) > 0) {
 
    # Find out how many edges each node has been proposed for
    nodes_i <- rows_prop_edges[which_non_self_edges]
    nodes_j <- cols_prop_edges[which_non_self_edges]
    N_edges <- length(nodes_i)
    nodes_i_j <- c(nodes_i, nodes_j)
    nodes_i_j_split <- split(seq_along(nodes_i_j), nodes_i_j)
    edges_presplit_order <- unlist(nodes_i_j_split, use.names = FALSE)
    N_i_j_sorted <- tabulate(nodes_i_j, N)
    nodes_i_j_idx <- as.integer(names(nodes_i_j_split))
    
    # For each edge, propose a community for both nodes assigned to the edge
    kappas_proposed <- propose_communities(nodes_i_j_idx, 
                                           N_i_j_sorted, edges_presplit_order, 
                                           Pi)
    kappas_proposed_i <- kappas_proposed[1:N_edges] 
    kappas_proposed_j <- kappas_proposed[(N_edges + 1) : (2 * N_edges)] 
    
    # Select across-community edges
    edges_selected <- which(kappas_proposed_i != kappas_proposed_j)
    edges_selected_long <- c(edges_selected, edges_selected + N_edges)
    kappas_selected <- kappas_proposed[edges_selected_long]
    nodes_i_j_selected <- nodes_i_j[edges_selected_long]
    
    # Update edge counts accordingly
    M_counts <- book_keep_selected_edges(M_counts, 
                                         nodes_i_j_selected,
                                         kappas_selected)
  }
  
  return(M_counts)
  
}

#' Samples across-community edges in the latent multigraph involving unobserved nodes
#' @param w: vector of scalars: nodes sociabilities at which likelihood is evaluated 
#' @param Pi: matrix of probability vectors: nodes community memberships at which likelihood is evaluated
#' @param w_nobs: vector of scalars: nodes sociabilities at which likelihood is evaluated (unobserved nodes)
#' @param Pi_nobs: matrix of probability vectors: nodes community memberships at which likelihood is evaluated (unobserved nodes)
#' @param M_counts: matrix indicating number of edges each node has in each community
#' @keywords internal
update_across_community_edges_unobs_nodes <- function(w, Pi, w_nobs, Pi_nobs, M_counts) {
  
  # Declare dimensions
  N <- length(w)
  N_nobs <- length(w_nobs)
  N_tot <- N + N_nobs
  K <- dim(Pi)[2]
  
  # Rates with observed nodes
  # Prepare (N_nobs x N) matrix R_obs where R_obs[i, j] = 2 * w_j * (1 - ∑\pi_ih \pi_jh)
  # R_obs <- matrix(w, nrow = N_nobs, ncol = N, byrow = TRUE)
  # R_obs <- 2 * R_obs * (1 - Pi_nobs %*% t(Pi))
  R_obs <- prepare_R_obs(N_nobs, N, w, Pi, Pi_nobs)
  tot_R_obs <- rowSums(R_obs)
  
  # Rates with unobserved nodes
  # Pi_component <- 1 - Pi_nobs %*% t(Pi_nobs)
  # diag(Pi_component) <- 1 # self edges
  # R_unobs <- matrix(w_nobs, nrow = N_nobs, ncol = N_nobs, byrow = TRUE)
  # R_unobs <- R_unobs * upper.tri(R_unobs, diag = TRUE) # turn off lower triang. 
  # R_unobs <- 2 * R_unobs * Pi_component
  # diag(R_unobs) <- diag(R_unobs)/2 # self edges
  R_unobs <- prepare_R_unobs(N_nobs, w_nobs, Pi_nobs)
  tot_R_unobs <- rowSums(R_unobs)
  
  # Next, need to:
  # 1) Sample total edges
  # 2) Assign sampled edges 
  # 3) Assign to communities
  
  joined_network <- rep(0L, N_tot)
  tot_R <- w_nobs * (tot_R_obs + tot_R_unobs)
  R_obs_unobs <- cbind(R_obs, R_unobs)
  Pi_all <- rbind(Pi, Pi_nobs)
  
  M_counts <- sample_unobs_nodes_edges(N_tot, N_nobs, K, 
                                      tot_R, joined_network,
                                      R_obs_unobs, Pi_all, M_counts)
  
  return(M_counts)
  
}


#' Samples number of edges that each node participates in, given current nodes'sociabilities
#' @param rows_obs_edges: vector of integers; row indices of 1s in lower triangular binary adjacency matrix of undirected network
#' @param cols_obs_edges: vector of integers; column indices of 1s in lower triangular binary adjacency matrix of undirected network
#' @param w: vector of scalars: nodes sociabilities at which likelihood is evaluated 
#' @keywords internal
update_N_counts_ggp <- function(rows_obs_edges, 
                                cols_obs_edges,
                                w) {
  
  N <- length(w)
  N_counts <- rep(0L, N)
  
  # Non-self edges ---------------------------------------------------------
  two_w_i_w_j <- 2 * w[rows_obs_edges] * w[cols_obs_edges]
  N_ij <- rztpois_vectorized(two_w_i_w_j)
  E <- length(rows_obs_edges)
  N_counts <- book_keep_edges_counts_ggp(N_ij,
                                         rows_obs_edges,
                                         cols_obs_edges,
                                         N_counts)
  # Self edges -------------------------------------------------------------
  w_sq <- w^2
  N_counts <- N_counts + 2*rpois_vectorized(1, w_sq)
  
  return(N_counts)
  
}

#' Zero-truncated Poisson rv generation
#' 
#' Generates random variables from a zero-truncated Poisson distribution
#' using sequential conditioning
#' https://en.wikipedia.org/wiki/Zero-truncated_Poisson_distribution
#'
#' @param lambda: strictly positive real number. Letting X ~ Poisson(lambda), we
#'                want to generate a rv from P(X|X>0).
#'
#' @return t: draw from a Zero-truncated Poisson distribution         
#' @keywords internal
rztpois <- function(lambda) {
  
  t <- 0
  k <- 1
  
  while (t == 0) {
    
    pdf_conditional_k <- dpois(k, lambda = lambda) / ppois(k-1, lambda, 
                                                           lower.tail = FALSE)
    
    if (runif(1) <= pdf_conditional_k) {
      
      t <- k
      
    } else {
      
      k <- k + 1
      
    }
  }
  
  return(t)
  
}

rztpois_vectorized <- Vectorize(rztpois, "lambda")

rpois_vectorized <- Vectorize(rpois, "lambda")