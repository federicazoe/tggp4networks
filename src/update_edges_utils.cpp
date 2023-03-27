#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector update_auxiliary_tables(NumericMatrix M_counts,
                                      NumericVector beta_vec,
                                      double alpha_0) {
  
  int K = beta_vec.size();
  int N_star_0 = M_counts(_, 0).size();
  NumericVector auxiliary_tables(K);
  
  for (int h = 0; h < K; ++h) {
    for (int i = 0; i < N_star_0; ++i){
      
      if (M_counts(i, h) >= 1) {
        auxiliary_tables(h) += 1;
      }
      
      if (M_counts(i, h) >= 2) {
        int n = 1;
        for (int c = 0; c < M_counts(i, h); ++c) {
          auxiliary_tables(h) += R::rbinom(1, beta_vec(h)*alpha_0 / (n + beta_vec(h)*alpha_0));
          n += 1;
        }
      }
      
    }
  }
  
  return(auxiliary_tables);
  
}

// [[Rcpp::export]]
NumericMatrix sample_book_keep_concordant(IntegerVector n_ij, 
                                          NumericMatrix Pi_ij,
                                          NumericMatrix M_counts, 
                                          IntegerVector rows_obs_edges,
                                          IntegerVector cols_obs_edges) {
  
  int N_edges = Pi_ij(_, 0).size();
  int K = Pi_ij(0, _).size();
  
  for (int e = 0; e < N_edges; ++e) {
    
    int i = rows_obs_edges[e] - 1;
    int j = cols_obs_edges[e] - 1;
    int n_e = n_ij[e];
    NumericVector p_ij = Pi_ij(e, _);
    IntegerVector k_ij = Rcpp::sample(K, n_e, TRUE, p_ij, TRUE);
    
    for (int ee = 0; ee < n_e; ++ee) {
      
      int k = k_ij[ee] - 1;
      M_counts(i, k)++;
      M_counts(j, k)++;
      
    }
    
  }
  
  return(M_counts);
  
}

// [[Rcpp::export]]
NumericMatrix sample_multiple_n_and_p(IntegerVector n_ii, 
                                      NumericMatrix p_ii) {
  
  int K = p_ii(0, _).size();
  int N_ii = p_ii(_, 0).size();
  NumericMatrix K_ii(N_ii, K);
  
  for (int e = 0; e < N_ii; ++e) {
    
    int n_e = n_ii[e];
    NumericVector p = p_ii(e, _);
    
    IntegerVector k_ii = Rcpp::sample(K, n_e, TRUE, p, TRUE);
    
    for (int ee = 0; ee < n_e; ++ee) {
      
      int k = k_ii[ee] - 1;
      K_ii(e, k)++;
      
    }
    
  }
  
  return(K_ii);
  
}

// [[Rcpp::export]]
IntegerVector propose_communities(IntegerVector nodes_i_j_idx,
                                  IntegerVector N_i_j_sorted,
                                  IntegerVector edges_presplit_order,
                                  NumericMatrix Pi){
  
  int K = Pi(0, _).size();
  int N_edges_2 = edges_presplit_order.size();
  IntegerVector kappas_proposed(N_edges_2);    
  int e_plc_hld = 0;
  
  int num_nodes = nodes_i_j_idx.size();
  
  for (int v = 0; v < num_nodes; ++v) {
    
    int idx_v = nodes_i_j_idx[v] - 1;
    int n_v = N_i_j_sorted[idx_v];
    
    if (n_v > 0) {
      
      NumericVector pi_v = Pi(idx_v, _);
      IntegerVector k_v = Rcpp::sample(K, n_v, TRUE, pi_v, TRUE);
      
      for (int t = 0; t < n_v; ++t){
        
        int idx_e = edges_presplit_order[e_plc_hld] - 1;
        kappas_proposed[idx_e] = k_v[t];
        ++e_plc_hld;
      }  
    }
  }
  
  return(kappas_proposed);
  
}

// [[Rcpp::export]]
NumericMatrix book_keep_selected_edges(NumericMatrix M_counts, 
                                       IntegerVector nodes_i_j_selected,
                                       IntegerVector kappas_selected){
  
  int N_edges_2 = nodes_i_j_selected.size();
  
  for (int e = 0; e < N_edges_2; ++e) {
    
    int idx_v = nodes_i_j_selected[e] - 1;
    int k_e = kappas_selected[e] - 1;
    ++M_counts(idx_v, k_e);
    
  }
  
  return(M_counts);
  
}

// [[Rcpp::export]]
NumericMatrix prepare_R_obs(int N_nobs,
                            int N,
                            NumericVector w_obs,
                            NumericMatrix Pi_obs,
                            NumericMatrix Pi_nobs){
  
  NumericMatrix R_obs(N_nobs, N);
  int K = Pi_obs(0, _).size();
  
  for (int j = 0; j < N; ++j) {
    
    for (int i = 0; i < N_nobs; ++i) {
      
      double sum_pi_i_pi_j_not = 0;
      for (int k = 0; k < K; ++k) {
        sum_pi_i_pi_j_not = sum_pi_i_pi_j_not + Pi_obs(j, k) * Pi_nobs(i, k);
      }
      
      R_obs(i, j) = 2 * w_obs[j] * (1 - sum_pi_i_pi_j_not);
      
    }
  }
  
  return(R_obs);
  
}

// [[Rcpp::export]]
NumericMatrix prepare_R_unobs(int N_nobs,
                              NumericVector w_nobs,
                              NumericMatrix Pi_nobs){
  
  NumericMatrix R_unobs(N_nobs, N_nobs);
  int K = Pi_nobs(0, _).size();
  
  for (int i = 0; i < (N_nobs - 1); ++i) {
    
    R_unobs(i, i) = w_nobs[i];
    
    for (int j = i + 1; j < N_nobs; ++j) {
      
      double sum_pi_i_pi_j_not = 0;
      for (int k = 0; k < K; ++k) {
        sum_pi_i_pi_j_not = sum_pi_i_pi_j_not + Pi_nobs(j, k) * Pi_nobs(i, k);
      }
      
      R_unobs(i, j) = 2 * w_nobs[j] * (1 - sum_pi_i_pi_j_not);
      
    }
  }
  
  R_unobs(N_nobs-1, N_nobs-1) = w_nobs[N_nobs-1];
  
  return(R_unobs);
  
}

// [[Rcpp::export]]
NumericMatrix sample_unobs_nodes_edges(int N_tot, int N_nobs, int K,
                                       NumericVector tot_R, 
                                       IntegerVector joined_network,
                                       NumericMatrix R_obs_unobs,
                                       NumericMatrix Pi_all,
                                       NumericMatrix M_counts){
  
  Function rztp("rztpois"); 
  IntegerVector N_i(N_nobs);
  int n_i_not;
  IntegerVector x_2;
  int N = N_tot - N_nobs;
  
  for (int i = 0; i < N_nobs; ++i) {
    
    int idx_i = N + i;
    
    if (joined_network[idx_i] == 0) {
      n_i_not = as<int>(rztp(tot_R[i]));
    } else {
      n_i_not = R::rpois(tot_R[i]);
    }
    
    if (n_i_not > 0) {
      
      NumericVector p = R_obs_unobs(i, _);
      x_2 = Rcpp::sample(N_tot, n_i_not, TRUE, p, TRUE);
      
      IntegerVector counts(N_tot);
      
      for (int e = 0; e < n_i_not; ++e) {
        
        ++counts[x_2[e] - 1];
        
      }
      
      NumericVector pi_i = Pi_all(idx_i, _);
      
      for (int j = 0; j < N_tot; ++j) {
        
        int n_j = counts[j];
        
        if (n_j > 0) {
          
          joined_network[j] = 1;
          
          if (j != idx_i) {
            
            NumericVector pi_j = Pi_all(j, _);
            NumericVector pi_i_adj = pi_i * (1 - pi_j);
            
            IntegerVector k_i = Rcpp::sample(K, n_j, TRUE, pi_i_adj, TRUE);
            
            IntegerVector counts_j(K);
            
            for (int e = 0; e < n_j; ++e) {  
              
              int k = k_i[e] - 1;
              M_counts(idx_i, k) = M_counts(idx_i, k) + 1;  
              ++counts_j[k];
              
            }
            
            for (int k = 0; k < K; ++k) {
              
              int n_j_k = counts_j[k];
              
              if (n_j_k > 0) {
                
                NumericVector pi_j_adj(K);
                pi_j_adj = Rcpp::clone(pi_j);
                pi_j_adj[k] = 0;
                IntegerVector k_j = Rcpp::sample(K, n_j_k, TRUE,  pi_j_adj, TRUE); 
                
                for (int ee = 0; ee < n_j_k; ++ee) {
                  
                  int k2 = k_j[ee] - 1;
                  M_counts(j, k2) = M_counts(j, k2) + 1;
                  
                }
              }
            }
          } else {
            
            int n_j_double = 2*n_j;
            IntegerVector k_i = Rcpp::sample(K, n_j_double, TRUE, pi_i, TRUE);
            
            for (int e = 0; e < n_j_double; ++e) {  
              
              int k = k_i[e] - 1;
              M_counts(idx_i, k) = M_counts(idx_i, k) + 1;  
              
            }
          }
        }
      }
    }
  }
  
  return(M_counts);
  
}



// [[Rcpp::export]]
IntegerVector book_keep_edges_counts_ggp(IntegerVector N_ij, 
                                         IntegerVector rows_obs_edges,
                                         IntegerVector cols_obs_edges,
                                         IntegerVector N_counts){
  int E = rows_obs_edges.size();
  for (int e = 0; e < E; ++e) {
    N_counts[rows_obs_edges[e] - 1] += N_ij[e];
    N_counts[cols_obs_edges[e] - 1] += N_ij[e];
  }
  return(N_counts);
}

