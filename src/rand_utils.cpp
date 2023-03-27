#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector rand_communities(int K, 
                               IntegerVector nodes_i_j_idx,
                               IntegerVector N_i_j_sorted,
                               IntegerVector edges_presplit_order,
                               NumericMatrix Pi_i_j){
  
  int N_edges_2 = edges_presplit_order.size();
  IntegerVector communities(N_edges_2);    
  int e_plc_hld = 0;
  
  int num_nodes = nodes_i_j_idx.size();
  
  for (int v = 0; v < num_nodes; ++v) {
    
    NumericVector pi_v = Pi_i_j(v, _);
    int n_v = N_i_j_sorted[v];
    IntegerVector k_v = Rcpp::sample(K, n_v, TRUE, pi_v, TRUE);
    
    for (int t = 0; t < n_v; ++t){
      
      int idx_e = edges_presplit_order[e_plc_hld] - 1;
      communities[idx_e] = k_v[t];
      ++e_plc_hld;
    }  
  }
  
  return(communities);
  
}

//' Draws `n_samples` random variates from a Dirichlet distribution with
//' `concentrations`= (concentrations[1], ..., concentrations[K])
// [[Rcpp::export]]
NumericMatrix rand_dirichlet_singleconc(NumericVector concentrations,
                                         int n_samples) {
   
   int K = concentrations.size();
   NumericMatrix Pi(n_samples, K);
   
   for (int k = 0; k < K; ++k) {
     Pi(_, k) = Rcpp::rgamma(n_samples, concentrations[k], 1);
   }
   
   NumericVector summed_by_row = rowSums(Pi);
   
   for (int i = 0; i < n_samples; ++i) {
     Pi(i, _) = Pi(i, _) / summed_by_row[i];
   }
  
   return Pi;
   
}

//' Provided with a KxN matrix `concentrations`, draws one random variate from
//' each i of N Dirichlet distributions with concentration parameters
//' (concentrations[1, i], ..., concentrations[K, i])
// [[Rcpp::export]]
NumericMatrix rand_dirichlet_multiconc(NumericMatrix concentrations) {
   
   int K = concentrations(0,_).size();
   int N = concentrations(_,0).size();
   NumericMatrix Pi(N, K);
   
   for (int i = 0; i < N; ++i) {
     for (int k = 0; k < K; ++k) {
       Pi(i, k) = R::rgamma(concentrations(i, k), 1.0);
     }
     
     Pi(i, _) = Pi(i, _) / sum(Pi(i, _));
     
   }
   
   return Pi;
   
}
