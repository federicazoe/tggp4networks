#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double compute_log_posterior_Pi_cpp(NumericMatrix Pi,
                                    NumericMatrix beta_vec) {
  
  int N_star_0 = Pi(_, 0).size();
  
  double log_lik_Pi_beta = 0;
  
  for (int i = 0; i < N_star_0; ++i) {
    
    double logD = sum(lgamma(beta_vec(i, _) + 1e-300)) - lgamma(sum(beta_vec(i, _)));
    double s = sum((beta_vec(i, _) - 1)  * log(Pi(i, _) + 1e-300));
    log_lik_Pi_beta = log_lik_Pi_beta + s - logD;
    
  }
  
  return(log_lik_Pi_beta);
  
}

// [[Rcpp::export]]
NumericVector compute_entropy_Pi_cpp(NumericMatrix Pi) {
  
  double one = 1;
  int N = Pi(_, 0).size();
  int K = Pi(0, _).size();
  NumericVector entropies (N, one);
  
  for (int i = 0; i < N ; ++i) {
    
    for (int k = 0; k < K; ++k) {
      
      entropies[i] = entropies[i] * exp(- log2( pow(Pi(i, k), Pi(i, k)) ));
      
    }
    
    entropies[i] = log(entropies[i]);
    
  }
  
  return(entropies);
  
}