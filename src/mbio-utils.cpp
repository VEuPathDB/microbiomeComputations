#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericMatrix jsd(Rcpp::NumericMatrix df) {
  int ncol = df.ncol();
  Rcpp::NumericMatrix disimilarityMat(ncol, ncol);
  
  for (int i = 0; i < ncol; i++) {
    Rcpp::NumericVector p = df(Rcpp::_, i);

    Rcpp::NumericVector log2p = Rcpp::log(2 * p);
    log2p[!Rcpp::is_finite(log2p)] = 0;

    double pterm = Rcpp::sum(p * log2p);

    for (int j=0; j<= i; j++) {
      Rcpp::NumericVector q = df(Rcpp::_, j);
      
      Rcpp::NumericVector log2q = Rcpp::log(2 * q);
      log2q[!Rcpp::is_finite(log2q)] = 0;

      Rcpp::NumericVector logpq = Rcpp::log(p+q);
      logpq[!Rcpp::is_finite(logpq)] = 0;


      disimilarityMat(i,j) = pterm - Rcpp::sum(p * logpq) + Rcpp::sum(q * log2q) - Rcpp::sum(q * logpq);
      disimilarityMat(j,i) = disimilarityMat(i,j);
    }
  }
  disimilarityMat = disimilarityMat/2;
  return disimilarityMat;
}