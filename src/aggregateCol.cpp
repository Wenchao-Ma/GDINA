#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]

arma::mat aggregateCol(arma::mat & mX, //N x L
                   arma::vec ind){
  ind--; //indictor-1: c++ style
  arma::vec uniq = arma::unique(ind);
  int N = mX.n_rows;
  int Lj = uniq.n_elem;
  arma::mat output = arma::zeros<arma::mat>(N,Lj);
  for (int l=0;l<Lj;++l){
    arma::uvec loc = arma::find(ind==l);
    output.col(l) = arma::sum(mX.cols(loc),1); //N x 1
  }


  return output;

}
