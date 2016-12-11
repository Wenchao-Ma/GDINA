#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]

arma::mat uP(arma::mat mloc,
        arma::mat mpar){
mloc--; //indictor-1: c++ style
int L = mloc.n_cols;
int J = mloc.n_rows;
arma::mat mP = arma::zeros<arma::mat>(J,L);
for (int l=0;l<L;++l){
  for (int j=0;j<J;++j){
    mP(j,l)=mpar(j,mloc(j,l));//probability of success of latent class l on item j
  }
}

return mP;
  
}