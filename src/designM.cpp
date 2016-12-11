#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]

arma::mat designM(arma::mat malpha,
                  int model) {
  int Lj = malpha.n_rows;
  arma::mat M;
  if (model==0){//GDINA
    M = arma::ones<arma::mat>(1,1);
  }else if(model==1){//DINA
    M = arma::ones<arma::mat>(Lj,2);
    M(arma::span(0,Lj-2),arma::span(1,1)).fill(0);
  }else if(model==2){//DINO
    M = arma::ones<arma::mat>(Lj,2);
    M(0,1)=0;
  }else if (model>=3){//ACDM
    M = join_rows(arma::ones<arma::mat>(Lj,1),malpha);
  }
   return (M);
}

