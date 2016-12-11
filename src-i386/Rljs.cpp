#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]

arma::mat Rljs(arma::mat mlogPost,
               arma::mat mX,
               arma::vec vC){
  //arma::mat mlogPost = Rcpp::as<arma::mat>(logpost);; //standarized log posterior N x L
  //arma::mat mX = Rcpp::as<arma::mat>(Y);
  //arma::vec vC = Rcpp::as<arma::vec>(C);
  vC++;
  int J = mX.n_cols;
  int S0 = sum(vC); //including 0
  int L = mlogPost.n_cols;
  arma::mat mPost = exp ( mlogPost ); //N x L
  arma::mat Rls = arma::zeros<arma::mat>(L,S0);
  int start_col;
  for (int j=0;j<J;++j){
    //find the starting column for item j
    if (j==0){
      start_col=0;
    }else{
      start_col+=vC[j-1];
    }
    int s=0;
    while (s<vC(j)){
      arma::uvec Xjs = arma::find(mX.col(j)==s); //the row number of examinees getting s score on item j
      if (!(Xjs.is_empty())){ //if some examinees get s score
        arma::mat postjs = mPost.rows(Xjs);
        Rls.col(start_col + s)=trans(sum(postjs,0));
      }
      ++s;
    }
  }
  Rls=trans(Rls);
  return (Rls); //Return L x S0 expected number of examinees getting sj score on item j

}
