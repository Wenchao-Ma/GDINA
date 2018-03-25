#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]

List NgRg(arma::mat & mlogPost,
                arma::mat & mX,
                arma::mat & mloc,
                arma::vec weights){
    int J = mX.n_cols;
    //int L = mlogPost.n_cols;
    int Ljmax = mloc.max();
    arma::mat mIndmiss = arma::zeros<arma::mat>(arma::size(mX));
    mIndmiss.elem( arma::find_finite(mX) ).ones(); //missing - 0;nonmissing - 1

      arma::mat Ng = arma::zeros<arma::mat>(J,Ljmax);
      arma::mat Rg = arma::zeros<arma::mat>(J,Ljmax);
      arma::mat mPost = exp ( mlogPost ); //N x L
      mPost.each_col() %= weights;
      --mloc;
    for (int j=0;j<J;++j){ //for each item
      int Kjmax = mloc.row(j).max()+1;
      //missing values in item j are removed from posterior and X
      arma::mat mPostmiss=mPost.rows(arma::find(mIndmiss.col(j)==1));
      arma::mat mXmiss=mX.rows(arma::find(mIndmiss.col(j)==1));
      for (int k=0;k<Kjmax;++k){
        arma::mat tmp = mPostmiss.cols(arma::find(mloc.row(j)==k));
        Ng(j,k) = accu(tmp);
        Rg(j,k) = accu((mXmiss.col(j)*arma::ones<arma::rowvec>(tmp.n_cols)) % tmp);
      }
    }
    return Rcpp::List::create(Rcpp::Named("Ng") = Ng,
    Rcpp::Named("Rg") = Rg);

}
