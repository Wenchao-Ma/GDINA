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


// [[Rcpp::export]]
arma::mat Rljs_DTM(arma::mat & mlogPost,//standarized log posterior N x L
                   arma::mat & mX,
                   arma::vec vC){
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
