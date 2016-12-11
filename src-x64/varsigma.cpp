#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]

arma::mat varsigma(arma::mat mloc, //L-1 x J
             arma::mat mP,
             arma::vec vw){
  mloc--; //indictor-1: c++ style
  int Q = mloc.n_cols;
  int J = mP.n_cols;
  //std::cout << Q << std::endl;
  arma::mat varsig = arma::zeros<arma::mat>(J,Q);
  //std::cout << varsig << std::endl;
  //std::cout<< "begin"<<std::endl;
  for (int j=0;j<J;++j){//each item
    arma::vec wp = mP.col(j)%vw; //w*p
    //std::cout << "wp" << std::endl;
    //wp.print();
    for(int q=0;q<Q;++q){//each possible q-vector
      arma::vec locq = mloc.col(q);
      arma::vec reducedp = arma::zeros<arma::vec>(max(locq)+1);
      arma::vec reducedw = arma::zeros<arma::vec>(max(locq)+1);

      for (int l=0;l<=max(locq);++l){//latent groups

        arma::uvec q1 = arma::find(locq==l);
        reducedw(l) = arma::accu(vw.elem(q1));
        //if (reducedw(l)>0.000001) {
          reducedp(l) = arma::accu(wp.elem(q1))/reducedw(l);
        //}
      }
      //reducedp.print();
      //reducedw.print();
      double pbar = arma::accu(reducedp%reducedw);
      double Sbar = arma::accu(reducedp%reducedp%reducedw);
      //std::cout << "pbar="<< pbar << std::endl;
      //std::cout << "Sbar="<< Sbar << std::endl;
      varsig(j,q) = Sbar - pbar*pbar;
      //varsig.print();
      }
  }

  return varsig;

}
