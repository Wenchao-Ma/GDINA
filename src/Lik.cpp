#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//using namespace Rcpp;

// [[Rcpp::export]]

Rcpp::List Lik(arma::mat mP,
          arma::mat mX,
          arma::mat vlogPrior, //a vector of log prior or a matrix of log prior: col 1 for group 1; col 2 for group 2
          arma::vec vgroup){ //vgroup vector: 1 for group1; 2 for group2,...

  int N = mX.n_rows;
  int J = mX.n_cols;
  mP=mP.t(); //L x J
  int L = mP.n_rows;
  int no_mg = vgroup.max();
  vgroup--;
  arma::mat mIndmiss = arma::zeros<arma::mat>(arma::size(mX));
  mIndmiss.elem( arma::find_finite(mX) ).ones(); //missing - 0;nonmissing - 1

  //std::cout << "vlogPrior: " << vlogPrior << std::endl;
  //for each examinee, the likelihood L(Xi) is calculated
  //Xi contain responses for examinee i
  //all L rows are identical
  arma::mat logLik = arma::zeros<arma::mat>(L,N);
  arma::mat mlogPrior = arma::zeros<arma::mat>(L,N);
  arma::mat Xi = arma::zeros<arma::mat>(L,J); //response vector for person i
  for (int i=0;i<N;++i){
    //std::cout << "i= " << i << std::endl;
    Xi.each_row() = mX.row(i); //L x J
    Xi = Xi % mP + (1-Xi) % (1-mP);  //element-wise product --> L x J
    arma::mat tmp = Xi.cols(arma::find(mIndmiss.row(i)==1));
    logLik.col(i) = arma::sum(log(tmp),1);//prod -> L x 1
    arma::uword ind = vgroup(i);
    //std::cout << "ind: " << ind << std::endl;
    mlogPrior.col(i) = vlogPrior.col(ind);
  }

  //Calculate posterior
  //arma::mat mlogPrior = arma::ones<arma::mat>(N,1)*arma::trans(vlogPrior); //N x 1 * 1 x L --> N x L
  arma::mat mPost = exp(logLik + mlogPrior); //unstandarized posterior L x N
  //std::cout << "size of mPost: " << arma::size(mPost) << std::endl;
  arma::mat msumPost = arma::ones<arma::mat>(L,1)*arma::sum(mPost,0);//L x N
  //std::cout << "size of msumPost: " << arma::size(msumPost) << std::endl;
  arma::mat mlogPost = log(mPost) - log(msumPost); //standarized log posterior L x N
  //std::cout << "size of mlogPost: " << arma::size(mlogPost) << std::endl;
  //arma::vec updlogPrior = log(arma::mean(exp(mlogPost),1)); //updated log priors - for all groups combination L x 1 (dim=1; each row)
  //std::cout << "size of updlogPrior: " << arma::size(updlogPrior) << std::endl;
  // log priors for each group
  arma::mat updlogPrior = arma::ones<arma::mat>(L,no_mg);
  //std::cout << "size of updlogPrior_g: " << arma::size(updlogPrior_g) << std::endl;

  for (int g=0;g<no_mg;g++){
    updlogPrior.col(g) = log(arma::mean(exp(mlogPost.cols(arma::find(vgroup==g))),1)); //updated log priors
  }
  //logLik = arma::trans(logLik);
  //mlogPrior = arma::trans(mlogPrior);
  //calculate Log Likelihood
  //arma::vec vpost=arma::sum(mPost,0);
  double LL = arma::accu(log(arma::sum(mPost,0)));
  //std::cout << "LL: " << LL << std::endl;
  //Log likelihood and log standarized posterior are returned for stability
  //updated log prior is also returned
  return Rcpp::List::create(Rcpp::Named("loglik") = arma::trans(logLik), //N x L
  Rcpp::Named("logpost") = arma::trans(mlogPost),
  Rcpp::Named("logprior") = updlogPrior,
  //Rcpp::Named("logprior.g") = updlogPrior_g,
  Rcpp::Named("LL") = LL);
}

