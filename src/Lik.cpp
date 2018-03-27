#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//using namespace Rcpp;

// [[Rcpp::export]]

arma::mat uP(const arma::mat & mloc, //J x L
             const arma::mat & mpar){
  // calculate mP J x L matrix P(Xij=1|alpha_c)
  const int L = mloc.n_cols;
  const int J = mloc.n_rows;
  const arma::vec Lj = arma::max(mloc,1);
  arma::mat mP = arma::zeros<arma::mat>(J, L);
  for (int j=0;j<J;++j){
    arma::vec mPj(L);
    mPj.fill(0);
    for (int l=0;l<Lj(j);++l){
      mPj.elem(arma::find(mloc.row(j)==(l+1)))+=mpar(j,l);
    }
    mP.row(j) = mPj.t();
  }
  return mP;
}


// [[Rcpp::export]]

double ObsLogLik(const arma::mat & mpar,
                  const arma::mat & mX,
                  arma::mat vlogPrior,
                  arma::vec  vgroup,
                  arma::mat  mloc,
                  arma::vec weights){


  arma::mat mP = uP(mloc, mpar);

  const int no_mg = vgroup.max();
  vgroup--;
  arma::mat mlogLik;
  arma::mat mlogPost;
  arma::uvec locg;
  if(mX.has_nan()){

    arma::mat mX0 = arma::ones<arma::mat>(arma::size(mX)); //missing -> 0
    arma::mat mX1 = arma::ones<arma::mat>(arma::size(mX)); //missing -> 1
    mX0.elem( arma::find_nonfinite(mX) ).zeros(); //missing - 0
    mX1.elem( arma::find_nonfinite(mX) ).ones(); //missing - 1
    mlogLik = mX0*log(mP) + (1-mX1)*log(1-mP); //N x L
  }else{
    mlogLik = mX*log(mP) + (1-mX)*log(1-mP); //N x L
  }

  //joint prob (X_i, alpha_c) = log P(X_i|alpha_c) + log p(alpha_c)
  if(no_mg==1){
    mlogPost = mlogLik.each_row() + vlogPrior.t(); // N x L
  }else{
    mlogPost = mlogLik;
    for (int g=0;g<no_mg;g++){
      locg=arma::find(vgroup==g);
      mlogPost.rows(locg) += arma::ones<arma::mat>(locg.n_elem,1)*arma::trans(vlogPrior.col(g)); // N x L
    }
  }

  //sum_i weight_i log[ sum_c P(X_i|alpha_c)p(alpha_c)]
  double LL = arma::accu(log(arma::sum(exp(mlogPost),1))%weights);
  return LL;
}



// [[Rcpp::export]]

Rcpp::List LikNR(const arma::mat & mpar,
                 const arma::mat & mX,
                  arma::mat vlogPrior, //a vector of log prior or a matrix of log prior: col 1 for group 1; col 2 for group 2
                  arma::vec  vgroup,
                  arma::mat  mloc,
                  arma::vec weights,
                  int simplify = 1){

  // calculate mP L x J matrix P(Xij=1|alpha_c)

  arma::mat mP = uP(mloc, mpar);

  const int Ljmax = mloc.max();
  mloc--; //indictor-1: c++ style
  const int L = mloc.n_cols;
  const int J = mloc.n_rows;
  const int N = mX.n_rows;
  const int no_mg = vgroup.max();
  vgroup--;
  arma::mat mlogLik;
  arma::mat mlogPost;
  arma::uvec locg;
  arma::mat mX0; //missing->0
  arma::mat mX1; //missing->1
  arma::mat mXMissing; //missing index: missing ->0; nonmissing ->1
  if(mX.has_nan()){
    mX0 = mX;
    mX1 = mX;
    mXMissing = mX;
    arma::uvec missingloc = arma::find_nonfinite(mX);
    mX0.elem( missingloc ).zeros(); //missing - 0
    mX1.elem( missingloc ).ones(); //missing - 1
    mXMissing.elem(missingloc).zeros();
    mXMissing.elem(arma::find_finite(mX)).ones();
    mlogLik = mX0*log(mP) + (1-mX1)*log(1-mP); //N x L
  }else{
    mlogLik = mX*log(mP) + (1-mX)*log(1-mP); //N x L
  }

  //joint prob (X_i, alpha_c) = log P(X_i|alpha_c) + log p(alpha_c)
  if(no_mg==1){
    mlogPost = mlogLik.each_row() + vlogPrior.t(); // N x L
  }else{
    mlogPost = mlogLik;
    for (int g=0;g<no_mg;g++){
      locg=arma::find(vgroup==g);
      mlogPost.rows(locg) += arma::ones<arma::mat>(locg.n_elem,1)*arma::trans(vlogPrior.col(g)); // N x L
    }
  }
  arma::mat mPost = exp(mlogPost);// P(X_i|alpha_c)p(alpha_c)
  double LL = arma::accu(log(arma::sum(mPost,1))%weights); //sum_i weight_i log[ sum_c P(X_i|alpha_c)p(alpha_c)]
  // normalized posterior log [P(alpha_c|X_i)]
  arma::mat msdPost = mPost;
  msdPost.each_col()/=arma::sum(mPost,1); //N x L
  msdPost.each_col() %= weights;

  arma::mat updlogPrior = arma::ones<arma::mat>(L,no_mg);
  if(no_mg==1){
    updlogPrior = log(arma::trans(arma::sum(msdPost,0)/arma::accu(weights)));//L x 1
  }else{
    for (int g=0;g<no_mg;g++){
      locg = arma::find(vgroup==g);
      updlogPrior.col(g) = log(arma::sum(msdPost.rows(locg),0).t()/arma::accu(weights(locg))); //updated log priors
    }
  }


  arma::mat Ng = arma::zeros<arma::mat>(J,Ljmax);
  arma::mat Rg = arma::zeros<arma::mat>(J,Ljmax);
  arma::mat expR;
  arma::mat expN;

  if(mX.has_nan()){
    //missing values -> 0
    expR = arma::trans(mX0)*msdPost;//JxN * NxL -> JxL
    expN = arma::trans(mXMissing)*msdPost;//JxN * NxL -> JxL
  }else{
    expR = arma::trans(mX)*msdPost;//JxN * NxL -> JxL
    expN = arma::ones<arma::mat>(J,N)*msdPost;//JxN * NxL -> JxL

  }
  for (int j=0;j<J;++j){ //for each item
    arma::rowvec expNj = expN.row(j);
    arma::rowvec expRj = expR.row(j);
      int Kjmax = mloc.row(j).max()+1;

        for (int k=0;k<Kjmax;++k){
          Ng(j,k) = arma::accu(expNj.elem(arma::find(mloc.row(j)==k)));
          Rg(j,k) = arma::accu(expRj.elem(arma::find(mloc.row(j)==k)));
        }

    }


  Rcpp::List ret;
  ret["LL"]=LL;
  ret["logprior"]=updlogPrior;
  ret["Ng"]=Ng;
  ret["Rg"]=Rg;
  if(simplify==0){
    ret["loglik"]=mlogLik;
    ret["logpost"]=log(msdPost);
  }

  return ret;
}


// [[Rcpp::export]]

Rcpp::List LikNR_LC(const arma::mat & mP,//J x L
                 const arma::mat & mX,
                 arma::mat vlogPrior, //a vector of log prior or a matrix of log prior: col 1 for group 1; col 2 for group 2
                 arma::vec  vgroup,
                 arma::vec weights,
                 int simplify = 1){

  const int L = mP.n_cols;
  const int J = mP.n_rows;
  const int N = mX.n_rows;
  const int no_mg = vgroup.max();
  vgroup--;
  arma::mat mlogLik;
  arma::mat mlogPost;
  arma::uvec locg;
  arma::mat mX0; //missing->0
  arma::mat mX1; //missing->1
  arma::mat mXMissing; //missing index: missing ->0; nonmissing ->1
  if(mX.has_nan()){
    mX0 = mX;
    mX1 = mX;
    mXMissing = mX;
    arma::uvec missingloc = arma::find_nonfinite(mX);
    mX0.elem( missingloc ).zeros(); //missing - 0
    mX1.elem( missingloc ).ones(); //missing - 1
    mXMissing.elem(missingloc).zeros();
    mXMissing.elem(arma::find_finite(mX)).ones();
    mlogLik = mX0*log(mP) + (1-mX1)*log(1-mP); //N x L
  }else{
    mlogLik = mX*log(mP) + (1-mX)*log(1-mP); //N x L
  }

  //joint prob (X_i, alpha_c) = log P(X_i|alpha_c) + log p(alpha_c)
  if(no_mg==1){
    mlogPost = mlogLik.each_row() + vlogPrior.t(); // N x L
  }else{
    mlogPost = mlogLik;
    for (int g=0;g<no_mg;g++){
      locg=arma::find(vgroup==g);
      mlogPost.rows(locg) += arma::ones<arma::mat>(locg.n_elem,1)*arma::trans(vlogPrior.col(g)); // N x L
    }
  }
  arma::mat mPost = exp(mlogPost);// P(X_i|alpha_c)p(alpha_c)
  double LL = arma::accu(log(arma::sum(mPost,1))%weights); //sum_i weight_i log[ sum_c P(X_i|alpha_c)p(alpha_c)]
  // normalized posterior log [P(alpha_c|X_i)]
  arma::mat msdPost = mPost;
  msdPost.each_col()/=arma::sum(mPost,1); //N x L
  msdPost.each_col() %= weights;

  arma::mat updlogPrior = arma::ones<arma::mat>(L,no_mg);
  if(no_mg==1){
    updlogPrior = log(arma::trans(arma::sum(msdPost,0)/arma::accu(weights)));//L x 1
  }else{
    for (int g=0;g<no_mg;g++){
      locg = arma::find(vgroup==g);
      updlogPrior.col(g) = log(arma::sum(msdPost.rows(locg),0).t()/arma::accu(weights(locg))); //updated log priors
    }
  }

  arma::mat expR;
  arma::mat expN;

  if(mX.has_nan()){
    //missing values -> 0
    expR = arma::trans(mX0)*msdPost;//JxN * NxL -> JxL
    expN = arma::trans(mXMissing)*msdPost;//JxN * NxL -> JxL
  }else{
    expR = arma::trans(mX)*msdPost;//JxN * NxL -> JxL
    expN = arma::ones<arma::mat>(J,N)*msdPost;//JxN * NxL -> JxL

  }


  Rcpp::List ret;
  ret["LL"]=LL;
  ret["logprior"]=updlogPrior;
  ret["N"]=expN;
  ret["R"]=expR;
  if(simplify==0){
    ret["loglik"]=mlogLik;
    ret["logpost"]=log(msdPost);
  }

  return ret;
}
