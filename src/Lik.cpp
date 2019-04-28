#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//using namespace Rcpp;
//using namespace arma;

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
                  arma::mat & mloc,
                  arma::vec & weights){


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
    mlogLik = mX0*arma::trunc_log(mP) + (1-mX1)*arma::trunc_log(1-mP); //N x L
  }else{
    mlogLik = mX*arma::trunc_log(mP) + (1-mX)*arma::trunc_log(1-mP); //N x L
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
  double LL = arma::accu(arma::trunc_log(arma::sum(arma::trunc_exp(mlogPost),1))%weights);
  return LL;
}





// [[Rcpp::export]]

Rcpp::List LikNR(const arma::mat & mpar,
                 const arma::mat & mX,
                  arma::mat vlogPrior, //a vector of log prior or a matrix of log prior: col 1 for group 1; col 2 for group 2
                  arma::vec  vgroup,
                  arma::mat  mloc,
                  arma::vec & weights,
                  int simplify = true){
  Rcpp::List ret;
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
    mlogLik = mX0*arma::trunc_log(mP) + (1-mX1)*arma::trunc_log(1-mP); //N x L
  }else{
    mlogLik = mX*arma::trunc_log(mP) + (1-mX)*arma::trunc_log(1-mP); //N x L
  }

  //joint prob (X_i, alpha_c) = log P(X_i|alpha_c) + log p(alpha_c)
  if(no_mg==1){
    mlogPost = mlogLik + arma::ones<arma::mat>(N,1) * arma::trans(vlogPrior); // N x L
  }else{
    mlogPost = mlogLik;
    for (int g=0;g<no_mg;g++){
      locg=arma::find(vgroup==g);
      mlogPost.rows(locg) += arma::ones<arma::mat>(locg.n_elem,1)*arma::trans(vlogPrior.col(g)); // N x L
    }
  }
  arma::mat msdPost = arma::trunc_exp(mlogPost);// P(X_i|alpha_c)p(alpha_c)
  arma::vec rowsum_msdPost = arma::sum(msdPost,1); //Nx1 sum_c P(X_i|alpha_c)p(alpha_c)
  double LL = arma::dot(arma::trunc_log(rowsum_msdPost),weights); //sum_i weight_i log[ sum_c P(X_i|alpha_c)p(alpha_c)]
  ret["LL"]=LL;

  // normalized posterior log [P(alpha_c|X_i)]

  msdPost/=rowsum_msdPost*arma::ones<arma::mat>(1,L); //N x L

  if(simplify==false){
    ret["loglik"]=mlogLik;
    ret["logpost"]=arma::trunc_log(msdPost);
  }
  //-----> weighted msdPost
  msdPost.each_col() %= weights;

  arma::mat updlogPrior = arma::ones<arma::mat>(L,no_mg);
  if(no_mg==1){
    updlogPrior = arma::trunc_log(arma::trans(arma::sum(msdPost,0)/arma::accu(weights)));//L x 1
  }else{
    for (int g=0;g<no_mg;g++){
      locg = arma::find(vgroup==g);
      updlogPrior.col(g) = arma::trunc_log(arma::sum(msdPost.rows(locg),0).t()/arma::accu(weights(locg))); //updated log priors
    }
  }
  ret["logprior"]=updlogPrior;


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

  ret["Ng"]=Ng;
  ret["Rg"]=Rg;


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


// [[Rcpp::export]]
Rcpp::List Lik_DTM(arma::mat & mP, //S0 x L
                   arma::mat & mX, //N x J
                   arma::vec vC, //J x 1
                   arma::vec vlogPrior){

  int N = mX.n_rows;
  int J = vC.n_rows;
  arma::mat mPt = mP.t(); //L x S0
  int S0 = mPt.n_cols;
  int L = mPt.n_rows;

  //mX contains polytomous responses--each item is in one column
  //X contains dichomotous responses--each item lies in Sj columns
  arma::mat X = arma::zeros<arma::mat>(N,S0);
  //This code converts mX to X, that is, N x J --> N x S0
  int start_col;
  vC++;
  //std::cout << "vC: " << vC << std::endl;
  for (int j=0;j<J;++j){
    if (j==0){
      start_col=0;
    }else{
      start_col+=vC[j-1];
    }
    for (int row=0;row<N;++row){
      if(arma::is_finite(mX(row,j))){
        X(row,start_col+mX(row,j))=1;
      }
    }
  } //X:N x S0
  //for each examinee, the likelihood L(Xi) is calculated
  //Xi contain responses for examinee i
  //all L rows are identical
  arma::mat Lik = arma::zeros<arma::mat>(N,L);
  arma::mat Xi = arma::zeros<arma::mat>(L,S0);
  for (int i=0;i<N;++i){
    Xi.each_row() = X.row(i);
    Xi = Xi % mPt;  //element-wise product --> L x S0
    Xi.elem(find(Xi==0)).ones(); //change all 0 to 1
    Lik.row(i) = trans(prod(Xi,1));//prod -> L x 1==>transpose is needed
  }

  //Calculate posterior
  arma::mat mlogPrior = arma::ones<arma::mat>(N,1)*trans(vlogPrior); //N x 1 * 1 x L --> N x L
  arma::mat mPost = exp(log(Lik) + mlogPrior); //unstandarized posterior N x L
  arma::mat msumPost = sum(mPost,1)*arma::ones<arma::mat>(1,L);//N x L
  arma::mat mlogPost = log(mPost) - log(msumPost); //standarized log posterior N x L
  arma::rowvec updlogPrior = log(mean(exp(mlogPost),0)); //updated log priors

  //calculate Log Likelihood
  arma::vec vpost=sum(mPost,1);
  double LL = sum(log(vpost));
  //Log likelihood and log standarized posterior are returned for stability
  //updated log prior is also returned
  return Rcpp::List::create(Rcpp::Named("loglik") = log(Lik), //N x L
                            Rcpp::Named("logpost") = mlogPost,
                            Rcpp::Named("logprior") = updlogPrior,
                            Rcpp::Named("LL") = LL);

}

