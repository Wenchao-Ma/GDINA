#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//using namespace Rcpp;
//using namespace arma;

void uP_byref(arma::mat & mP, //J x L -> to be modified
              arma::mat & mloc, //J x L
              arma::mat & mpar){
  // calculate mP J x L matrix P(Xij=1|alpha_c)
  const int L = mloc.n_cols;
  const int J = mloc.n_rows;
  arma::vec mPj(L);
  const arma::vec Lj = arma::max(mloc,1);
  for (int j=0;j<J;++j){
    mPj.fill(0);
    for (int l=0;l<Lj(j);++l){
      mPj.elem(arma::find(mloc.row(j)==(l+1)))+=mpar(j,l);
    }
    mP.row(j) = mPj.t();
  }
}


arma::mat uP2(const arma::mat & mloc, //J x L
             const arma::mat & mpar){
  // calculate mP J x L matrix P(Xij=1|alpha_c)
  const int L = mloc.n_cols;
  const int J = mloc.n_rows;
  const arma::vec Lj = arma::max(mloc,1);
  arma::mat mP = arma::zeros<arma::mat>(J, L);
  arma::vec mPj(L);
  for (int j=0;j<J;++j){

    mPj.fill(0);
    for (int l=0;l<Lj(j);++l){
      mPj.elem(arma::find(mloc.row(j)==(l+1)))+=mpar(j,l);
    }
    mP.row(j) = mPj.t();
  }
  return mP;
}

// [[Rcpp::export]]


Rcpp::List fast_GDINA_EM(arma::mat  mloc,
                         arma::mat & mpar,
                         arma::mat & mX,
                         arma::vec vlogPrior, //a vector of log prior or a matrix of log prior: col 1 for group 1; col 2 for group 2
                         arma::vec model_numeric,
                         arma::uvec maxitr,
                         arma::vec lP,
                         arma::vec uP,
                         arma::vec smallNcorrection,
                         arma::vec vbeta,
                         bool prior = false,
                         double crit = 0.0001){

  // calculate mP L x J matrix P(Xij=1|alpha_c)
  int Ljmax = mloc.max();
  //arma::mat mloc_copy = mloc;
  mloc--; //indictor-1: c++ style
  int J = mloc.n_rows;
  int N = mX.n_rows;
  arma::mat mlogLik;
  arma::mat mlogPost;
  arma::mat msdPost;
  arma::mat mPost;
  arma::uvec locg;
  arma::mat mX0; //missing->0
  arma::mat mX1; //missing->1
  arma::mat mXMissing; //missing index: missing ->0; nonmissing ->1
  arma::mat expR;
  arma::mat expN;
  arma::mat Ng;
  arma::mat Rg;
  arma::mat mP;
  arma::vec vlogPrior0;
  arma::uword itr = 0;
  arma::uword maxmaxitr = arma::max(maxitr);
  while(itr < maxmaxitr){
    vlogPrior0 = vlogPrior;
     mP = uP2(mloc+1, mpar);

    //mP.print("mP=");

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
    mlogPost = mlogLik.each_row() + vlogPrior.t(); // N x L

    mPost = exp(mlogPost);// P(X_i|alpha_c)p(alpha_c)
    msdPost = mPost;
    msdPost.each_col()/=arma::sum(mPost,1); //N x L

    //double LL = arma::accu(log(arma::sum(mPost,1)));
    // update prior dist.
    vlogPrior = log(arma::trans(arma::sum(msdPost,0)/N));//L x 1

    Ng = arma::zeros<arma::mat>(J,Ljmax);
    Rg = arma::zeros<arma::mat>(J,Ljmax);


    if(mX.has_nan()){
      //missing values -> 0
      expR = arma::trans(mX0)*msdPost;//JxN * NxL -> JxL
      expN = arma::trans(mXMissing)*msdPost;//JxN * NxL -> JxL
    }else{
      expR = arma::trans(mX)*msdPost;//JxN * NxL -> JxL
      expN = arma::ones<arma::mat>(J,N)*msdPost;//JxN * NxL -> JxL

    }
    //expR.print();
    arma::mat mpar0 = mpar;
    double maxchange_itr = 0;
    for (int j=0;j<J;++j){ //for each item
      if(itr<maxitr(j)){
        double sum_Ng = 0;
        double sum_Rg = 0;

        arma::rowvec expNj = expN.row(j);
        arma::rowvec expRj = expR.row(j);
        int Kjmax = mloc.row(j).max()+1;
        if(prior==false){
          for (int k=0;k<Kjmax;++k){
            Ng(j,k) = arma::accu(expNj.elem(arma::find(mloc.row(j)==k)));
            Rg(j,k) = arma::accu(expRj.elem(arma::find(mloc.row(j)==k)));
            if(model_numeric(j) == 0){ //G-DINA
              mpar(j,k) = (Rg(j,k) + smallNcorrection(0))/(Ng(j,k) + smallNcorrection(1));
              if(mpar(j,k)<lP(j)){
                mpar(j,k) = lP(j);
              }else if(mpar(j,k)>uP(j)){
                mpar(j,k) = uP(j);
              }
            }else if(model_numeric(j) == 1){ //DINA
              if(k < Kjmax-1){
                sum_Ng += Ng(j,k);
                sum_Rg += Rg(j,k);
              }
              if(k == Kjmax-2){
                for (int k = 0; k < Kjmax-1; ++k){
                  mpar(j,k) = (sum_Rg + smallNcorrection(0))/(sum_Ng + smallNcorrection(1));
                  if(mpar(j,k)<lP(j)){
                    mpar(j,k) = lP(j);
                  }else if(mpar(j,k)>uP(j)){
                    mpar(j,k) = uP(j);
                  }
                }
              }else if(k == Kjmax-1){
                mpar(j,k) = (Rg(j,k) + smallNcorrection(0))/(Ng(j,k) + smallNcorrection(1));
                if(mpar(j,k)<lP(j)){
                  mpar(j,k) = lP(j);
                }else if(mpar(j,k)>uP(j)){
                  mpar(j,k) = uP(j);
                }
              }
            }else if(model_numeric(j) == 2){ //DINO
              if(k > 0){
                sum_Ng += Ng(j,k);
                sum_Rg += Rg(j,k);
              }
              if(k == Kjmax-1){
                for (int k = 1; k < Kjmax; ++k){
                  mpar(j,k) = (sum_Rg + smallNcorrection(0))/(sum_Ng + smallNcorrection(1));
                  if(mpar(j,k)<lP(j)){
                    mpar(j,k) = lP(j);
                  }else if(mpar(j,k)>uP(j)){
                    mpar(j,k) = uP(j);
                  }
                }
              }else if(k == 0){
                mpar(j,k) = (Rg(j,k) + smallNcorrection(0))/(Ng(j,k) + smallNcorrection(1));
                if(mpar(j,k)<lP(j)){
                  mpar(j,k) = lP(j);
                }else if(mpar(j,k)>uP(j)){
                  mpar(j,k) = uP(j);
                }
              }
            }

            if(maxchange_itr < std::abs(arma::as_scalar(mpar(j,k))-arma::as_scalar(mpar0(j,k)))){
              maxchange_itr = std::abs(arma::as_scalar(mpar(j,k))-arma::as_scalar(mpar0(j,k)));
            }
            //std::cout << "j = " << j << "change = " << arma::as_scalar(mpar(j,k))-arma::as_scalar(mpar0(j,k)) << std::endl;

          }
        }else{
          for (int k=0;k<Kjmax;++k){
            Ng(j,k) = arma::accu(expNj.elem(arma::find(mloc.row(j)==k)));
            Rg(j,k) = arma::accu(expRj.elem(arma::find(mloc.row(j)==k)));
            if(model_numeric(j) == 0){ //G-DINA
              mpar(j,k) = (Rg(j,k) + vbeta(0) - 1)/(Ng(j,k) + arma::accu(vbeta) - 2);
              if(mpar(j,k)<lP(j)){
                mpar(j,k) = lP(j);
              }else if(mpar(j,k)>uP(j)){
                mpar(j,k) = uP(j);
              }
            }else if(model_numeric(j) == 1){ //DINA
              if(k < Kjmax-1){
                sum_Ng += Ng(j,k);
                sum_Rg += Rg(j,k);
              }
              if(k == Kjmax-2){
                for (int k = 0; k < Kjmax-1; ++k){
                  mpar(j,k) = (sum_Rg + vbeta(0) - 1)/(sum_Ng + arma::accu(vbeta) - 2);
                  if(mpar(j,k)<lP(j)){
                    mpar(j,k) = lP(j);
                  }else if(mpar(j,k)>uP(j)){
                    mpar(j,k) = uP(j);
                  }
                }
              }else if(k == Kjmax-1){
                mpar(j,k) = (Rg(j,k) + vbeta(0) - 1)/(Ng(j,k) + arma::accu(vbeta) - 2);
                if(mpar(j,k)<lP(j)){
                  mpar(j,k) = lP(j);
                }else if(mpar(j,k)>uP(j)){
                  mpar(j,k) = uP(j);
                }
              }
            }else if(model_numeric(j) == 2){ //DINO
              if(k > 0){
                sum_Ng += Ng(j,k);
                sum_Rg += Rg(j,k);
              }
              if(k == Kjmax-1){
                for (int k = 1; k < Kjmax; ++k){
                  mpar(j,k) = (sum_Rg + vbeta(0) - 1)/(sum_Ng + arma::accu(vbeta) - 2);
                  if(mpar(j,k)<lP(j)){
                    mpar(j,k) = lP(j);
                  }else if(mpar(j,k)>uP(j)){
                    mpar(j,k) = uP(j);
                  }
                }
              }else if(k == 0){
                mpar(j,k) = (Rg(j,k) + vbeta(0) - 1)/(Ng(j,k) + arma::accu(vbeta) - 2);
                if(mpar(j,k)<lP(j)){
                  mpar(j,k) = lP(j);
                }else if(mpar(j,k)>uP(j)){
                  mpar(j,k) = uP(j);
                }
              }
            }

            if(maxchange_itr < std::abs(arma::as_scalar(mpar(j,k))-arma::as_scalar(mpar0(j,k)))){
              maxchange_itr = std::abs(arma::as_scalar(mpar(j,k))-arma::as_scalar(mpar0(j,k)));
            }
            //std::cout << "j = " << j << "change = " << arma::as_scalar(mpar(j,k))-arma::as_scalar(mpar0(j,k)) << std::endl;

          }
        }

      }


    }
    itr++;
    maxchange_itr = std::max(maxchange_itr,arma::as_scalar(arma::max(arma::abs(arma::exp(vlogPrior)-arma::exp(vlogPrior0)))));
    std::printf("\nIter =%4d",itr);
    std::printf("  Max. abs. change = %.5f",maxchange_itr);
    std::printf("  Deviance = %.2f",-2 * arma::accu(log(arma::sum(mPost,1))));
    if(maxchange_itr < crit) break;

  }


  Rcpp::List ret;
  ret["ip"]=mpar;
  ret["logprior"]=vlogPrior;
  ret["itr"]=itr;

  return ret;
}
