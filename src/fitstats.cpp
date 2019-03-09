#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List fitstats(arma::mat & mX,
                     arma::mat & Xfit,
                     bool cor = true) {
  int nr = Xfit.n_rows;//N*
  int nc = Xfit.n_cols;//J
  int NE = mX.n_rows;//N
  int rep = nr/NE;
  Rcpp::List ret;
  ret["pfit"] = arma::mean(Xfit,0);
  if(cor){
    //correlation
    ret["r"]=arma::cor(mX);
    ret["rfit"] = arma::cor(Xfit);
  }

  //l
  arma::mat l(nc,nc, arma::fill::zeros);
  arma::mat lfit(nc,nc, arma::fill::zeros);
  arma::mat sefit(nc,nc, arma::fill::zeros);
  double n11,n01,n10,n00,N11,N01,N10,N00;
  arma::vec x,X,y,Y;
  bool na,na2;
  for(int c1=0;c1<nc-1;c1++){
    x = mX.col(c1);
    X = Xfit.col(c1);
    na = x.has_nan();
    for(int c2=c1+1;c2<nc;c2++){
      y = mX.col(c2);
      na2 = y.has_nan();
      Y = Xfit.col(c2);
      if(na||na2){
        arma::uvec nonna = arma::intersect(arma::find_finite(x),arma::find_finite(y));
        n11 = dot(x.elem(nonna),y.elem(nonna));
        n00 = dot(1-x.elem(nonna),1-y.elem(nonna));
        n10 = dot(x.elem(nonna),1-y.elem(nonna));
        n01 = dot(1-x.elem(nonna),y.elem(nonna));
      }else{
        n11 = dot(x,y);
        n00 = dot(1-x,1-y);
        n10 = dot(x,1-y);
        n01 = dot(1-x,y);
      }

      N11 = dot(X,Y);
      N00 = dot(1-X,1-Y);
      N10 = dot(X,1-Y);
      N01 = dot(1-X,Y);
      if(n11==0) n11 = 0.5;
      if(n01==0) n01 = 0.5;
      if(n10==0) n10 = 0.5;
      if(n00==0) n00 = 0.5;
      if(N11==0) N11 = 0.5;
      if(N01==0) N01 = 0.5;
      if(N10==0) N10 = 0.5;
      if(N00==0) N00 = 0.5;
      l(c1,c2) = n11*n00/(n01*n10);
      l(c2,c1) = l(c1,c2);
      lfit(c1,c2) = N11*N00/(N01*N10);
      lfit(c2,c1) = lfit(c1,c2);
      sefit(c1,c2) = (1/N11 + 1/N00 + 1/N10 + 1/N01)*rep;
      sefit(c2,c1) = sefit(c1,c2);
    }
  }

  ret["l"] = l;
  ret["lfit"] = lfit;
  ret["sefit"] = sefit;
  return ret;
}




Rcpp::List fitstats2(arma::mat & mX,
                  arma::mat & LCprob,
                  arma::uvec attgroup) {
  attgroup--;
  arma::mat fullLC = LCprob.rows(attgroup);//N* x J
  int nr = fullLC.n_rows;//N*
  int nc = fullLC.n_cols;//J
  int NE = mX.n_rows;//N
  int rep = attgroup.n_elem/NE;

  arma::mat Xfit(nr,nc, arma::fill::zeros);
  arma::mat unif = arma::randu<arma::mat>(nr,nc);
  for (int nc1=0;nc1<nc;nc1++){
    for (int nr1=0;nr1<nr;nr1++){
      if(fullLC(nr1,nc1)>unif(nr1,nc1)){
        Xfit(nr1,nc1) = 1;
      }
    }
  }


  //correlation
  arma::mat r = arma::cor(mX);
  arma::mat rfit = arma::cor(Xfit);
  //l
  arma::mat l(nc,nc, arma::fill::zeros);
  for(int c1=0;c1<nc;c1++){
    arma::vec x = mX.col(c1);
    for(int c2=0;c2<nc;c2++){
      arma::vec y = mX.col(c2);
      l(c1,c2) = arma::accu(x % y) * arma::accu((1-x) % (1-y))/(arma::accu((x) % (1-y)) * arma::accu((1-x) % (y)));
      l(c2,c1) = l(c1,c2);
    }
  }

  arma::mat lfit(nc,nc, arma::fill::zeros);
  for(int c1=0;c1<nc;c1++){
    arma::vec x = Xfit.col(c1);
    for(int c2=0;c2<nc;c2++){
      arma::vec y = Xfit.col(c2);
      lfit(c1,c2) = arma::accu(x % y) * arma::accu((1-x) % (1-y))/(arma::accu((x) % (1-y)) * arma::accu((1-x) % (y)));
      lfit(c2,c1) = lfit(c1,c2);
    }
  }

  // se for l
  arma::mat sefit(nc,nc, arma::fill::zeros);
  for(int c1=0;c1<nc;c1++){
    arma::vec x = Xfit.col(c1);
    for(int c2=0;c2<nc;c2++){
      arma::vec y = Xfit.col(c2);
      sefit(c1,c2) = (1/arma::accu(x % y) + 1/arma::accu((1-x) % (1-y))
        + 1/arma::accu((x) % (1-y)) + 1/arma::accu((1-x) % (y)))*rep;
      sefit(c2,c1) = sefit(c1,c2);
    }
  }

  // column mean for Xfit
  arma::rowvec p = arma::mean(Xfit,0);

  return Rcpp::List::create(Rcpp::Named("r") = r,
                            Rcpp::Named("rfit") = rfit,
                            Rcpp::Named("l") = l,
                            Rcpp::Named("lfit") = lfit,
                            Rcpp::Named("sefit") = sefit,
                            Rcpp::Named("pfit") = p);
}

