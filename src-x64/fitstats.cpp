#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]

Rcpp::List fitstats(arma::mat mX,
                  arma::mat LCprob,
                  arma::uvec attgroup) {
  attgroup--;
  arma::mat fullLC = LCprob.rows(attgroup);
  int nr = fullLC.n_rows;
  int nc = fullLC.n_cols;
  int NE = mX.n_rows;
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

