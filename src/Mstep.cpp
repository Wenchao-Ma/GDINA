// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;

using namespace arma;



// [[Rcpp::export]]
arma::vec Calc_Pj(const arma::vec  par,
                  const arma::mat  designMj,
                  const int & linkfunc,
                  int boundary = 0,
                  const double  eps = 1e-16){
  arma::vec Pj;
  if(linkfunc==1){ //identity
    Pj = designMj*par;
  }else if(linkfunc==2){//logit
    Pj = exp(designMj*par)/(1 + exp(designMj*par));
  }else if(linkfunc==3){//log
    Pj = exp(designMj*par);
  }
  if(boundary==1){
    Pj.elem(find(Pj<eps)).fill(eps);
    Pj.elem(find(Pj>1-eps)).fill(1-eps);
  }
  return Pj;
}

// [[Rcpp::export]]
arma::vec Calc_Dj(arma::vec  par,
                  const arma::mat  designMj,
                  const int & linkfunc,
                  int  boundary = 0,
                  const double  eps = 1e-16){
  arma::vec fPj, Dj;
  if(boundary==1){
    par.elem(find(par<eps)).fill(eps);
    par.elem(find(par>1-eps)).fill(1-eps);
  }
  if(linkfunc==1){ //identity
    fPj = par;
  }else if(linkfunc==2){//logit
    fPj = log(par)-log(1-par);
  }else if(linkfunc==3){//log
    fPj = log(par);
  }
  Dj = arma::inv(designMj.t()*designMj)*arma::trans(designMj)*fPj;
  return Dj;
}

// [[Rcpp::export]]
arma::mat Calc_Pj_jac(arma::vec  par,
                  arma::mat designMj,
                  int & linkfunc,
                  int boundary = 0,
                  double eps = 1e-16){
  arma::mat ret = designMj;
  if(linkfunc>1){
    arma::vec Pj = Calc_Pj(par = par, designMj = designMj, linkfunc = linkfunc, boundary = boundary, eps = eps);
    if(linkfunc==2){//logit
      ret.each_col() %= Pj%(1-Pj);
    }else if(linkfunc==3){//log
      ret.each_col() %= Pj;
    }
  }
  return ret;
}

// [[Rcpp::export]]
double Mstep_obj_fn(arma::vec par,
                     arma::vec & Nj,
                     arma::vec & Rj,
                     arma::mat designMj,
                     arma::vec & uPj,
                     arma::vec & lPj,
                     int linkfunc,
                     Rcpp::Nullable<Rcpp::IntegerMatrix> ConstrMatrix = R_NilValue, //not used
                     double eps = 1e-16,
                     const int ConstrType = 0,
                     const bool greaterthan0 = true) {

  int boundary = 1;
  arma::vec Pj = Calc_Pj(par = par, designMj = designMj, linkfunc = linkfunc, boundary = boundary, eps = eps);

return(-1*arma::accu(Rj%log(Pj)+(Nj-Rj)%log(1-Pj)));

}

// [[Rcpp::export]]
double Mstep_obj_fn_prior(arma::vec par,
                    arma::vec & Nj,
                    arma::vec & Rj,
                    arma::mat designMj,
                    arma::vec & uPj,
                    arma::vec & lPj,
                    int linkfunc,
                    Rcpp::Nullable<Rcpp::IntegerMatrix> ConstrMatrix = R_NilValue, //not used
                    double eps = 1e-16,
                    const int ConstrType = 0,
                    const bool greaterthan0 = true,
                    double m = 0,    //1st normal prior param
                    double sd = 5) {  //2nd normal prior param


  int boundary = 1;
  arma::vec Pj = Calc_Pj(par = par, designMj = designMj, linkfunc = linkfunc, boundary = boundary, eps = eps);

  // adding normal priors
  double sumprior = arma::accu(arma::normpdf(par, m, sd));

  return(-1*arma::accu(Rj%log(Pj)+(Nj-Rj)%log(1-Pj))-sumprior);

}

// [[Rcpp::export]]
double Mstep_obj_fn_max(arma::vec par,
                    arma::vec & Nj,
                    arma::vec & Rj,
                    arma::mat designMj,
                    arma::vec & uPj,
                    arma::vec & lPj,
                    int linkfunc,
                    Rcpp::Nullable<Rcpp::IntegerMatrix> ConstrMatrix = R_NilValue, //not used
                    double eps = 1e-16,
                    const int ConstrType = 0,
                    const bool greaterthan0 = true) {

  int boundary = 1;
  arma::vec Pj = Calc_Pj(par = par, designMj = designMj, linkfunc = linkfunc, boundary = boundary, eps = eps);
  return(arma::accu(Rj%log(Pj)+(Nj-Rj)%log(1-Pj)));

}
// [[Rcpp::export]]
arma::vec Mstep_obj_gr(arma::vec & par,
                    arma::vec & Nj,
                    arma::vec & Rj,
                    arma::mat  designMj,
                    arma::vec & uPj,
                    arma::vec & lPj,
                    int linkfunc,
                    Rcpp::Nullable<Rcpp::IntegerMatrix> ConstrMatrix = R_NilValue, //not used
                    double eps = 1e-16,
                    const int ConstrType = 0,
                    const bool greaterthan0 = true) {
  int boundary = 0;

  arma::vec Pj = Calc_Pj(par = par, designMj = designMj, linkfunc = linkfunc, boundary = boundary, eps = eps);
  if(linkfunc==1){
    designMj.each_col() %= Rj/Pj-(Nj-Rj)/(1-Pj);
  }else if(linkfunc==2){
    designMj.each_col() %= Rj-(Nj%Pj);
  }else if(linkfunc==3){
    designMj.each_col() %= Rj-Pj%(Nj-Rj)/(1-Pj);
  }
  return(-1*arma::trans(sum(designMj,0)));
}

// [[Rcpp::export]]
arma::vec Mstep_ineq_fn(arma::vec par,
                        const arma::vec Nj,
                        const arma::vec Rj,
                        arma::mat designMj,
                        const double uPj,
                        const double lPj,
                        int linkfunc,
                        Rcpp::Nullable<Rcpp::IntegerMatrix> ConstrMatrix = R_NilValue,
                        double eps = 1e-16,
                        const int ConstrType = 0,
                        const bool greaterthan0 = true) {
  arma::mat constr2;
  if (ConstrMatrix.isNotNull()) {
    constr2 = as<arma::mat>(ConstrMatrix);
  }

  int boundary = 0;
  arma::vec Pj, ret;
  Pj = Calc_Pj(par = par, designMj = designMj, linkfunc = linkfunc, boundary = boundary, eps = eps);
    arma::mat A = arma::join_cols(Pj-lPj,uPj-Pj);
    if(ConstrType==1){ //constr for P and 1-P
      ret = A;
    }else if(ConstrType==2){ //constr for mono constraint
      ret = constr2*Pj;
    }else if(ConstrType==3){//constr for P, 1-P and mono
      ret = arma::join_cols(A,constr2*Pj);
    }
    if(!greaterthan0) ret = -1*ret;

return ret;
}

// [[Rcpp::export]]
NumericMatrix Mstep_ineq_jac(arma::vec par,
                            const arma::vec Nj,
                            const arma::vec Rj,
                            arma::mat designMj,
                            const double uPj,
                            const double lPj,
                            int linkfunc,
                            Rcpp::Nullable<Rcpp::IntegerMatrix> ConstrMatrix = R_NilValue,
                            double eps = 1e-16,
                            const int ConstrType = 0,
                            const bool greaterthan0 = true) {
  arma::mat constr2;
  if (ConstrMatrix.isNotNull()) {
    constr2 = as<arma::mat>(ConstrMatrix);
  }
  Rcpp::NumericMatrix ret;
  if(ConstrType==0){ //no constraint
    ret = R_NilValue;
  }else{
    int boundary = 0;
    arma::mat Pj_jac = Calc_Pj_jac(par = par, designMj = designMj, linkfunc = linkfunc, boundary = boundary, eps = eps);
    arma::mat A = arma::join_cols(Pj_jac,-1*Pj_jac);
    if(ConstrType==1){ //constr for P and 1-P
      ret = wrap(A);
    }else if(ConstrType==2){ //constr for mono constraint
      ret = wrap(constr2*Pj_jac);
    }else if(ConstrType==3){//constr for P, 1-P and mono
      ret = wrap(arma::join_cols(A,constr2*Pj_jac));
    }
    if(!greaterthan0) ret = -1*ret;
  }
  return ret;
}
