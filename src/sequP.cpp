#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]

List sequP(arma::mat mloc,
        arma::mat mpar,
        arma::vec vC){
  mloc--;
  int L = mloc.n_cols;
  int S = mloc.n_rows;
  int J = vC.n_elem;
  arma::mat mcP = arma::zeros<arma::mat>(S,L);
  arma::mat muP = arma::zeros<arma::mat>(S,L);
  for (int l=0;l<L;++l){
    for (int s=0;s<S;++s){
      mcP(s,l)=mpar(s,mloc(s,l));
    }
  }

  int first_row;
  int last_row=-1;
  arma::mat sumP;
  arma::mat muP0; //complete Probability matrix including 0
  for (int j=0;j<J;++j){
    if (j==0){
      first_row = 0;
    }else{
      first_row += vC[j-1];
    }

    last_row += vC(j);
    int current_row=last_row;
    while(current_row>=first_row){
      if (current_row==last_row){
        muP.row(current_row)=prod(mcP.rows(first_row,current_row),0);
      }
      else
      {
        muP.row(current_row)=prod(mcP.rows(first_row,current_row),0)%(1-mcP.row(current_row+1));
      }
      --current_row;
    }
    //calculating prob 0
    sumP=arma::join_cols(1-sum(muP.rows(first_row,last_row),0),muP.rows(first_row,last_row));
    muP0=arma::join_cols(muP0,sumP);  //S0 x L including category 0
  }

  //unconditional probability includes category 0
  //conditional probability does not include category 0

  return Rcpp::List::create(Rcpp::Named("uPr") = muP0,
                            Rcpp::Named("cPr") = mcP);

}
