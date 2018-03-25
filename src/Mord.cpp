#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]

Rcpp::List Mord(arma::vec item_no, //first col in Qc matrix
                    arma::mat & LCprob, //S x L
                    arma::vec prior) { //prior is a col vector
  //arma::vec item_no = Qc.col(0);
  arma::vec unique_item_no = arma::unique(item_no);
  int nitem = unique_item_no.n_elem;
  arma::rowvec post = prior.t();
  // matrices for final asymptotic covariance matrix
  arma::mat Xi11 = arma::zeros<arma::mat>(nitem,nitem);
  arma::mat Xi21 = arma::zeros<arma::mat>(nitem*(nitem-1)/2,nitem);
  arma::mat Xi22 = arma::zeros<arma::mat>(nitem*(nitem-1)/2,nitem*(nitem-1)/2);

  // univariate or bivariate expected values - E(Yi), E[YiYj]
  arma::vec uni_exp = arma::zeros<arma::vec>(nitem);
  arma::mat bi_exp = arma::zeros<arma::mat>(nitem,nitem);

  // E[Yi]
  for (int i = 0;i < nitem; ++i){ //item i
    arma::mat lci = LCprob.rows(arma::find(item_no==(i+1))); //Si x L
    double mp = 0;
    for (arma::uword si = 0;si < lci.n_rows; si++){ // sum over all categories si
      mp += arma::accu(lci.row(si)%post)*(si+1); //\sum_c(si+1)P(Yi=si+1|alpha_c)p(alpha_c)
    }
    uni_exp(i) = mp;
  }

  for (int j = 0;j < nitem; ++j){ //item j
    arma::mat lcj = LCprob.rows(arma::find(item_no==(j+1))); //Sj x L
    for (int k = 0; k < nitem; ++k){ //item k
      double bmp = 0;
      arma::mat lck = LCprob.rows(arma::find(item_no==(k+1))); //Sk x L
      if(j==k){
        for (arma::uword sk = 0;sk < lck.n_rows; sk++){ //category sk
          bmp += arma::accu(lck.row(sk)%post)*(sk+1)*(sk+1);

        }
      }else{
        for (arma::uword sj = 0;sj < lcj.n_rows; sj++){ // category sj
          for (arma::uword sk = 0;sk < lck.n_rows; sk++){ //category sk
            bmp += arma::accu(lcj.row(sj)%lck.row(sk)%post)*(sj+1)*(sk+1);

          }
        }
      }

    bi_exp(j,k) = bmp;
    bi_exp(k,j) = bmp;
    }
  }

  // calculating Xi11
  for (int j = 0; j < nitem; ++j){
    for (int k = 0; k < nitem; ++k){
      Xi11(j,k) = Xi11(k,j) = bi_exp(j,k) - uni_exp(j)*uni_exp(k);
    }
  }

  // calculating Xi12
  int ind1 = 0;
  for (int i = 0; i < nitem - 1; ++i){
    arma::mat lci = LCprob.rows(arma::find(item_no==(i+1))); //Si x L
    int ir = lci.n_rows;
    for (int j = i + 1; j < nitem; ++j){
      arma::mat lcj = LCprob.rows(arma::find(item_no==(j+1))); //Sj x L
      int jr = lcj.n_rows;
      for (int k = 0; k < nitem; ++k){
        arma::mat lck = LCprob.rows(arma::find(item_no==(k+1))); //Sk x L
        int kr = lck.n_rows;
        double tmp=0;
        if(i==k){
          for(int jj = 0; jj < jr; ++jj){
            for (int kk = 0; kk < kr; ++kk){
              tmp += arma::accu(lcj.row(jj)%lck.row(kk)%post)*arma::as_scalar(jj+1)*arma::as_scalar(kk+1)*arma::as_scalar(kk+1);
            }
          }
        }else if(k==j){
          for(int ii = 0; ii < ir; ++ii){
            for (int kk = 0; kk < kr; ++kk){
              tmp += arma::accu(lci.row(ii)%lck.row(kk)%post)*arma::as_scalar(ii+1)*arma::as_scalar(kk+1)*arma::as_scalar(kk+1);
            }
          }
        }else{
          for(int jj = 0; jj < jr; ++jj){
            for (int kk = 0; kk < kr; ++kk){
              for(int ii = 0; ii < ir; ++ii){
                tmp += arma::accu(lcj.row(jj)%lck.row(kk)%lci.row(ii)%post)*arma::as_scalar(jj+1)*arma::as_scalar(kk+1)*arma::as_scalar(ii+1);
              }
            }
          }
        }
      Xi21(ind1,k) = tmp - bi_exp(i,j)*uni_exp(k);

      }
      ind1++;
    }
  }

  // calculating Xi22---------------------------
  int ind2 = 0;
  for (int i = 0; i < nitem - 1; ++i){
    arma::mat lci = LCprob.rows(arma::find(item_no==(i+1))); //Si x L
    int ir = lci.n_rows;
    for (int j = i + 1; j < nitem; ++j){
      arma::mat lcj = LCprob.rows(arma::find(item_no==(j+1))); //Sj x L
      int jr = lcj.n_rows;
      int ind3 = 0;
      while(ind3<=ind2){
        for (int k = 0; k < nitem - 1; ++k){
          arma::mat lck = LCprob.rows(arma::find(item_no==(k+1))); //Sk x L
          int kr = lck.n_rows;
          for (int l = k + 1; l < nitem; ++l){
            arma::mat lcl = LCprob.rows(arma::find(item_no==(l+1))); //Sm x L
            int lr = lcl.n_rows;
            double tmp=0;

            if((i==k)&(j==l)){
              for(int ii = 0; ii < ir; ++ii){
                for (int jj = 0; jj < jr; ++jj){
                  tmp += arma::accu(lci.row(ii)%lcj.row(jj)%post)*arma::as_scalar(jj+1)*arma::as_scalar(jj+1)*arma::as_scalar(ii+1)*arma::as_scalar(ii+1);
                }
              }
            }else if(i==k){
              for (int jj = 0; jj < jr; ++jj){
                for(int kk = 0; kk < kr; ++kk){
                  for(int ll = 0; ll < lr; ++ll){
                    tmp += arma::accu(lcj.row(jj)%lck.row(kk)%lcl.row(ll)%post)*arma::as_scalar(jj+1)*arma::as_scalar(kk+1)*arma::as_scalar(kk+1)*arma::as_scalar(ll+1);
                  }
                }
              }
            }else if(j==k){
              for (int ii = 0; ii < ir; ++ii){
                for(int kk = 0; kk < kr; ++kk){
                  for(int ll = 0; ll < lr; ++ll){
                    tmp += arma::accu(lci.row(ii)%lck.row(kk)%lcl.row(ll)%post)*arma::as_scalar(ii+1)*arma::as_scalar(kk+1)*arma::as_scalar(kk+1)*arma::as_scalar(ll+1);
                  }
                }
              }
            }else if(i==l){
              for (int jj = 0; jj < jr; ++jj){
                for(int kk = 0; kk < kr; ++kk){
                  for(int ll = 0; ll < lr; ++ll){
                    tmp += arma::accu(lcj.row(jj)%lck.row(kk)%lcl.row(ll)%post)*arma::as_scalar(jj+1)*arma::as_scalar(kk+1)*arma::as_scalar(ll+1)*arma::as_scalar(ll+1);
                  }
                }
              }
            }else if(j==l){
              for (int ii = 0; ii < ir; ++ii){
                for(int kk = 0; kk < kr; ++kk){
                  for(int ll = 0; ll < lr; ++ll){
                    tmp += arma::accu(lci.row(ii)%lck.row(kk)%lcl.row(ll)%post)*arma::as_scalar(ii+1)*arma::as_scalar(kk+1)*arma::as_scalar(ll+1)*arma::as_scalar(ll+1);
                  }
                }
              }
            }else{
              for(int ii = 0; ii < ir; ++ii){
                for (int jj = 0; jj < jr; ++jj){
                  for(int kk = 0; kk < kr; ++kk){
                    for(int ll = 0; ll < lr; ++ll){
                      tmp += arma::accu(lci.row(ii)%lcj.row(jj)%lck.row(kk)%lcl.row(ll)%post)*arma::as_scalar(jj+1)*arma::as_scalar(kk+1)*arma::as_scalar(ii+1)*arma::as_scalar(ll+1);
                    }
                  }
                }
              }
            }

            Xi22(ind3,ind2) = Xi22(ind2,ind3) = tmp - bi_exp(i,j)*bi_exp(k,l);
            ind3++;
          }
        }
      }

      ind2++;
    }
  }


  return Rcpp::List::create(Rcpp::Named("uni") = uni_exp,
                            Rcpp::Named("bi") = bi_exp,
                            Rcpp::Named("Xi11") = Xi11,
                            Rcpp::Named("Xi21") = Xi21,
                            Rcpp::Named("Xi22") = Xi22);
}

