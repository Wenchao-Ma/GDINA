#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//using namespace Rcpp;

// [[Rcpp::export]]

Rcpp::List scorefun(arma::mat mX,
              arma::mat mlogPost,
              arma::mat itmpar,
              arma::mat parloc,
              arma::vec model){
  // model must be a vector; its elements must be 0, 1 or 2
  // SE for A-CDM cannot be caclulated now
  int N = mX.n_rows;
  int J = mX.n_cols;
  int Lj; //the number of latent classes for item j
  --parloc; //latent groups # starting from 0

  //-----Prepare index

  arma::vec c2;
  arma::vec c3;
  arma::rowvec uqloc;

  for (int j=0;j<J;++j){//for each item
    uqloc = arma::unique(parloc.row(j));
    if (model(j)==0){ //G-DINA
      Lj=uqloc.n_elem;
    }else if (model(j)>=3){ //A-CDM
      Lj=uqloc.n_elem;
    }else{//DINA or DINO
      Lj=2;
    }

    //item number: 0011223333...
    arma::vec cc2=arma::ones<arma::vec>(Lj);
    c2=arma::join_vert(c2,j*cc2);
    //kj index: 0101010123...
    arma::vec cc3=arma::linspace<arma::vec>(0,(Lj-1),Lj);
    c3=arma::join_vert(c3,cc3);

  }
  // the number of parameters
  arma::uword npar = c3.n_elem;
  //0,1,2,3,...,npar-1
  arma::vec c1 = arma::linspace<arma::colvec>(0,(c2.n_rows-1),c2.n_rows);
  arma::mat score = arma::zeros<arma::mat>(N,npar);
  arma::rowvec parlocj;

  //-----modify parloc based on different models used
  for (int j=0;j<J;++j){
    parlocj = parloc.row(j);
    if(model(j)==1)
    {//DINA
      parlocj.elem(arma::find(parlocj!=arma::max(parlocj))).zeros(); //all latent groups but 1 -> 0
      parlocj.elem(arma::find(parlocj==arma::max(parlocj))).ones(); // 1 latent group -> 1
      parloc.row(j)=parlocj;
    }
    else if(model(j)==2)
    { //DINO
      parlocj.elem(arma::find(parlocj!=0)).ones(); //all latent groups but 0 -> 1
      parloc.row(j)=parlocj;
    }
  }



  arma::mat mstdPost = exp(mlogPost); //standarized posterior N x L


  for (arma::uword m=0;m<npar;++m){//for each parameter k
    int j = c2(m);//item number for par m
    parlocj = parloc.row(j);

    arma::uvec loc1 = find(parlocj==c3(m));
    arma::vec P1 = itmpar(j,c3(m))*arma::ones<arma::colvec>(N); //N x 1
    score.col(m) = sum(mstdPost.cols(loc1),1) % (mX.col(j)-P1)/(P1%(1-P1)); //N x 1
  }
  arma::mat c = arma::join_horiz(c1,arma::join_horiz(c2,c3));

  return Rcpp::List::create(Rcpp::Named("score") = score,
                            Rcpp::Named("index") = c);

}
