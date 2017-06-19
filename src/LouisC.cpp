#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//using namespace Rcpp;

// [[Rcpp::export]]

Rcpp::List LouisC(arma::mat mX,
                  arma::vec np,
                    arma::mat mlogPost,
                    arma::mat itemparmLC,
                    arma::mat parloc, //J x L
                    arma::vec weight,
                    int SEtype){
  //std::cout << "size of mlogPost: " << arma::size(mlogPost) << std::endl;
  // model must be a vector; its elements must be 0, 1 or 2
//  arma::vec np = arma::ones<arma::vec>(3);
//np = 2*np;
//np(3) = 4;
  int N = mX.n_rows;
  //N = 1;
  int J = mX.n_cols; //test length - J
  int L = itemparmLC.n_cols; //2^K
  int NP = arma::accu(np); // total number of parameters
  int NiP = NP;
  parloc--;
  arma::vec postp = arma::ones<arma::vec>(L);
  if (SEtype==3) {
    //std::cout << "size of postp: " << arma::trans(arma::sum(exp(mlogPost)%(weight*arma::ones<arma::mat>(1,L)),0)) << std::endl;
    postp = arma::trans(arma::sum(exp(mlogPost)%(weight*arma::ones<arma::mat>(1,L)),0))/arma::accu(weight);
    NP = NP + L -1;
  }
  //postp.print();
  arma::mat term1 = arma::zeros<arma::mat>(NP,NP);
  arma::mat term2 = term1;
  arma::mat term3 = term1;
  //std::cout << "size of term1: " << arma::size(term1) << std::endl;
  for (int i=0;i<N;++i){
    //std::cout << "i = " << i << "\n" << std::endl;
    arma::vec yi = arma::trans(mX.row(i));
   // std::cout << "size of y1: " << arma::size(yi) << std::endl;
    arma::rowvec post = exp(mlogPost.row(i));
    arma::mat yiL = yi*arma::ones<arma::mat>(1,L); //J x L
    //std::cout << "size of yiL: " << arma::size(yiL) << "\n" << std::endl;
    arma::mat der1 = (yiL-itemparmLC)/(itemparmLC%(1-itemparmLC)); //J x L first deriv
    //der1.print();
    //std::cout << "\nsize of der1: " << arma::size(der1) << std::endl;
    arma::mat der2 = -1*yiL/arma::pow(itemparmLC,2)-(1-yiL)/arma::pow(1-itemparmLC,2);
    arma::mat sumcp = arma::zeros<arma::mat>(NP,NP);
    arma::mat sum2deriv = sumcp;
    arma::vec sumsco = arma::zeros<arma::vec>(NP);
    //std::cout << "L = " << L << "\n" << std::endl;
    for(int l=0;l<L;++l){
      //std::cout << "l = " << l << "\n" << std::endl;
      int np0 = np(0);
      arma::vec score = arma::zeros<arma::vec>(np0);
      arma::vec diaghess = arma::zeros<arma::vec>(np0);
      score(parloc(0, l)) = der1(0,l);
      diaghess(parloc(0,l)) = der2(0,l);
      for (int j=1;j<J;++j){
        //std::cout << "j = " << j << "\n" << std::endl;
        int npj = np(j);
        arma::vec tmpscj = arma::zeros<arma::vec>(npj);
        arma::vec tmphessj = tmpscj;
        tmpscj(parloc(j,l)) = der1(j,l);
        tmphessj(parloc(j,l)) = der2(j,l);
        score=arma::join_vert(score,tmpscj);
        diaghess=arma::join_vert(diaghess,tmphessj);
        //score.print();
      }
      arma::mat mhess = arma::zeros<arma::mat>(NP,NP);
      if(SEtype==3){
        arma::vec scp = arma::zeros<arma::vec>(L-1);
        arma::vec hessp = scp;
        if(l<L-1){
          scp(l) = 1/postp(l);
          hessp(l) = -1/(postp(l)*postp(l));
        }else{
          scp.fill(-1/postp(l));
          hessp.fill(-1/(postp(l)*postp(l)));
        }
        //scp.print();
        score=arma::join_vert(score,scp);
        diaghess=arma::join_vert(diaghess,hessp);
        mhess = arma::diagmat(diaghess);
        if(l==L-1) {
          mhess.submat(NiP,NiP,NP-1,NP-1).fill(-1/(postp(l)*postp(l)));
        }
      }else{
        mhess = arma::diagmat(diaghess);
      }
      sumcp = sumcp + post(l)*score*arma::trans(score);
      sumsco = sumsco + post(l)*score;
      sum2deriv = sum2deriv + post(l)*mhess;

    }
    term1 = term1 + (-1)*sum2deriv*weight(i);
    term2 = term2 + (-1)*sumcp*weight(i);
    term3 = term3 + weight(i)*(sumsco*arma::trans(sumsco));
  }
  arma::mat An = term1 + term2 + term3;
  //std::cout << "\n An = \n" << std::endl;
  //An.print();
  arma::mat invAn = arma::inv(An);
arma::mat robust = invAn*term3*invAn;
  arma::vec SE=sqrt(arma::diagvec(invAn));
  return Rcpp::List::create(Rcpp::Named("invAn") = invAn,
                            Rcpp::Named("term1") = term1,
                            Rcpp::Named("term2") = term2,
                            Rcpp::Named("term3") = term3,
                            Rcpp::Named("robust") = robust,
                            Rcpp::Named("SE") = SE);

}
