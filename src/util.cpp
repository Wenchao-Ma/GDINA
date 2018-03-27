// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;




// [[Rcpp::export]]

arma::umat combnCpp(double n, double k) {

  double n_subsets = Rf_choose(n, k);

  arma::umat out = arma::zeros<arma::umat>(k, n_subsets);

  arma::uvec a = arma::linspace<arma::uvec>(1, k, k);

  out.col(0) = a;

  int m = 0;

  int h = k;

  arma::uvec j;



  for(long int i = 1; i < n_subsets; i++){

    if(m < (n - h)){

      h = 1;

      m = a(k - 1);

      j = arma::linspace<arma::uvec>(1, 1, 1);

    }

    else{

      m = a(k - h - 1);

      ++h;

      j = arma::linspace<arma::uvec>(1, h, h);

    }

    a.elem(k - h - 1 + j) = m + j;

    out.col(i) = a;

  }

  return out;

}

// [[Rcpp::export]]

arma::mat rowProd(arma::mat & m,
                  arma::vec & v) {

  arma::rowvec vv =  arma::trans(v);
  arma::mat mm = m.each_row()%vv;
  return mm;
}


// [[Rcpp::export]]

arma::umat alpha2(const int K) {
  arma::umat D;
  arma::umat A;
    A.eye(K,K);

    arma::umat B = arma::ones<arma::umat>(K,std::pow(2.0,K)-K-1);
    int j = 0;
    for(int l=2;l<K;++l){
      arma::umat cob = combnCpp(K, l);
      cob--;
      int tot = Rf_choose(K, l);
      for(int i = 0; i< tot; ++i){
        B.col(j) = arma::sum(A.cols(cob.col(i)),1);
        ++j;
      }
    }
    arma::umat C = arma::join_rows(A,B);
    D = arma::trans(arma::join_rows(arma::zeros<arma::umat>(K,1),C));

  return D;

}

// [[Rcpp::export]]
arma::umat alphap(const arma::uvec maxlevel) {
  IntegerVector maxKj = wrap(maxlevel+1);
  const int K = maxKj.size();
  IntegerVector reptimes = cumprod(rev(maxKj));
  int pd = max(reptimes);
  reptimes.push_front(1);
  IntegerMatrix ret(pd,K);
  int kk = 0;
  for(int k=K-1;k>=0;--k){
    IntegerVector A = seq_len(maxKj(k));
    A = A - 1;
    IntegerVector B = rep_each(A,reptimes(kk));
    IntegerMatrix::Column retcol = ret( _, k);
    retcol = rep_len(B,pd);
    ++kk;
  }
  arma::umat ret2 = as<arma::umat>(ret);
  return ret2;
}

// [[Rcpp::export]]

arma::mat ColNormalize(arma::mat & X){
  arma::mat Y = X;
  arma::rowvec denom = arma::sum(X,0);
  Y.each_row() /= denom;
  return Y;
}

// [[Rcpp::export]]

arma::mat RowNormalize(arma::mat & X){
  arma::mat Y = X;
  arma::vec denom = arma::sum(X,1);
  Y.each_col() /= denom;
  return Y;
}

// [[Rcpp::export]]

double Pr_2PL(double theta, double a, double b) {

  double L = a*theta+b;
  L = 1/(1+std::exp(-1.0*L));
  return L;

}

// [[Rcpp::export]]

arma::mat Pr_2PL_vec(const arma::vec & theta, //N x 1
                     const arma::vec & a,     //K x 1
                     const arma::vec & b,     //K x 1
                     const double minvalue = 1e-16,
                     const double maxvalue= 1 - 1e-16) {
  int N = theta.n_elem;
  arma::mat P = 1/(1+1/arma::exp(theta*a.t()+arma::ones<arma::mat>(N,1)*b.t()));
  P.elem(arma::find(P<minvalue)).fill(minvalue);
  P.elem(arma::find(P>maxvalue)).fill(maxvalue);
  return P; // N x K
}

// [[Rcpp::export]]

arma::mat logLikPattern(arma::mat AlphaPattern, //2^K x K
                       arma::vec theta, //quand point
                       arma::vec a,
                       arma::vec b){
  arma::mat P = arma::trans(Pr_2PL_vec(theta,a,b)); //K x nnodes
  arma::mat logP = AlphaPattern*arma::log(P) + (1-AlphaPattern)*arma::log(1-P);
  return logP; // 2^K x nnodes log P(AlphaPattern|theta_q,a,b)
}

// [[Rcpp::export]]

arma::mat PostTheta(arma::mat AlphaPattern, //2^K x K
                        arma::vec theta, //quand point
                        arma::vec f_theta, // weights
                        arma::vec a,
                        arma::vec b){

  int N = AlphaPattern.n_rows; //2^K
  arma::mat logP = logLikPattern(AlphaPattern, theta, a, b); // 2^K x nnodes
  arma::mat jointP = arma::exp(logP+arma::ones<arma::mat>(N,1)*arma::log(arma::trans(f_theta)));// 2^K x nnodes
  arma::vec denom = arma::sum(jointP,1);
  arma::mat post = jointP.each_col() / denom; // 2^K x nnodes P(theta_q|AlphaPattern)
  return post;
}

// [[Rcpp::export]]
Rcpp::List expectedNR(arma::mat AlphaPattern, //2^K x K
                      arma::vec nc, //2^K x 1
                    arma::vec theta, //quand point
                    arma::vec f_theta, // weights
                    arma::vec a,
                    arma::vec b){
  int Q = f_theta.n_elem;
  int N = AlphaPattern.n_rows; //2^K
  int K = log2(N);
  arma::mat post = PostTheta(AlphaPattern, theta, f_theta, a,b); // 2^K x nnodes P(theta_q|AlphaPattern)
  post.each_col()%=nc;
  arma::vec n = arma::sum(post,0).t();
  arma::mat r = arma::zeros<arma::mat>(K,Q);
  for(int k=0;k<K;++k){
    r.row(k) = arma::sum(post.rows(arma::find(AlphaPattern.col(k)==1)),0);
  }
  Rcpp::List ret;
  ret["n"]=n;
  ret["r"]=r.t();
  return ret;
}

// [[Rcpp::export]]

arma::vec logP_AlphaPattern(arma::mat AlphaPattern, //2^K x K
                       arma::vec theta, //quand point
                       arma::vec f_theta, // weights
                       arma::vec a,
                       arma::vec b){
  int N = AlphaPattern.n_rows;
  arma::mat logP = logLikPattern(AlphaPattern, theta, a,b); // 2^K x nnodes
  arma::vec lP = arma::log(arma::sum(arma::exp(logP + arma::ones<arma::mat>(N,1)*log(arma::trans(f_theta))),1));
  return lP; //log pi_c 2^K x 1
}

// [[Rcpp::export]]

double HoIRTlogLik(arma::mat AlphaPattern, //2^K x K
                   arma::vec ns,
                   arma::vec theta, //quand point
                   arma::vec f_theta, // weights
                   arma::vec a,
                   arma::vec b){
  arma::vec logL = logP_AlphaPattern(AlphaPattern,theta,f_theta,a,b);
  double L = arma::accu(ns%logL); //sum_c^{2^K} n_c log(pi_c)

  return L;
}

// [[Rcpp::export]]

double HoIRTlogLik3(arma::vec & ns,
                   arma::mat & mX,
                   arma::vec & theta, //quand point
                   arma::vec & f_theta, // weights
                   arma::vec a,
                   arma::vec b){
  //int N = ns.n_elem;
  arma::mat P1 = arma::trans(Pr_2PL_vec(theta,a,b)); //K x nnodes
  double L = arma::accu(ns%log(sum(exp(mX*log(P1) + (1-mX)*log(1-P1)+arma::ones<arma::mat>(ns.n_elem,1)*log(arma::trans(f_theta))),1)));

  return L;
}
// [[Rcpp::export]]
double incomplogL(arma::vec a,
                  arma::vec b,
                  arma::mat & logL,//N x 2^K
                  arma::mat & AlphaPattern, //2^K x K
                  arma::vec theta, //quand point
                  arma::vec f_theta // weights
                  ){
  arma::mat lP = exp(logLikPattern(AlphaPattern, theta, a, b)); //2^K x nnodes P(alpha_c|theta_s)
  double L = arma::accu(log(sum(exp(logL)*rowProd(lP,f_theta),1))); //N x 2^K * 2^K x nnodes
  return L;
}

// [[Rcpp::export]]

arma::umat designM(const int Kj,
                  const int model) {
  arma::umat M;
  arma::umat malpha;
  int Lj = std::pow(2.0,Kj);
  if (model==0){//GDINA
      malpha = alpha2(Kj);
    arma::umat M0 = join_rows(arma::ones<arma::umat>(Lj,1),malpha);
    if(Kj>=2){
      arma::umat Mr = arma::ones<arma::umat>(Lj,Lj-M0.n_cols);
      int m = 0;
      for(int j=2;j<=Kj;++j){
        arma::umat comb = combnCpp(Kj, j);
        comb--;
        double tot = Rf_choose(Kj, j);
        for(int i=0;i<tot;++i){
          Mr.col(m) = arma::prod(malpha.cols(comb.col(i)),1);
          m++;
        }

      }
      M = join_rows(M0,Mr);
    }else{
      M = M0;
    }

  }else if(model==1){//DINA
    M = arma::ones<arma::umat>(Lj,2);
    M(arma::span(0,Lj-2),arma::span(1,1)).fill(0);
  }else if(model==2){//DINO
    M = arma::ones<arma::umat>(Lj,2);
    M(0,1)=0;
  }else if ( (model>=3) & (model<=5) ){//ACDM

    malpha = alpha2(Kj);
    M = join_rows(arma::ones<arma::umat>(Lj,1),malpha);
  }
  return M;
}

// [[Rcpp::export]]

arma::uvec matchMatrix(arma::umat A,
                  arma::umat B) {
  //Caveat: rows of A need to be a subset of rows of B
  // this function return a vector of length equal to the number of rows of B
  // if element i = j, row i of B is the same as row j of A
  // if a row of B does not match with any row of A, it has an element of k = nrow(A) + 1
  int Ar = A.n_rows;
  int Br = B.n_rows;
  arma::uvec R(Br);
  R.fill(0);
  arma::umat one = arma::ones<arma::umat>(Br,1);
  for(int r=0;r<Ar;++r){
    arma::uvec bo = all(B == one * A.row(r),1);
    if(any(bo)){
      R.elem(arma::find(bo)).fill(r);
    }

  }
  R++;
  return R;
}

// [[Rcpp::export]]

arma::umat eta(arma::umat & Q) {
  // can be used to replace eta.loc function
  int K = Q.n_cols;
  int J = Q.n_rows;
  int maxQ = Q.max();
  arma::uvec maxcolQ = arma::trans(arma::max(Q,0));
  arma::umat patt;
  int L = arma::prod(maxcolQ+1);
  arma::umat parloc = arma::zeros<arma::umat>(L,J);
  if(maxQ==1){
    patt = alpha2(K);
  }else{
    patt = alphap(maxcolQ);
  }
  for(int j=0;j<J;++j){
    arma::umat pattj;
    arma::urowvec Qj = Q.row(j);
    arma::uvec locj = arma::find(Qj>=1);

    arma::umat transform_patt = patt.cols(locj);
    pattj = alpha2(locj.n_elem);
    if(maxQ>1){
      arma::mat repQj2 = arma::ones<arma::mat>(transform_patt.n_rows,1)*arma::trans(Qj.elem(locj));
      transform_patt.elem(arma::find(transform_patt<repQj2)).zeros();
      transform_patt.elem(arma::find(transform_patt>=repQj2)).ones();
    }
    parloc.col(j) = matchMatrix(pattj, transform_patt);
  }
  return arma::trans(parloc); //J x L
}

