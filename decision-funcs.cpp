// #include <tr1/random>
#include <Rcpp.h>
#include <vector>
#include <boost/math/special_functions/beta.hpp>

using namespace Rcpp;
using std::vector;

// utility functions
double u_identity(double&x, NumericVector& u){
  return x;
}

double u_pow(double&x, NumericVector& u){
  return pow(x, u[0]);
}

double u_exp(double&x, NumericVector& u){
  double t = u[0];
  if(t > 0)
    return (1 - exp(-t*x))/(1 - exp(-t));
  else
    return (exp(-t*x) - 1)/(exp(-t) - 1);
}

double u_exp_log(double&x, NumericVector& u){
  // log parametrization 
  double t = log(u[0]);
  if(t < 0)
    return (1 - exp(t*x))/(1 - exp(t));
  else
    return (exp(t*x) - 1)/(exp(t) - 1); 
}

// weighting functions

double w_identity(double& p, NumericVector& r){
  return p;
}

double w_prelec(double& p, NumericVector& r){
  if(r.length() > 1){
    return exp(-r[1]*pow(-log(p), r[0]));
  } else {
    return exp(-pow(-log(p), r[0]));
  }
}

double w_ibeta(double& p, NumericVector& r) {
  return boost::math::ibeta(r[0], r[1], p);
}

double w_TK(double& p, NumericVector& r){
  double num;
  num = pow(p, r[0]);
  return num/pow(num + pow(1 + p, r[0]), 1/r[0]);
}

// double w_linlog(double& p, NumericVector& r){
//   double num;
//   if(r.length() > 1){
//     num = r[1]*pow(p, r[0]);
//   } else {
//     num = pow(p, r[0]);
//   }
//   return num/(num + pow(1 + p, r[0]));
// }

double w_pow(double& p, NumericVector& r){

  if(r.length() > 1){
    double H = 1 - pow(1 - p, r[1]);
    double L = pow(p, r[0]);
    return L*p + H*(1 - p);
  } else {
    return pow(p, r[0]);
  }
}

double pt(vector<double>& xx, vector<double>& pp, double& rp,
          double& l, NumericVector& u, NumericVector& r,
          double (*fU)(double& x, NumericVector& u),
          double (*fW)(double& p, NumericVector& r)){
  
  // map over Xs and compute PT(rpx)
  // note: probabilities must be normalized to 1
  // note: NAs in Xs are critical! don't put 0s
  // note: Xs are assumed sorted

  int nR = r.length();

  // negative side (from bottom)
  double accum=0, psum=0, wpH=0, wpL=0;
  for( int j=0; j < xx.size(); j++ ){
    double x=xx[j], p = pp[j];
    // assume a sorted vector
    if( rp <= x ) break;
    if( p > 0 && !NumericVector::is_na(x) ){
      psum += p;
      wpH = fW(psum, r);
      double dx = rp - x;
      accum += -fU(dx, u)*(wpH - wpL)*l;
      wpL = wpH;
    }
  }

  // positive side starting from the top
  psum=0; wpH=0, wpL=0;
  for( int j=(xx.size() - 1); j >= 0; j-- ){
    double x=xx[j], p = pp[j];
    // assume a sorted vector
    if( rp >= x ) break;
    if( p > 0 && !NumericVector::is_na(x) ){
      psum += p;
      wpH = fW(psum, r);
      double dx = x - rp;
      accum += fU(dx, u)*(wpH - wpL);
      wpL = wpH;
    }
  }
  return accum;
}

NumericVector SPT(NumericMatrix& X, NumericMatrix& P,
                  NumericMatrix& rpX, NumericMatrix& rpP,
                  NumericVector& L, NumericMatrix& U, NumericMatrix& R,
                  double (*fU)(double& x, NumericVector& u),
                  double (*fW)(double& p, NumericVector& r)){
  /*
    STOCHASTIC PT:
    X - matrix of outcomes
    P - matrix of probabilities (with 0s)
    rpX - matrix of RPs
    rpP - matrix of coresponding RP probs
    L - vector of loss aversion parameters
    U - vector of power parameters of utility functions
    R - matrix of weighting function parameters

    Notes: p vector need not be normalized. NAs in X and rpX are skipped and
    coresponding Ps ignored
  */
  
  int N = X.nrow(), ncol=X.ncol();
  vector<double> x(ncol), p(ncol);
  if(N != P.nrow() || N != L.length() || N != U.length() || N != R.nrow() ||
     N != rpX.nrow() || N != rpP.nrow())
    Rf_error("nrows of X, P, L, U and R must be the same");
  if(X.ncol() != P.ncol())
    Rf_error("ncol of X and P is not the same");
  if(rpX.ncol() != rpP.ncol())
    Rf_error("ncol of rpX and rpP is not the same");

  NumericVector out(N), r(R.ncol());
  
  for( int row=0; row < N; row++ ){
    double accum = 0, psum=0, l=L[row];
    // map over RPs and compute mixture PT
    // note: probabilities need not be normalized

    NumericVector u=U.row(row), r = R.row(row);

    for( int j=0; j < ncol; j++ ){
      int ixj = row + j*N;
      x[j] = X[ixj];
      p[j] = P[ixj];
    }
    
    for( int i=0; i < rpX.ncol(); i++ ){
      int cixi = row + i*N;
      double rpx = rpX[cixi], rpp=rpP[cixi];
      if( rpp > 0 && !NumericVector::is_na(rpx) ){
        psum += rpp;
        accum += rpp*pt(x, p, rpx, l, u, r, fU, fW);
      }
    }
    if(psum < 0.000001)
      // all RPs are NAs
      out[row] = NA_REAL;
    else 
      out[row] = accum/psum;
  }
  return(out);
}

// [[Rcpp::export]]
NumericVector IBeta(NumericVector& x, double& a, double& b){
  NumericVector out(x.length()), pars = NumericVector::create(a, b);
  for(int i = 0; i < x.length(); i++){
    out[i] = w_ibeta(x[i], pars);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector c_SPT_pow(NumericMatrix& X, NumericMatrix& P,
                        NumericMatrix& rpX, NumericMatrix& rpP,
                        NumericVector& L, NumericMatrix& U){
  NumericMatrix R(U.length()); // this is a waste but what can we do
  return SPT(X, P, rpX, rpP, L, U, R, u_pow, w_identity);
}

// [[Rcpp::export]]
NumericVector c_SPT(NumericMatrix& X, NumericMatrix& P,
                    NumericMatrix& rpX, NumericMatrix& rpP,
                    NumericVector& L, NumericMatrix& U, NumericMatrix& R,
                    vector<std::string> util_type,  vector<std::string> weight_type){
  if(util_type.size() != 1)
    Rf_error("'util_type' parameter must be a character vector of length 1");
  if(weight_type.size() != 1)
    Rf_error("'weight_type' parameter must be a character vector of length 1");

  std::string util = util_type[0];
  std::string weight = weight_type[0];
  if(util == "pow"){
    if(weight == "prelec") return SPT(X, P, rpX, rpP, L, U, R, u_pow, w_prelec);
    if(weight == "pow") return SPT(X, P, rpX, rpP, L, U, R, u_pow, w_pow);
    if(weight == "TK") return SPT(X, P, rpX, rpP, L, U, R, u_pow, w_TK);
    if(weight == "ibeta") return SPT(X, P, rpX, rpP, L, U, R, u_pow, w_ibeta);
    Rf_error("weight_type parameter must be one of prelec, pow, TK, ibeta");
  } 

  if(util == "exp_log"){
    if(weight == "prelec") return SPT(X, P, rpX, rpP, L, U, R, u_exp_log, w_prelec);
    if(weight == "pow") return SPT(X, P, rpX, rpP, L, U, R, u_exp_log, w_pow);
    if(weight == "TK") return SPT(X, P, rpX, rpP, L, U, R, u_exp_log, w_TK);
    if(weight == "ibeta") return SPT(X, P, rpX, rpP, L, U, R, u_exp_log, w_ibeta);
    Rf_error("weight_type parameter must be one of prelec, pow, TK, ibeta");
  }

  Rf_error("util_type parameter must be one of pow, exp_log");
}

//// DIRECT CALLS 

// [[Rcpp::export]]
NumericVector c_SPT_pow_prelec(NumericMatrix& X, NumericMatrix& P,
                               NumericMatrix& rpX, NumericMatrix& rpP,
                               NumericVector& L, NumericMatrix& U, NumericMatrix& R){
  return SPT(X, P, rpX, rpP, L, U, R, u_pow, w_prelec);
}

// // [[Rcpp::export]]
// NumericVector c_SPT_pow_linlog(NumericMatrix& X, NumericMatrix& P,
// 			       NumericMatrix& rpX, NumericMatrix& rpP,
// 			       NumericVector& L, NumericMatrix& U, NumericMatrix& R){
//   return SPT(X, P, rpX, rpP, L, U, R, u_pow, w_linlog);
// }

// [[Rcpp::export]]
NumericVector c_SPT_pow_pow(NumericMatrix& X, NumericMatrix& P,
                            NumericMatrix& rpX, NumericMatrix& rpP,
                            NumericVector& L, NumericMatrix& U, NumericMatrix& R){
  return SPT(X, P, rpX, rpP, L, U, R, u_pow, w_pow);
}

// [[Rcpp::export]]
NumericVector c_SPT_pow_TK(NumericMatrix& X, NumericMatrix& P,
                           NumericMatrix& rpX, NumericMatrix& rpP,
                           NumericVector& L, NumericMatrix& U, NumericMatrix& R){
  return SPT(X, P, rpX, rpP, L, U, R, u_pow, w_TK);
}

// [[Rcpp::export]]
NumericVector c_SPT_pow_ibeta(NumericMatrix& X, NumericMatrix& P,
                              NumericMatrix& rpX, NumericMatrix& rpP,
                              NumericVector& L, NumericMatrix& U, NumericMatrix& R){
  return SPT(X, P, rpX, rpP, L, U, R, u_pow, w_ibeta);
}
