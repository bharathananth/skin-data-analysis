#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::cx_vec S(arma::cx_vec x, double lambda){
  return x/abs(x) % (abs(x) - lambda)  % (arma::abs(x) >= lambda);
}

double psi(arma::vec x, double lambda){
  arma::vec y = (x-lambda) % (x >= lambda);
  return arma::norm(y, 1)/arma::norm(y, 2);
}

// [[Rcpp::export]]
Rcpp::List circSVD(arma::cx_mat const &X, double t=1e1, double rtol=1e-5, int maxit=50) { //t is to constrain the number of 'rhy' genes to be determined (decrease) -> based on sum of weights, now it can be very large
  //default: maxit=50, double t=1e6
  int nrow = X.n_rows, ncol = X.n_cols;
  arma::cx_vec u (nrow, arma::fill::randn);
  arma::cx_vec v (ncol, arma::fill::randn);
  v = X.t() * u;
  v = v / norm(v);
  int iter = 0;
  arma::cx_double d, d_last;
  d = arma::as_scalar(u.t() * X * v);
  while (iter++<maxit && std::abs(d/d_last - 1.0)>rtol) {
    u = X * v;
    u = u / arma::abs(u);
    u = u / sqrt(nrow);
    v = X.t() * u;
    if (t < sqrt(ncol)) {
      arma::vec vtilde = sort(arma::abs(v), "descend");
      int l;
      for (l=1; psi(vtilde, vtilde[l])<t && l<ncol; l++);
      double delta = norm((vtilde-vtilde[l]) % (vtilde >= vtilde[l]))/l*(t*sqrt((l - pow(psi(vtilde, vtilde[l-1]), 2))/(l-pow(t,2)))-psi(vtilde, vtilde[l-1]));
      v = S(v, vtilde[l-1] - delta);
    }
    v = v / arma::norm(v);
    d_last = d;
    d = arma::as_scalar(u.t() * X * v);
  }
  
  return Rcpp::List::create(Rcpp::Named("U")=u, Rcpp::Named("V")=v, Rcpp::Named("D")=d);
}