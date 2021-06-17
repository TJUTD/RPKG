#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' @noRd
// [[Rcpp::export]]
Rcpp::List fLM(Rcpp::NumericMatrix Xr, Rcpp::NumericVector yr) {
  int n = Xr.nrow(), k = Xr.ncol();
  arma::mat X(Xr.begin(), n, k, false);
  arma::colvec y(yr.begin(), yr.size(), false);
  int df = n - k;

  // fit model y ~ X, extract residuals
  arma::colvec coef = arma::solve(X, y);
  arma::colvec fit = X*coef;
  arma::colvec res = y - X*coef;

  return Rcpp::List::create(Rcpp::Named("coefficients") = coef,
                            Rcpp::Named("fitted.values") = fit,
                            Rcpp::Named("residuals") = res);
}
