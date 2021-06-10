#' Title
#'
#' @param y
#' @param x
#' @param a.order
#' @param sig.level
#' @param crit.type
#' @param nonp.method
#' @param num.resamp
#' @param num.discard
#' @param CPP
#'
#' @return
#' @export
#'
#' @examples
rlmtest <- function(y,
                    x = NA,
                    a.order = 1,
                    sig.level = 0.05,
                    crit.type = c("asymptotic", "nonparametric"),
                    nonp.method = c("bootstrap","permutation"),
                    num.resamp = 1000,
                    num.discard = 0,
                    CPP = T)
{
  # the length of time series y
  n <- length(y)

  # regression
  if (is.na(x[1])) {
    if (CPP) {
      ls.fit <- fLm(matrix(rep(1,n),ncol=1),y)
    } else {
      ls.fit <- lm(y ~ 1)
    }
  } else {
    if (CPP) {
      ls.fit <- fLm(cbind(1,x),y)
    } else {
      ls.fit <- lm(y ~ x)
    }
  }

  # extract residuals
  res <- ls.fit$residuals

  # demean residuals from the regression above ?
  res.demean <- res - mean(res)

  # square the residuals
  res.sq <- res^2
  res.sq.demean <- res.sq-mean(res.sq)

  # run autoregressive ls.fit of a.order on res.sq
  if (CPP) {
    ar.fit <- fLm(cbind(1,res.sq.demean[1:(n-a.order)]), res.sq.demean[(a.order+1):n])
  } else {
    ar.fit <- lm(res.sq.demean[(a.order+1):n]~res.sq.demean[1:(n-a.order)])
  }

  # calculate residuals based on AR(a.order)
  ee <- ar.fit$residuals

  # calculate tau_hat (see the CJS paper)
  SSR0.hat <- sum((res.sq.demean[(a.order+1):n])^2)
  # SSR1.hat <- sum(ee[(a.order+1):n]^2)
  SSR1.hat <- sum(ee^2)
  tau.hat <- (n-a.order)*(1-SSR1.hat/SSR0.hat) # this the observed statistic

  if (crit.type == "asymptotic") {
    pValue <- pchisq(tau.hat, df = a.order, lower.tail = FALSE)
    crit <- qchisq(sig.level, df = a.order, lower.tail = FALSE)
    return(list(stat = tau.hat, crit.val = crit, p.value = pValue))
  } else {
    # resampling test
    tau.bs <- numeric(num.resamp)

    for (j in 1:num.resamp) {

      ### resample residuals
      if(nonp.method == "permutation") {
        # permutation
        e.bs <- sample(res.demean,n,replace=FALSE)
      } else {
        # bootstrap
        e.bs <- sample(res.demean,n+num.discard,replace=TRUE)
      }

      # create a new time series y.bs that mimicks y
      y.bs <- fitted(ls.fit) + e.bs[(1+num.discard):(n+num.discard)]
      # regression ls.fit
      if (is.na(x[1])) {
        if (CPP) {
          ls.fit.bs <- fLm(matrix(rep(1,n),ncol=1),y.bs)
        } else {
          ls.fit.bs <- lm(y.bs ~ 1)
        }
      } else {
        if (CPP) {
          ls.fit.bs <- fLm(cbind(1,x),y.bs)
        } else {
          ls.fit.bs <- lm(y.bs ~ x)
        }
      }

      # extract residuals
      res.bs <- ls.fit.bs$residuals

      # square the residuals
      res.sq.bs <- (res.bs-mean(res.bs))^2

      res.sq.demean.bs <- res.sq.bs-mean(res.sq.bs)

      # run autoregressive of a.order on res.sq
      if (CPP) {

      } else {

      }
      # ar.fit.bs <- lm(res.sq.demean.bs[(a.order+1):n]~res.sq.demean.bs[1:(n-a.order)])
      ar.fit.bs <- fLm(cbind(1,res.sq.demean.bs[1:(n-a.order)]), res.sq.demean.bs[(a.order+1):n])
      # AR(a.order) residuals
      ee.bs <- ar.fit.bs$residuals

      # calculate tau_hat (see the CJS paper)
      SSR0.bs <- sum(res.sq.demean.bs[(a.order+1):n]^2)
      # SSR1.bs <- sum(ee.bs[(a.order+1):n]^2)
      SSR1.bs <- sum(ee.bs^2)
      # this the observed statistic
      tau.bs[j] <- (n-a.order)*(1-SSR1.bs/SSR0.bs) # this the observed statistic
    }

    # associated p-value to report
    pValue.bs = mean(tau.bs > tau.hat)
    # associated critical value to report
    crit.bs = quantile(tau.bs, 1-sig.level)

    return(list(stat = tau.hat, crit.val = crit.bs, p.value = pValue.bs))
  }
}



if (!require(inline)) install.packages('inline')
library(inline)
src <- '
Rcpp::NumericMatrix Xr(Xs);
Rcpp::NumericVector yr(ys);
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
'
fLm <- cxxfunction(signature(Xs="numeric", ys="numeric"), src, plugin="RcppArmadillo")
