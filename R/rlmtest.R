#' Robust Lagrange multiplier test for detecting ARCH/GARCH effect
#'
#' The function performs two resampling techniques to find critical values of the LM test, namely permutation and bootstrap.
#'
#' @param y  A numeric time series.
#' @param x  A time series of regressors.
#' @param a.order Order of the autoregressive ls.fit which must be a nonnegative integer number.
#' @param crit.type A string parameter allowing to choose "asymptotic" or "nonparametric" options.
#' @param nonp.method A string parameter allowing to choose "bootstrap" or "permutation" method of resampling
#' "bootstrap" - resampling of the estimated residuals (with replacement).
#' "permutation" - resampling of the estimated residuals (without replacement).
#' @param num.resamp Number of bootstrap replications if \code{"crit.type"}="nonparametric". Default number is 1,000.
#'
#' @return
#' @export
#'
#' @examples
#'
#' 1. ARCH effect
#' # model y_t = y_tb + e_t
#' # e_t ~ ARCH(1) coef
#' a.0 <- 0.6
#' a.1 <- 0.4
#' t-distribution
#' # z ~ t(5)        var = 5/(5-2)
#' n <- 500
#' z <- rt(n, df = 5)
#' e <- numeric(n)
#' e[1] <- rt(1, df = 5) * sqrt(a.0/(1.0-a.1) / (5/(5-2)))
#' # Generate ARCH(1) process
#' for(t in 2:n) {
#'   e[t] <- z[t]*sqrt(a.0+a.1*e[t-1]^2)
#' }
#' x <- runif(n,0,10)
#' y <- 0.25 + 0.5*x + e
#'
#' 2. independent series
#' # e_t ~ t(5) no ARCH effect
#' e <- rt(n, df = 5) * sqrt(a.0 / (5/(5-2)))
#' x <- runif(n,0,10)
#' y <- 0.25 + 0.5*x + e
#'
#' 3. real data
#' #' library(rbenchmark)
#' benchmark(rlmtest(y,x,crit.type = "nonparametric",nonp.method = "bootstrap",CPP = T),rlmtest(y,x,crit.type = "nonparametric",nonp.method = "bootstrap",CPP = F))
#'
#' data(GoldBTC)
#' x <- GoldBTC$Gold
#' y <- GoldBTC$BTC
#'
#' # test ARCH in Gold
#' rlmtest(y,x,crit.type = "nonparametric",nonp.method = "bootstrap",CPP = T)
#'
#' # test ARCH in bitcoin
#' rlmtest(y,x,crit.type = "nonparametric",nonp.method = "bootstrap",CPP = T)
#'
#' rlmtest(y,x,crit.type = "nonparametric",nonp.method = "bootstrap",CPP = T)
#'
# $`stat`
# [1] 0.6353227
#
# $crit.val
# 95%
# 3.109865
#
# $p.value
# [1] 0.262
#
# > rlmtest(x,y,crit.type = "nonparametric",nonp.method = "bootstrap",CPP = T)
# $`stat`
# [1] 41.62449
#
# $crit.val
# 95%
# 2.006489
#
# $p.value
# [1] 0.001
#
# > rlmtest(x,crit.type = "nonparametric",nonp.method = "bootstrap",CPP = T)
# $`stat`
# [1] 32.4635
#
# $crit.val
# 95%
# 2.469525
#
# $p.value
# [1] 0.001
#
# > rlmtest(y,crit.type = "nonparametric",nonp.method = "bootstrap",CPP = T)
# $`stat`
# [1] 0.6009142
#
# $crit.val
# 95%
# 2.672682
#
# $p.value
# [1] 0.246
#
# > plot(1:361,x)
# > plot(1:361,y)
# > plot(1:361,x)
# > plot(1:361,y)
# > rlmtest(y,crit.type = "asymptotic",CPP = T)
# $`stat`
# [1] 0.6009142
#
# $crit.val
# [1] 3.841459
#
# $p.value
# [1] 0.4382294
#
# > rlmtest(x,crit.type = "asymptotic",CPP = T)
# $`stat`
# [1] 32.4635
#
# $crit.val
# [1] 3.841459
#
# $p.value
# [1] 1.214525e-08

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
      ls.fit <- fLM(matrix(rep(1,n),ncol=1),y)
    } else {
      ls.fit <- lm(y ~ 1)
    }
  } else {
    if (CPP) {
      ls.fit <- fLM(cbind(1,x),y)
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
    ar.fit <- fLM(cbind(1,res.sq.demean[1:(n-a.order)]), res.sq.demean[(a.order+1):n])
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
          ls.fit.bs <- fLM(matrix(rep(1,n),ncol=1),y.bs)
        } else {
          ls.fit.bs <- lm(y.bs ~ 1)
        }
      } else {
        if (CPP) {
          ls.fit.bs <- fLM(cbind(1,x),y.bs)
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
      if (1) {
        # ar.fit.bs <- lm(res.sq.demean.bs[(a.order+1):n]~res.sq.demean.bs[1:(n-a.order)])
        ar.fit.bs <- fLM(cbind(1,res.sq.demean.bs[1:(n-a.order)]), res.sq.demean.bs[(a.order+1):n])
      } else {
        ar.fit.bs <- lm(res.sq.demean.bs[(a.order+1):n]~res.sq.demean.bs[1:(n-a.order)])
      }

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

