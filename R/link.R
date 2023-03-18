#' Integrated tri-weight kernel function
#'
#' Integrated tri-weight kernel function:
#' \deqn{\int_{-\infty}^x k(u)du,}
#' where \eqn{k(u) = \frac{35}{32}(1-u^2)^3I_{[-1,1]}(u)}.
#'
#' @param x A numeric vector of points for which the value of
#' `integrated.kern` is computed.
#'
#' @return A numeric vecotr of the computed `integrated.kern`.
#' @keywords internal
integrated.kern <- function(x){
  K <- ((35/32)*(x-x^3+3*x^5/5-x^7/7)+0.5)*ifelse(x>=-1,1,0)*ifelse(x<=1,1,0)+ifelse(x>1,1,0)
  return(K)
}


#' Kernel smoothing of monotone link function using tri-weight kernel
#'
#' @param x A numeric vector of points at which the smoothing is computed.
#' @param jump.x A numeric vector of points at which there is a jump.
#' @param jump.size A numeric vector of jump heights.
#' @param bw A positive number of the bandwidth parameter.
#'
#' @return A numeric vector of the smoothed link values at \code{x}.
#' @keywords internal
smoothing <- function(x, jump.x, jump.size, bw) {
  xt <- outer(x, jump.x, "-") / bw
  k.xt <- integrated.kern(xt)
  x.smoothed <- -1 * drop(jump.size %*% diff(t(k.xt)))
  return(x.smoothed)
}


#' Smooth monotone link function
#'
#' @description Smooth monotone link function, \eqn{\hat{\varphi}_{n,\theta}^s}
#'
#' @param index A numeric vector (length \code{n}) of index \eqn{\gamma^TX}.
#' @param mlink A \code{\link[stats]{stepfun}} of the monotone link estimate, see \code{\link{em.link}}.
#' @param bw A positive number of the bandwidth parameter.
#'
#' @return A numeric vector of the smoothed link values at \code{index}.
#' @keywords internal
#' @importFrom stats knots
smooth.monotone.link <- function(index, mlink, bw) {
  index <- drop(index)
  index.knots <- knots(mlink)
  p <- c(environment(mlink)$y, environment(mlink)$yright)

  # Extending the index range for smoothing
  index.knots.extended <- c(min(index) - 3 * bw, index.knots, max(index) + 3 * bw)
  smoothed.link <- smoothing(index, index.knots.extended, p, bw)
  if (any(smoothed.link < sqrt(.Machine$double.eps)) ||
      any(smoothed.link > (1 - sqrt(.Machine$double.eps)))) {
    writeLines("Smoothed link exceeding (0, 1)")
  }
  return(smoothed.link)
}


#' EM algorithm for monotone link estimator
#'
#' @description EM algorithm for monotone link estimator, \eqn{\hat{\varphi}_{n, \theta}}.
#'
#' @param y A numeric vector (length \code{n}) of the observed event time.
#' @param delta A numeric vector (length \code{n}) of the censoring indicator.
#' @param index A numeric vector (length \code{n}) of index, \eqn{\gamma^TX}.
#' @param uncured.survival A numeric vector (length \code{n}) of the
#' conditional survival function values, see \code{\link{estimate.survival}}.
#' @param p.init A numeric vector (length \code{n}) of the
#' initial uncured probability values.
#' @param truncation.lower A positive number of the lower truncation parameter
#' for the link. Default is \code{1e-6}.
#' @param truncation.upper A positive number of the upper truncation parameter
#' for the link. Default is \code{1 - 1e-6}.
#' @param y.max A positive number of the cure threshold.
#' Default is the maximum uncensored event time.
#' @param tol A positive number of convergence tolerance.
#' @param max.iter A positive integer of maximum number of iterations.
#'
#'
#' @return A list with components:
#' \describe{
#' \item{link}{A \code{\link[stats]{stepfun}} of the monotone link estimate.}
#' \item{convergence}{A logical if the EM algorithm converges.}
#' }
#' @keywords internal
#' @importFrom stats stepfun
#' @export
em.link <- function(y, delta, index, uncured.survival, p.init,
                    y.max = max(y[delta == 1]),
                    truncation.lower = 1e-6, truncation.upper = 1 - 1e-6,
                    tol = 1e-6, max.iter = 50L) {
  index <- drop(index)
  index.sorted <- sort(index)
  monotone.link <- NA
  p.previous <- rep(NA, length(index))
  p.current <- p.init

  for (k in seq(max.iter)) {
    w <- em.estep(y, delta, p.current, uncured.survival, y.max)
    l <- isotonized.link(length(index), index, w)

    knots.idx <- l$iKnots
    index.knots <- index.sorted[knots.idx]  # non-standardized index
    p.knot <- l$yf[knots.idx]
    # Truncation step
    p.knot <- pmax(truncation.lower, pmin(p.knot, truncation.upper))

    # Remove knots with duplicated slopes
    p.duplicated <- duplicated(p.knot, fromLast = TRUE)
    index.knots <- index.knots[!p.duplicated]
    p.knot <- p.knot[!p.duplicated]

    d <- p.knot
    monotone.link <- stepfun(index.knots, c(d, d[length(d)]),
                             right = TRUE, f = 1)  # non-standardized index
    p.previous <- p.current
    p.current <- monotone.link(index)
    p.succ.change <- sqrt(sum((p.current - p.previous) ^ 2))
    if (p.succ.change < tol) {
      return(list(link = monotone.link, convergence = TRUE))
    }
  }
  return(list(link = monotone.link, convergence = FALSE))
}


#' Expected prediction error for cure probability (EPECP)
#'
#' @param y A numeric vector (length \code{n}) of the observed event time.
#' @param delta A numeric vector (length \code{n}) of the censoring indicator.
#' @param p.est A numeric vector (length \code{n}) of the
#' estimated uncured probability.
#' @param revkm.stepfun A \code{\link[stats]{stepfun}} of the reversed Kaplan-Meier estimate of
#'                      the survivor function, see \code{\link{estimate.km}}.
#' @param y.max A positive number of the cure threshold.
#' Default is the maximum uncensored event time.
#' @param n An integer of the sample size.
#'
#' @references Jiang, W., Sun, H. and Peng, Y. (2017)
#' \emph{Prediction accuracy for the cure probabilities in mixture cure models,
#' Statistical Methods in Medical Research, Vol. 26, 2029-2041},
#' \doi{10.1177/0962280217708673}.
#'
#' @return A numeric of the computed EPECP.
#' @keywords internal
epecp <- function(y, delta, p.est, km.stepfun, y.max = max(y[delta == 1]),
                  n = length(y)) {
  sel <- as.integer((y >= y.max) | (delta == 1)) # Selected observations
  u <- as.integer((delta == 0) & (y >= y.max)) # Partially observed cure status
  epecp.val <- sel / km.stepfun(pmin(y, y.max))
  epecp.val <- epecp.val * (u - p.est) ^ 2
  epecp.val <- sum(epecp.val) / n
  return(epecp.val)
}


#' Construct a step function from \code{isoreg}
#'
#' @param l An object of class \code{\link[stats]{isoreg}}.
#' @param index.sorted A numeric vector (length \code{n}) of the sorted index,
#' \eqn{\gamma^TX}.
#' @param a A positive number of the estimated lower truncation,
#' see \code{\link{bounded.isoreg}}.
#' @param b A positive number of the estimated upper truncation,
#' see \code{\link{bounded.isoreg}}.
#' @param truncation.lower A positive number of the lower truncation parameter
#' for the link. Default is \code{1e-6}.
#' @param truncation.upper A positive number of the upper truncation parameter
#' for the link. Default is \code{1 - 1e-6}.
#'
#' @return A \code{\link[stats]{stepfun}} of the monotone link function estimate.
#' @keywords internal
isoreg.to.stepfun <- function(l, index.sorted, a, b, truncation.lower = 1e-6,
                              truncation.upper = 1 - 1e-6) {
  knots.idx <- l$iKnots
  index.knots <- index.sorted[knots.idx]  # non-standardized index
  p.knot <- l$yf[knots.idx]
  # Truncation step
  truncation.lower <- max(a, truncation.lower)
  truncation.upper <- min(b, truncation.upper)
  p.knot <- pmax(truncation.lower, pmin(p.knot, truncation.upper))

  # Remove knots with duplicated slopes
  p.duplicated <- duplicated(p.knot, fromLast = TRUE)
  index.knots <- index.knots[!p.duplicated]
  p.knot <- p.knot[!p.duplicated]

  d <- p.knot
  monotone.link <- stepfun(index.knots, c(d, d[length(d)]), right = TRUE, f = 1)
  return(monotone.link)
}


#' Bounded isotonic regression
#'
#' @param y.isoreg A numeric vector of the isotonic regressor.
#' @param y.mean A number of the sample mean of the observations.
#'
#' @references Luss, R. and Rosset, S. (2017)
#' \emph{Bounded isotonic regression, Electronic Journal of Statistics, Vol. 11, 4488-4514},
#' \doi{10.1214/17-EJS1365}.
#'
#' @return A list with components:
#' \describe{
#' \item{a}{A numeric vector of the lower range estimates along
#' the regularization path.}
#' \item{b}{A numeric vector of the upper range estimates along
#' the regularization path.}
#' \item{lambda}{A numeric vector of the regularization parameters
#' along the regularization path.}
#' }
#' @export
bounded.isoreg <- function(y.isoreg, y.mean) {
  a.k <- b.k <- y.mean
  y.k <- rep(y.mean, length(y.isoreg))
  lm.lower <- 2 * pmax(0, y.mean - y.isoreg) # Dual variable for lower range constraint
  lm.upper <- 2 * pmax(0, y.isoreg - y.mean) # Dual variable for upper range constraint
  lambda.k <- sum(lm.lower) # Regularization parameter
  y.isoreg.min <- min(y.isoreg)
  y.isoreg.max <- max(y.isoreg)
  a <- a.k; b <- b.k; lambda <- lambda.k
  while (a.k > y.isoreg.min && b.k < y.isoreg.max) {
    # Check active upper and lower constraints using complementary slackness.
    lower.active <- which(y.k == a.k & lm.lower > 0)
    upper.active <- which(y.k == b.k & lm.upper > 0)
    lm.lower.min <- min(lm.lower[lower.active])
    lm.upper.min <- min(lm.upper[upper.active])
    if (lm.lower.min * length(lower.active) < lm.upper.min * length(upper.active)) {
      a.k <- a.k - lm.lower.min / 2
      lm.lower <- 2 * pmax(0, a.k - y.isoreg)
      lambda.k <- sum(lm.lower)
      b.k <- (sum(y.isoreg[upper.active]) - lambda.k / 2) / length(upper.active)
      lm.upper <- 2 * pmax(0, y.isoreg - b.k)
    } else if (lm.lower.min * length(lower.active) > lm.upper.min * length(upper.active)) {
      b.k <- b.k + lm.upper.min / 2
      lm.upper <- 2 * pmax(0, y.isoreg - b.k)
      lambda.k <- sum(lm.upper)
      a.k <- (sum(y.isoreg[lower.active]) + lambda.k / 2) / length(lower.active)
      lm.lower <- 2 * pmax(0, a.k - y.isoreg)
    } else {
      a.k <- a.k - lm.lower.min / 2
      b.k <- b.k + lm.upper.min / 2
      lm.lower <- 2 * pmax(0, a.k - y.isoreg)
      lm.upper <- 2 * pmax(0, y.isoreg - b.k)
      lambda.k <- sum(lm.lower)
    }
    y.k <- pmax(a.k, pmin(y.isoreg, b.k))
    a <- c(a, a.k); b <- c(b, b.k); lambda <- c(lambda, lambda.k)
  }
  return(list(a = a, b = b, lambda = lambda))
}


#' Fit a isotonic regression for the monotone link estimator
#'
#' @param n An integer of the sample size.
#' @param index A numeric vector (length \code{n}) of index computed using all
#' observations for a fixed incidence regression parameter, \eqn{\gamma}.
#' @param w A numeric vector (length \code{n}) of the
#' conditional expected uncured status from the E-step, see \code{\link{em.estep}}.
#'
#' @return An object of class \code{\link[stats]{isoreg}}.
#' @keywords internal
#' @importFrom stats isoreg
isotonized.link <- function(n, index, w) {
  ordering <- order(index)
  x <- seq(1, n)
  y <- w[ordering]
  isoreg.fitted <- isoreg(x, y)
  return(isoreg.fitted)
}


#' EM algorithm for monotone link estimator using bounded isotonic regression
#'
#' @description EM algorithm for monotone link estimator,
#'     \eqn{\hat{\varphi}_{n, \theta}}, using bounded isotonic regression to
#'     determine the truncation.
#'
#' @param y A numeric vector (length \code{n}) of the observed event time.
#' @param delta A numeric vector (length \code{n}) of the censoring indicator.
#' @param index A numeric vector (length \code{n}) of index, \eqn{\gamma^TX}.
#' @param uncured.survival A numeric vector (length \code{n}) of the
#' conditional survival function values, see \code{\link{estimate.survival}}.
#' @param rev.km.fit A \code{\link[stats]{stepfun}} of the reverse Kaplan-Meier estimate from
#' \code{\link{estimate.km}}.
#' @param p.init A numeric vector (length \code{n}) of the
#' initial uncured probability values.
#' @param y.max A positive number of the cure threshold.
#' Default is the maximum uncensored event time.
#' @param n.fold An integer of the number folds used in the cross-validation for
#' determining the truncation. Default is \code{10}.
#' @param truncation.lower A positive number of the lower truncation parameter for the link.
#' Default is \code{1e-6}.
#' @param truncation.upper A positive number of the upper truncation parameter for the link.
#' Default is \code{1 - 1e-6}.
#' @param tol A positive number of convergence tolerance. Default is \code{1e-6}
#' @param max.iter A positive integer the of maximum number of iterations of the EM algorithm.
#' Default is \code{1L}.
#' @param shuffle A logical if data is shuffled. Default is \code{FALSE}.
#'
#'
#' @return A list with components:
#' \describe{
#' \item{link}{A \code{\link[stats]{stepfun}} of the link estimate.}
#' \item{convergence}{A logical if the EM algorithm converged.}
#' }
#' @importFrom stats stepfun
#' @export
em.link.bounded.isoreg <- function(y, delta, index, uncured.survival, rev.km.fit, p.init,
                                   y.max = max(y[delta == 1]), n.fold = 10,
                                   truncation.lower = 1e-6, truncation.upper = 1 - 1e-6,
                                   tol = 1e-6, max.iter = 1L, shuffle = FALSE) {
  index <- drop(index)
  n <- length(index)
  index.sorted <- sort(index)
  monotone.link <- NA
  p.previous <- rep(NA, length(index))
  p.current <- p.init

  folds <- kfold(n, n.fold, shuffle)
  for (k in seq(max.iter)) {
    # E-step
    w <- em.estep(y, delta, p.current, uncured.survival, y.max)

    # CV step
    epecp.list <- lapply(folds, function(idx) {
      # Bounded isotonic regression (BIR)
      l <- isotonized.link(n - length(idx), index[-idx], w[-idx])
      bisoreg.res <- bounded.isoreg(l$yf, mean(w[-idx]))
      bisoreg.res.df <- data.frame(bisoreg.res)

      # Compute EPECP for each fold
      epecp.val <- lapply(seq(nrow(bisoreg.res.df)), function(j) {
        monotone.link <- isoreg.to.stepfun(l, sort(index[-idx]),
                                           bisoreg.res.df[j, ][["a"]],
                                           bisoreg.res.df[j, ][["b"]])
        p.est <- 1 - monotone.link(index)
        return(c(
          lambda = bisoreg.res.df[j, ][["lambda"]],
          a = bisoreg.res.df[j, ][["a"]], b = bisoreg.res.df[j, ][["b"]],
          epecp = epecp(y[idx], delta[idx], p.est[idx], rev.km.fit,
                        max(y[delta == 1]), length(y))))
      })
      df <- data.frame(do.call(rbind, epecp.val))
      lambdas.order <- order(df[["lambda"]])
      epecp.vals <- df[["epecp"]][lambdas.order]
      epecp.vals <- c(epecp.vals[1], epecp.vals)
      stepfunc <- stepfun(df[["lambda"]][lambdas.order], epecp.vals)
      return(list(df = df, stepfunc = stepfunc))
    })

    # Combine all piece-wise EPECP against lambda
    lambdas.unique <- sort(unique(do.call(c, lapply(epecp.list, function(l) {
      return(l[["df"]][["lambda"]])
    }))))
    lambda.epecp.df <- lapply(lambdas.unique, function(lambda) {
      epecp.lambda <- sapply(epecp.list, function(l) {
        return(l[["stepfunc"]](lambda))
      })
      return(c(lambda = lambda, epecp = sum(epecp.lambda)))
    })
    lambda.epecp.df <- data.frame(do.call(rbind, lambda.epecp.df))
    # Select the left most lambda with minimum EPECP
    lambda.epecp.min <- lambda.epecp.df[which.min(lambda.epecp.df[["epecp"]]), "lambda"]

    # Plot of combined EPECP against lambda
    lambdas.order <- order(lambda.epecp.df[["lambda"]])
    epecp.vals <- lambda.epecp.df[["epecp"]][lambdas.order]
    epecp.vals <- c(epecp.vals[1], epecp.vals)
    stepfunc <- stepfun(lambda.epecp.df[["lambda"]][lambdas.order], epecp.vals)

    # BIR using the "optimal" lambda
    l <- isotonized.link(n, index, w)
    bisoreg.res <- bounded.isoreg(l$yf, mean(w))
    bisoreg.res.df <- data.frame(bisoreg.res)

    idx.selected <- findInterval(lambda.epecp.min, sort(bisoreg.res.df$lambda))
    idx.selected <- ifelse(idx.selected == nrow(bisoreg.res.df),
                           1L, ifelse(idx.selected == 0, nrow(bisoreg.res.df),
                                      nrow(bisoreg.res.df) - idx.selected + 1L))

    monotone.link <- isoreg.to.stepfun(l, sort(index),
                                       bisoreg.res.df[idx.selected, ][["a"]],
                                       bisoreg.res.df[idx.selected, ][["b"]],
                                       truncation.lower, truncation.upper)
    a <- max(bisoreg.res.df[idx.selected, ][["a"]], truncation.lower)
    b <- min(bisoreg.res.df[idx.selected, ][["b"]], truncation.upper)

    if (max.iter > 1L) {
      p.previous <- p.current
      p.current <- monotone.link(index)
      p.succ.change <- sqrt(sum((p.current - p.previous) ^ 2))
      if (p.succ.change < tol) {
        return(list(link = monotone.link, a = a, b = b))
      }
    }
  }
  return(list(link = monotone.link, a = a, b = b))
}
