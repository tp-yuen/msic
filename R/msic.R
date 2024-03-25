#' Fit a monotone single-index mixture cure model
#'
#' Fit a monotone single-index mixture cure model
#'
#' @param X A numeric matrix (dimension \code{n} x \code{d};
#' each row is an observation vector) of covariates for the incidence (no intercept term).
#' @param Z A numeric matrix (dimension \code{n} x \code{q};
#' each row is an observation vector) of covariates for the latency (no intercept term).
#' @param y A numeric vector (length \code{n}) of the observed event time.
#' @param delta A numeric vector (length \code{n}) of the censoring indicator.
#' @param bw A positive number of bandwidth, \code{NULL} indicates using
#' \eqn{rn^{-1/5}}, where \eqn{r} is the range of the index \eqn{\gamma^TX}.
#' @param truncation.lower A positive number of the lower truncation parameter for the link.
#' Default is \code{1e-6}. See Details below.
#' @param truncation.upper A positive number of the upper truncation parameter for the link.
#' Default is \code{1 - 1e-6}. See Details below.
#' @param max.iter An integer of the maximum number of iterations.
#' Default is \code{50L}.
#' @param tol A positive number of convergence tolerance. Default is \code{1e-5}.
#' @param truncation.est A logical if the upper and lower truncation are to be
#' determined by the bounded isotonic regression,
#' see \code{\link{em.link.bounded.isoreg}}. See also Details below.
#' @param gamma.init A numeric vector (length \code{d}) of the
#' initial incidence regression parameter value.
#' \code{NULL} indicates it is initialized by fitting a logistic regression of
#' the censoring indicator (\code{delta}) against the covariates \code{X}.
#' Default is \code{NULL}.
#' @param beta.init A numeric vector (length \code{q}) of the
#' initial latency regression parameter value.
#' \code{NULL} indicates it is initialized by fitting a Cox PH model
#' \code{\link[survival]{coxph}} to the uncensored observations with covariates \code{Z}.
#' Default is \code{NULL}.
#'
#' @details If \code{truncation.lower} or \code{truncation.upper} is \code{NULL},
#' then the two truncation parameters are determined from
#' the logistic regression of the censoring indicator, \code{delta}, against
#' the covariates, \code{X}, in the initialization step.
#' If \code{truncation.est} is \code{TRUE}, the two truncation parameters are
#' determined by the bounded isotonic regression, see \code{\link{em.link.bounded.isoreg}}.
#'
#' @examples
#' # Using default truncation and simulated Data 1 ----------------------------
#' data("data1")
#' X <- cbind(data1$X1, data1$X2, data1$X3, data1$X4)
#' Z <- cbind(data1$Z1, data1$Z2)
#' # Model estimation
#' m.fit.1 <- msic(X, Z, data1$y, data1$delta, truncation.est = FALSE)
#' # Print the estimated regression parameters
#' print(m.fit.1)
#' # Plot the estimated link functions
#' plot(m.fit.1)
#' # Predict uncured probability
#' X.new <- cbind(
#'   matrix(runif(5), ncol = 1),
#'   matrix(rnorm(5, 0, 1), ncol = 1),
#'   matrix(rbinom(5, 1, 0.3), ncol = 1),
#'   matrix(rbinom(5, 1, 0.6), ncol = 1)
#' )
#' p.hat.1 <- predict(m.fit.1, X.new)
#'
#' # Using truncation being determined by the bounded isotonic regression
#' # method and simulated Data 2 ----------------------------
#' data("data2")
#' X <- cbind(data2$X1, data2$X2, data2$X3, data2$X4)
#' Z <- cbind(data2$Z1, data2$Z2)
#' # Model estimation
#' m.fit.2 <- msic(X, Z, data2$y, data2$delta, truncation.est = TRUE)
#' # Print the estimated regression parameters
#' print(m.fit.2)
#' # Plot the estimated link functions
#' plot(m.fit.2)
#' # Predict uncured probability
#' X.new <- cbind(
#'   matrix(runif(5), ncol = 1),
#'   matrix(rnorm(5, 0, 1), ncol = 1),
#'   matrix(rbinom(5, 1, 0.3), ncol = 1),
#'   matrix(rbinom(5, 1, 0.6), ncol = 1)
#' )
#' p.hat.2 <- predict(m.fit.2, X.new)
#'
#' # Using pre-specified truncation and simulated Data 2 ----------------------------
#' data("data2")
#' X <- cbind(data2$X1, data2$X2, data2$X3, data2$X4)
#' Z <- cbind(data2$Z1, data2$Z2)
#' # Model estimation with pre-specified truncation
#' m.fit.3 <- msic(X, Z, data2$y, data2$delta, truncation.est = FALSE,
#'                 truncation.lower = 0.5, truncation.upper = 0.75)
#' # Print the estimated regression parameters
#' print(m.fit.3)
#' # Plot the estimated link functions
#' plot(m.fit.3)
#' # Predict uncured probability
#' X.new <- cbind(
#'   matrix(runif(5), ncol = 1),
#'   matrix(rnorm(5, 0, 1), ncol = 1),
#'   matrix(rbinom(5, 1, 0.3), ncol = 1),
#'   matrix(rbinom(5, 1, 0.6), ncol = 1)
#' )
#' p.hat.3 <- predict(m.fit.3, X.new)
#'
#' @references Musta, E. and Yuen, T. P. (2022)
#' \emph{Single-index mixture cure model under monotonicity constraints, arXiv preprint},
#' \doi{10.48550/ARXIV.2211.09464}.
#'
#' @return An object of class \code{msic}.
#' @importFrom survival coxph Surv
#' @importFrom stats glm.fit binomial
#' @importFrom stats as.formula
#' @export
msic <- function(X, Z, y, delta,
                 bw = NULL,
                 truncation.lower = 1e-6, truncation.upper = 1 - 1e-6,
                 max.iter = 50L, tol = 1e-5, truncation.est = FALSE,
                 gamma.init = NULL, beta.init = NULL) {
  y <- drop(y)
  delta <- drop(delta)
  n <- length(y)
  y.max <- max(y[delta == 1])

  # Initialization
  logistic.fitted <- numeric(ncol(X) + 1L)
  if (is.null(gamma.init)) {
    ## Link initialization using logistic regression on censoring indicator
    logistic.fitted <- glm.fit(cbind(rep(1, n), X), delta,
                               family = binomial(link = logit), intercept = FALSE)
    logistic.fitted <- unname(logistic.fitted$coefficients)
  } else {
    if (length(gamma.init) == ncol(X)) {
      logistic.fitted <- c(0, drop(gamma.init))
    } else {
      logistic.fitted <- drop(gamma.init)
    }
  }
  ## Compute uncured probability
  p.current <- numeric(n)
  p.current <- drop(X %*% logistic.fitted[-1] + logistic.fitted[1])
  p.current <- exp(p.current) / (1 + exp(p.current))
  gamma.current <- unname(logistic.fitted[-1])

  ## Set the upper and lower truncation parameters if
  ## both truncation.lower and truncation.upper are NULL
  if (any(is.null(c(truncation.lower, truncation.upper)))) {
    truncation.lower <- max(min(p.current), 1e-6)
    truncation.upper <- min(max(p.current), 1 - 1e-6)
  }

  beta.current <- numeric(ncol(Z))
  if (is.null(beta.init)) {
    ## Latency regression parameter initialization using
    ## Cox PH model with the uncensored data
    beta.df <- as.data.frame(cbind(y, delta, Z))
    colnames(beta.df) <- c("y", "delta", paste("Z", seq(ncol(Z)), sep = ""))
    cox.ph.formula <- as.formula(
      paste0("Surv(y, delta) ~ ", paste0("Z", seq(ncol(Z)), collapse = " + ")))
    beta.current <- unname(coxph(cox.ph.formula, subset = delta != 0,
                                 data = beta.df, method = "breslow")$coef)
  } else {
    beta.current <- beta.init
  }
  ## Conditional survival function evaluated at all observations
  s.current <- estimate.survival(y, delta, Z, beta.current, delta)
  s.u.current <- s.current$s.u
  s.baseline.current <- s.current$s.baseline

  ## Monotone link estimation 0-th iteration
  ### Normalize gamma
  gamma.current <- normalize.gamma(gamma.current)
  index.current <- drop(X %*% gamma.current)

  l.current <- NULL
  if (truncation.est) {
    # Using Bounded isotonic regression to determine the thresholds
    rev.km.fit <- estimate.km(y, 1L - delta)
    em.bisoreg.link <- em.link.bounded.isoreg(y, delta, index.current,
                                              s.u.current, rev.km.fit, p.current)
    l.current <- em.bisoreg.link
    truncation.lower <- em.bisoreg.link$a
    truncation.upper <- em.bisoreg.link$b
  } else {
    # Using previous link and current gamma
    l.current <- em.link(y, delta, index.current, s.u.current, p.current, y.max,
                         truncation.lower, truncation.upper)
  }
  monotone.link.current <- l.current$link

  # Smoothing the monotone link
  smoothed.link.current <- function(x, b) {
    return(smooth.monotone.link(x, monotone.link.current, b))
  }
  bw.current <- bw
  if (is.null(bw)) {
    bw.current <- diff(range(index.current)) / (length(index.current) ^ 0.2)
  }
  # Using current link and current gamma
  p.current <- smoothed.link.current(index.current, bw.current)

  # Estimates in the previous iteration
  gamma.previous <- rep(NA, length(gamma.current))
  p.previous <- rep(NA, length(p.current))
  beta.previous <- rep(NA, length(beta.current))
  s.previous <- s.current
  s.u.previous <- rep(NA, length(s.u.current))
  s.baseline.previous <- rep(NA, length(s.baseline.current))
  bw.previous <- bw.current
  monotone.link.previous <- monotone.link.current
  smoothed.link.previous <- smoothed.link.current

  # Trace of variables
  gamma.trace <- list(gamma.current)
  p.trace <- list(p.current)
  beta.trace <- list(beta.current)
  s.u.trace <- list(s.u.current)
  s.baseline.trace <- list(s.baseline.current)
  monotone.link.trace <- list(monotone.link.current)
  link.trace <- list(smoothed.link.current)
  succ.change.trace <- NULL
  bw.trace <- c(bw.current)

  output <- NULL
  # Model estimation iterations
  for (m in seq(max.iter)) {
    writeLines(sprintf("Iteration: %d", m))
    gamma.previous <- gamma.current
    p.previous <- p.current
    beta.previous <- beta.current
    s.previous <- s.current
    s.u.previous <- s.u.current
    s.baseline.previous <- s.baseline.current
    bw.previous <- bw.current
    monotone.link.previous <- monotone.link.current
    smoothed.link.previous <- smoothed.link.current

    # EM algorithm for \theta
    ## E-step of EM for \theta
    w.current <- em.estep(y, delta, p.previous, s.u.previous, y.max)

    ## M-step of EM for \theta
    ### M-step for \gamma
    incidence <- estimate.incidence(X, y, delta, w.current, s.u.current,
                                    p.previous, gamma.previous, bw.previous,
                                    y.max, truncation.lower, truncation.upper)
    if (incidence$convergence != 0) {
      stop(sprintf(paste("Algorithm does not converge when estimating gamma",
                         "with error code: %d"), incidence$convergence))
    }
    gamma.current <- incidence$gamma
    l.current <- incidence$l
    monotone.link.current <- l.current$link
    smoothed.link.current <- incidence$smoothed.link

    ### M-step for \beta
    beta.est <- estimate.beta(y, delta, Z, w.current, beta.previous)
    beta.current <- drop(beta.est$beta)

    ### M-step for \Lambda
    s.current <- estimate.survival(y, delta, Z, beta.current, w.current)
    s.u.current <- s.current$s.u
    s.baseline.current <- s.current$s.baseline

    # Update the link
    ## Compute uncured probabilities using previous link and current gamma
    index.current <- drop(X %*% gamma.current)

    ## Update the cure probabilities using current link and current gamma
    bw.current <- bw
    if (is.null(bw)) {
      bw.current <- diff(range(index.current)) / (length(index.current) ^ 0.2)
    }
    p.current <- smoothed.link.current(index.current, bw.current)

    # Compute successive changes
    succ.change <- sum(c(beta.current - beta.previous,
                         s.baseline.current - s.baseline.previous,
                         gamma.current - gamma.previous) ^ 2)

    # Update trace of variables
    gamma.trace <- c(gamma.trace, list(gamma.current))
    p.trace <- c(p.trace, list(p.current))
    beta.trace <- c(beta.trace, list(beta.current))
    s.u.trace <- c(s.u.trace, list(s.u.current))
    s.baseline.trace <- c(s.baseline.trace, list(s.baseline.current))
    monotone.link.trace <- c(monotone.link.trace, list(monotone.link.current))
    link.trace <- c(link.trace, list(smoothed.link.current))
    bw.trace <- c(bw.trace, bw.current)
    succ.change.trace <- c(succ.change.trace, succ.change)

    output <- list(gamma = gamma.current, beta = beta.current,
                   s.baseline = s.baseline.current,
                   monotone.link = monotone.link.current,
                   link = smoothed.link.current,
                   p = p.current, bw = bw.current,
                   X = X, Z = Z, y = y, delta = delta, y.max = y.max)
    if (succ.change < tol) {
      writeLines(sprintf("EM algorithm converges at the %d-th iteration", m))
      output <- c(output, list(convergence = TRUE))
      break()
    }
  }
  if (!("convergence" %in% names(output))) {
    writeLines(sprintf("EM algorithm reaches maximum number of iterations."))
    output <- c(output, list(convergence = FALSE))
  }
  class(output) <- "msic"
  return(output)
}


#' Print a msic object
#'
#' Print a summary of the fitted monotone single-index cure model.
#'
#' @details
#' The estimated incidence regression parameter \eqn{\hat{\gamma}} is printed,
#' followed by the estimated latency regression parameter \eqn{\hat{\beta}}.
#'
#' @param x An \code{msic} object.
#' @param digits An integer of significant digits in printout.
#' @param \dots Additional print arguments.
#' @return A printout mentioned in Details.
#' @method print msic
#' @export
print.msic <- function (x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("Estimated incidence regression parameter:\n")
  x.colnames <- colnames(x$X)
  if (is.null(x.colnames)) {
    x.colnames <- paste("X", seq(ncol(x$X)), sep = "")
  }
  gamma.out <- x$gamma
  names(gamma.out) <- x.colnames
  print(format(gamma.out, digit = digits))
  cat("\n")

  cat("Estimated latency regression parameter:\n")
  z.colnames <- colnames(x$Z)
  if (is.null(z.colnames)) {
    z.colnames <- paste("Z", seq(ncol(x$Z)), sep = "")
  }
  beta.out <- x$beta
  names(beta.out) <- z.colnames
  print(format(beta.out, digit = digits))
  cat("\n")
}


#' Plot a msic object
#'
#' Plot the estimated monotone link functions.
#'
#' @details
#' The estimated monotone link function \eqn{\hat{\varphi}_{n,\theta}} and
#' the smoothed monotone link function \eqn{\hat{\varphi}_{n,\theta}^s} are plotted.
#'
#' @param x An \code{msic} object.
#' @param \dots Additional plot arguments.
#' @return A plot mentioned in Details.
#' @importFrom graphics lines legend
#' @method plot msic
#' @export
plot.msic <- function(x, ...) {
  index.range <- range(drop(x$X %*% x$gamma))
  index.plot <- seq(index.range[1] - 3 * x$bw, index.range[2] + 3 * x$bw, 0.01)
  plot.title <- "Estimated link functions"
  plot(x = index.plot, y = x$monotone.link(index.plot), type = 'l',
       col = 2, lwd = 2, lty = 1, main = plot.title,
       xlab = "Index", ylab = "Uncured Probability", ylim = c(0, 1))
  lines(x = index.plot, y = x$link(index.plot, x$bw), col = 3, lwd = 2, lty = 2)
  legend(x = "bottomright", legend = c("Monotone", "Smoothed"),
         lty = c(1, 2), col = c(2, 3), lwd = 2, xpd = TRUE, bty = "n")
}


#' Predict uncure probability
#'
#' Predict the uncure probability given covariate \eqn{X}.
#'
#'
#' @param object An \code{msic} object.
#' @param X A numeric matrix of incidence covariates; each row is an observation vector.
#' @param \dots Additional arguments.
#' @return A numeric vector of the probability.
#' @method predict msic
#' @export
predict.msic <- function(object, X, ...) {
  index <- drop(X %*% object$gamma)
  p.hat <- object$link(index, object$bw)
  return(p.hat)
}
