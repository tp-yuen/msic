#' Log-likelihood function of the incidence part
#'
#' @description Log-likelihood function of the incidence part:
#'     \deqn{\sum_{i=1}^{n}w_i\log(p_i) + (1-w_i)\log(1-p_i).}
#'
#' @param p A numeric vector (length \code{n}) of the uncured probabilities.
#' @param w A numeric vector (length \code{n}) of the
#' conditional expected uncured status from the E-step, see \code{\link{em.estep}}.
#'
#' @return The computed log-likelihood value.
#' @keywords internal
llh.incidence <- function(p, w) {
  w.is.zero <- (w == 0)
  w.is.one <- (w == 1)
  w.is.in.zero.one <- !(w.is.zero | w.is.one)

  llh.single <- c(log(1 - p[w.is.zero]),
                  log(p[w.is.one]),
                  w[w.is.in.zero.one] * log(p[w.is.in.zero.one]) +
                    (1 - w[w.is.in.zero.one]) * log(1 - p[w.is.in.zero.one]))
  llh <- sum(llh.single)
  return(llh)
}


#' Maximize the log-likelihood of the incidence part w.r.t. regression parameter
#' using the Augmented Lagrangian method
#'
#' @description Maximize the log-likelihood of the incidence part w.r.t. \eqn{\gamma}.
#'
#' @param X A numeric matrix (dimension \code{n} x \code{d};
#' each row is an observation vector) of covariates for the incidence
#' (no intercept term).
#' @param y A numeric vector (length \code{n}) of the observed event time.
#' @param delta A numeric vector (length \code{n}) of the censoring indicator.
#' @param w A numeric vector (length \code{n}) of the
#' conditional expected uncured status from the E-step, see \code{\link{em.estep}}.
#' @param s.u.current A numeric vector (length \code{n}) of the
#' conditional survival function values, see \code{\link{estimate.survival}}.
#' @param p.init A numeric vector (length \code{n}) of the
#' initial uncured probability values.
#' @param gamma.init A numeric vector (length \code{d}) of the
#' initial incidence regression parameter value.
#' @param bw A positive number of the bandwidth parameter for smoothing.
#' @param y.max A positive number of the cure threshold.
#' @param truncation.lower A positive number of the lower truncation parameter for the link.
#' @param truncation.upper A positive number of the upper truncation parameter for the link.
#'
#' @return A list with components:
#' \describe{
#' \item{gamma}{A numeric vector (length \code{d}) of the incidence
#' regression parameter estimate.}
#' \item{convergence}{A integer code returned by \code{\link[Rsolnp]{solnp}}.}
#' \item{optim}{A list returned by \code{\link[Rsolnp]{solnp}}}
#' \item{l}{A \code{\link[stats]{stepfun}} of the monotone link estimate from the EM-algorithm,
#' see \code{\link{em.link}}.}
#' \item{smoothed.link}{A function of the smoothed monotone link estimate,
#' see \code{\link{smooth.monotone.link}}.}
#' }
#'
#' @importFrom Rsolnp solnp
#' @keywords internal
estimate.incidence <- function(X, y, delta, w, s.u.current, p.init, gamma.init,
                               bw, y.max, truncation.lower, truncation.upper) {
  n <- length(w)
  ### Normalize gamma
  gamma.init <- normalize.gamma(gamma.init)
  ### Environment to store the link estimates
  est.env <- new.env()

  # Define the objective function
  obj.func <- function(g) {
    if (any(abs(g) > 10)) {
      return(1e8)
    }
    index <- drop(X %*% g)
    l.current <- em.link(y, delta, index, s.u.current, p.init, y.max,
                         truncation.lower, truncation.upper)
    monotone.link.current <- l.current$link

    ### Smoothing the monotone link
    smoothed.link.current <- function(x, b) {
      return(smooth.monotone.link(x, monotone.link.current, b))
    }
    p <- smoothed.link.current(index, bw)
    assign("l", l.current, envir = est.env)
    assign("smoothed.link", smoothed.link.current, envir = est.env)
    return(-1.0 * llh.incidence(p, w) / length(w))
  }
  # Define the identifiability constraint
  eq.con <- function(g) sum(g ^ 2)

  res <- NA
  # Solve the equality constrained maximization problem
  tryCatch({
    control <- list(delta = 1e-8, tol = 1e-10)
    res <- solnp(gamma.init, obj.func, eqfun = eq.con, eqB = 1.0, control = control)
    if (any(abs(diff(res$values)) >= 0.01) || res$convergence != 0) {
      assign("bw", bw * 1.5, envir = environment(fun = obj.func))
      res <- solnp(gamma.init, obj.func, eqfun = eq.con, eqB = 1.0)
    }
    return(list(gamma = res$par, convergence = res$convergence, optim = res,
                l = get("l", envir = est.env),
                smoothed.link = get("smoothed.link", envir = est.env)))
  }, error = function(e) {
    print(e)
    stop(e)
  })
}
