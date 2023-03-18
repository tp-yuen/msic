#' Compute the conditional expectation of uncured status
#'
#' Compute \eqn{\mathbb{E}[B | Y, \Delta, X, Z; \hat{\theta}, \hat{\varphi}]},
#' where \eqn{\hat{\theta}} and \eqn{\hat{\varphi}} are estimates from the previous iteration.
#'
#' @param delta A numeric vector (length \code{n}) of the censoring indicator.
#' @param p A numeric vector (length \code{n}) of the uncured probability values.
#' @param s.u A numeric vector (length \code{n}) of the
#' conditional survival function values, see \code{\link{estimate.survival}}.
#' @param y A numeric vector (length \code{n}) of the observed event time.
#' @param y.max y.max A positive number of the cure threshold.
#' Default is the maximum uncensored event time.
#'
#' @return A numeric vector of the conditional expected uncured status.
#' @keywords internal
uncure.status <- function(delta, p, s.u, y, y.max = max(y[delta == 1])) {
  b <- p * s.u
  b <- delta + (1 - delta) * b / (1 - p + b)
  b[delta == 0 & y > y.max] <- 0
  return(b)
}


#' Compute the conditional expectation of uncured status in the E-step of
#' the EM-algorithm
#'
#' @param y A numeric vector (length \code{n}) of the observed event time.
#' @param delta A numeric vector (length \code{n}) of the censoring indicator.
#' @param p A numeric vector (length \code{n}) of the uncured probability values.
#' @param s.u A numeric vector (length \code{n}) of the
#' conditional survival function values, see \code{\link{estimate.survival}}.
#' @param y.max A positive number of the cure threshold.
#' Default is the maximum uncensored event time.
#'
#' @return A numeric vector of the conditional expected uncured status.
#' @keywords internal
em.estep <- function(y, delta, p, s.u, y.max = max(y[delta == 1])) {
  w <- uncure.status(delta, p, s.u, y, y.max = y.max)
  return(w)
}
