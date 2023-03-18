#' Estimate the latency regression parameter
#'
#' @description Estimate the latency regression parameter, \eqn{\beta}.
#'
#' @param y A numeric vector (length \code{n}) of the observed event time.
#' @param delta A numeric vector (length \code{n}) of the censoring indicator.
#' @param Z A numeric matrix (dimension \code{n} x \code{q};
#' each row is an observation vector) of covariates for the latency (no intercept term).
#' @param w A numeric vector (length \code{n}) of the
#' conditional expected uncured status from the E-step, see \code{\link{em.estep}}.
#' @param beta.init A numeric vector (length \code{q}) of the
#' initial latency regression parameter value. Default is a zero vector.
#'
#' @return A list with components: \describe{
#'   \item{beta}{A vector (length \code{q}) of the latency regression parameter estimate.}
#'   \item{fit}{An object returned from \code{\link[survival]{coxph}}.}
#' }
#'
#' @importFrom survival coxph
#' @keywords internal
estimate.beta <- function(y, delta, Z, w, beta.init = rep(0, ncol(Z))) {
  coxph.fitted <- coxph(Surv(y, event = delta) ~ Z + offset(log(w)), init = beta.init,
                        subset = (w != 0), method = "breslow")
  return(list(beta = as.matrix(unname(coxph.fitted$coef)), fit = coxph.fitted))
}

#' Estimate the conditional baseline cumulative hazard function
#'
#' Estimate the conditional baseline cumulative hazard function, \eqn{\Lambda(t)}.
#'
#' @param y A numeric vector (length \code{n}) of the observed event time.
#' @param delta A numeric vector (length \code{n}) of the censoring indicator.
#' @param Z A numeric matrix (dimension \code{n} x \code{q};
#' each row is an observation vector) of covariates for the latency (no intercept term).
#' @param beta A numeric vector (length \code{q}) of the
#' initial latency regression parameter value.
#' @param w  A numeric vector (length \code{n}) of the
#' conditional expected uncured status from the E-step, see \code{\link{em.estep}}.
#'
#' @return A list with components: \describe{
#'   \item{y.star}{A vector of the distinct uncensored event times in ascending order.}
#'   \item{hazard}{A vector of the conditional baseline hazard estimates.}
#'   \item{cum.hazard}{A vector of the conditional baseline cumulative hazard estimates.}
#' }
#'
#' @keywords internal
baseline.hazard <- function(y, delta, Z, beta, w) {
  # Distinct uncensored event times in ascending order
  y.star <- sort(unique(y[delta == 1]))
  # Risk set at t_{(j)} at the j-th row (n x length(y.star))
  risk.set <- sapply(y.star, function(x) c(y >= x, sum(y == x)))
  # No. of events at t_{(j)} at the j-th entry
  d.star <- drop(risk.set[nrow(risk.set), ])
  risk.set <- risk.set[-nrow(risk.set), ]

  # w_{i} * \exp(Z_{i} * \beta) at the i-th entry
  w.exp.zbeta <- w * exp(Z %*% beta)
  # \sum_{i \in R_{t_{(j)}}} w_{i} * \exp(Z_{i} * \beta)
  w.exp.zbeta.sum <- drop(crossprod(risk.set, w.exp.zbeta))
  # Baseline cumulative hazard estimates at the distinct times.
  hazard <- d.star / w.exp.zbeta.sum
  cum.hazard <- drop(risk.set %*% hazard)
  cum.hazard[y < min(y.star)] <- 1
  cum.hazard[y > max(y.star)] <- 0

  y.interval <- findInterval(y, y.star, left.open = TRUE) + 1L
  hazard.extend <- c(hazard, NA)
  hazard.out <- hazard.extend[y.interval]
  return(list(y.star = y.star, hazard = hazard.out, cum.hazard = cum.hazard))
}


#' Estimate the conditional survival function for the latency
#'
#' @description Estimate the conditional survival function for the latency, \eqn{S_u(t|Z)}.
#'
#' @param A numeric vector (length \code{n}) of the observed event time.
#' @param delta A numeric vector (length \code{n}) of the censoring indicator.
#' @param Z A numeric matrix (dimension \code{n} x \code{q};
#' each row is an observation vector) of covariates for the latency (no intercept term).
#' @param beta A numeric vector (length \code{q}) of the
#' initial latency regression parameter value.
#' @param w A numeric vector (length \code{n}) of the
#' conditional expected uncured status from the E-step, see \code{\link{em.estep}}.
#'
#' @return A list with components:
#' \describe{
#'     \item{s.u}{A vector of the conditional survival function values.}
#'     \item{s.baseline}{A vector of the conditional baseline survival function values.}
#'     \item{baseline}{The list return from \code{\link{baseline.hazard}}}
#' }
#' @keywords internal
estimate.survival <- function(y, delta, Z, beta, w) {
  baseline <- baseline.hazard(y, delta, Z, beta, w)
  s.baseline <- baseline$cum.hazard
  s.baseline[y > max(baseline$y.star)] <- Inf
  s.baseline[y < min(baseline$y.star)] <- 0
  s.baseline <- exp(-1 * s.baseline)
  s.u <- s.baseline ^ (exp(Z %*% beta))
  return(list(s.u = drop(s.u), s.baseline = drop(s.baseline), baseline = baseline))
}
