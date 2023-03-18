#' Normalizing incidence regression parameter
#'
#' @param g A numeric vector (length \code{d}) of the incidence regression parameter.
#'
#' @return A numeric vector (length \code{d}) of the normalized regression parameter.
#' @keywords internal
normalize.gamma <- function(g) {
  g.norm <- sqrt(sum(g ^ 2))
  if (abs(abs(g.norm) - 1) > sqrt(.Machine$double.eps)) {
    g <- g / g.norm
  }
  return(g)
}


#' Split data into K folds
#'
#' @param n An integer of sample size.
#' @param k An integer of number of folds.
#' @param shuffle A logical if data is shuffled.
#'
#' @return A list of indices of each fold.
#' @keywords internal
kfold <- function(n, k, shuffle = FALSE) {
  folds <- ceiling(seq(n) * k / n)
  if (shuffle) {
    folds <- sample(folds)
  }
  return(unname(split(seq(n), folds)))
}


#' Kaplan-Meier estimate of survivor function
#'
#' @param y A numeric vector (length \code{n}) of the observed event time.
#' @param delta  A numeric vector (length \code{n}) of the censoring indicator.
#'
#' @return A \code{\link[stats]{stepfun}} of the estimate.
#' @importFrom stats stepfun
#' @export
estimate.km <- function(y, delta) {
  # Distinct observed event times in ascending order
  y.star <- sort(unique(y))
  # Risk set at t_{(j)} at the j-th row (n x length(y.star))
  risk.set <- sapply(y.star, function(x) {
    c(y >= x, # Risk set
      sum(y == x) - sum(delta[y == x] == 0) # No. of events excluding censored
    )
  })
  # No. of events at t_{(j)} at the j-th entry
  d.star <- drop(risk.set[nrow(risk.set), ])
  risk.set <- risk.set[-nrow(risk.set), ]

  # Hazard components
  hazard <- d.star / drop(colSums(risk.set))

  product.component <- cumprod(c(1, 1 - hazard))

  km.stepfun <- stepfun(y.star, product.component, f = 0, right = FALSE)
  return(km.stepfun)
}
