% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{data1}
\alias{data1}
\title{Simulated Data 1}
\format{
A data frame with 250 rows and 2 columns:
\describe{
  \item{X1,X2,X3,X4}{Incidence covariates}
  \item{Z1,Z2}{Latency covariates}
  \item{y}{Observed event time}
  \item{delta}{Censoring indicator}
}
}
\usage{
data("data1")
}
\description{
A simulated data of a mixture cure model with the following settings:
\describe{
  \item{Incidence regression parameter, \eqn{\gamma}}{
    \eqn{\gamma=(1.2, -0.2383, 0.7423, 0.3156, 0.5409)^T},
    where the first term is the intercept parameter;
  }
  \item{Latency regression parameter, \eqn{\beta}}{
    \eqn{\beta=(-0.8, 0.5)^T};
  }
  \item{Link function, \eqn{\varphi}}{
    \eqn{\varphi(x)=\frac{\exp(x)}{1 + \exp(x)}};
  }
  \item{Incidence covariates, \eqn{X}}{
    \eqn{(X_1, X_2, X_3, X_4)^T},
    where \eqn{
          X_1\sim U[0,1], X_2\sim N(0, 1), X_3\sim B(1, 0.3)} and
          \eqn{X_4\sim B(1, 0.6)}, and they are independent;
  }
  \item{Latency covariates, \eqn{Z}}{
    \eqn{(Z_1, Z_2)^T}, where \eqn{Z_1 = X_1, Z_2 = X_4};
  }
  \item{Survival time of the uncured subjects, \eqn{T_u}}{
    \eqn{T_u} follows the Weibull distribution with shape and scale
    parameters 2.2 and 1.5, respectively;
  }
  \item{Censoring time, \eqn{C}}{
    \eqn{C} follows the exponential distribution with rate parameter 0.1.
  }
}
}
\keyword{datasets}
