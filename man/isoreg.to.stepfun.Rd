% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/link.R
\name{isoreg.to.stepfun}
\alias{isoreg.to.stepfun}
\title{Construct a step function from \code{isoreg}}
\usage{
isoreg.to.stepfun(
  l,
  index.sorted,
  a,
  b,
  truncation.lower = 1e-06,
  truncation.upper = 1 - 1e-06
)
}
\arguments{
\item{l}{An object of class \code{\link[stats]{isoreg}}.}

\item{index.sorted}{A numeric vector (length \code{n}) of the sorted index,
\eqn{\gamma^TX}.}

\item{a}{A positive number of the estimated lower truncation,
see \code{\link{bounded.isoreg}}.}

\item{b}{A positive number of the estimated upper truncation,
see \code{\link{bounded.isoreg}}.}

\item{truncation.lower}{A positive number of the lower truncation parameter
for the link. Default is \code{1e-6}.}

\item{truncation.upper}{A positive number of the upper truncation parameter
for the link. Default is \code{1 - 1e-6}.}
}
\value{
A \code{\link[stats]{stepfun}} of the monotone link function estimate.
}
\description{
Construct a step function from \code{isoreg}
}
\keyword{internal}
