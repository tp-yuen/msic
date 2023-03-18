# Monotone single-index mixture cure model

Mixture cure models are widely utilized to analyze time-to-event data in the presence of a cure fraction. Such models consist of two sub-models: one for the probability of being uncured (incidence) and one for the survival of the uncured subjects (latency). One of the most commonly used model is the logistic/Cox mixture cure model, which assumes a logistic model for the incidence and a Cox proportional hazard model for the latency. The R package [`smcure`](https://cran.r-project.org/web/packages/smcure/) implements the Expectation-Maximization (EM) algorithm for the maximum likelihood estimation for the logistic/Cox mixture cure model (Cai et al., 2012). Amico et al. (2019) introduce a single-index model for the incidence, which relaxes the parametric assumption for the incidence and circumvents the curse-of-dimensionality that non-parametric estimators (Xu and Peng, 2014) possess when the covariates are multidimensional.

We propose a monotone single-index model for the incidence and introduce a new estimation method that is based on the maximum likelihood approach and techniques from isotonic regression. Specifically, the estimation procedure consists of three main steps:

1.  Maximizing the likelihood function w.r.t. the link function parameter, when the remaining parameters are fixed, to obtain the monotone link estimator.

2.  Kernel smoothing the monotone link estimator to obtain the smooth monotone link estimator.

3.  Plugging the smooth monotone link in the likelihood function and maximizes it w.r.t. the remaining parameters.

The maximization problems in steps 1 and 3 are solved via the EM algorithm. Detailed algorithms of the implementation and the theoretical properties of the method can be found in Musta and Yuen (2022). The implementation is packed as an R package `msic` in this repository.

## Installation
Install the package from GitHub via the [**remotes**](https://remotes.r-lib.org) package:
```R
remotes::install_github('tp-yuen/msic')
```

## Example
After installing the package, one can run the following code based on a simulated dataset as a demonstration. More examples are provided in the package documentation.
```R
# Using default truncation and simulated Data 1 ----------------------------
library(msic)
data("data1")
X <- cbind(data1$X1, data1$X2, data1$X3, data1$X4)
Z <- cbind(data1$Z1, data1$Z2)
# Model estimation
m.fit.1 <- msic(X, Z, data1$y, data1$delta, truncation.est = FALSE)
# Print the estimated regression parameters
print(m.fit.1)
# Plot the estimated link functions
plot(m.fit.1)
# Predict uncured probability
X.new <- cbind(
  matrix(runif(5), ncol = 1),
  matrix(rnorm(5, 0, 1), ncol = 1),
  matrix(rbinom(5, 1, 0.3), ncol = 1),
  matrix(rbinom(5, 1, 0.6), ncol = 1)
)
p.hat.1 <- predict(m.fit.1, X.new)
```

## References

Amico, M., Van Keilegom, I., & Legrand, C. (2019). The single‐index/Cox mixture cure model. Biometrics, 75(2), 452--462. 

Cai, C., Zou, Y., Peng, Y., & Zhang, J. (2012). smcure: An R-package for estimating semiparametric mixture cure models. Computer Methods and Programs in Biomedicine, 108(3), 1255--1260. 

Musta, E., & Yuen, T. P. (2022). Single-index mixture cure model under monotonicity constraints. 

Xu, J., & Peng, Y. (2014). Nonparametric cure rate estimation with covariates. Canadian Journal of Statistics, 42, 1--17.
