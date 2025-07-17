# simIC

<!-- badges: start -->
<!-- badges: end -->

The `simIC` package provides flexible tools for simulating **interval-censored survival data** from a variety of parametric distributions. It is useful for teaching, model development, and testing methods in survival analysis.

## ✨ Features

- Supports commonly used distributions:
  - Weibull
  - Exponential
  - Log-Normal
  - Logistic
  - Normal
  - Log-Logistic
  - Gamma
  - Gompertz
  - EMV (Extreme Minimum Value / Gumbel)
- Generates interval-censored event times based on user-defined visit schedules
- Includes direct and imputation-based MLE functions: `mle_int()` and `mle_imp()`

## 📦 Installation

You can install the development version of `simIC` from GitHub:

```r
# Install from GitHub
install.packages("remotes")
remotes::install_github("jayarasan/simIC")

library(simIC)

# Simulate interval-censored data from a Weibull distribution
data <- simIC(n = 100, dist = "weibull", shape = 1.5, scale = 5, width = 2)

# Fit model using imputation (midpoint)
fit <- mle_imp(data$left, data$right, dist = "weibull", impute = "midpoint")
print(fit$estimates)

