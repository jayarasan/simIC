# simIC

<!-- badges: start -->
<!-- badges: end -->

The `simIC` package provides tools for simulating and analyzing **interval-censored survival data**, including left-, right-, and uncensored observations, using a variety of parametric distributions. It is useful for teaching, model development, and method evaluation in survival analysis.

## ✨ Features

* Supports commonly used parametric distributions:
  * Weibull  
  * Exponential  
  * Log-Normal  
  * Logistic  
  * Normal  
  * Log-Logistic  
  * Gamma  
  * Gompertz  
  * EMV (Extreme Minimum Value / Gumbel)

* Simulates survival data with **interval**, **left**, **right**, and **uncensored** observations using user-defined visit schedules (`start_time`, `end_time`) and an optional tolerance (`uncensored_tol`) for detecting exact event times.

* Provides two estimation functions:

  * `mle_int()`  
    - Performs **direct maximum likelihood estimation**  
    - Automatically detects and handles:
      - **Interval-censored**: contribution from `F(Ri) - F(Li)`  
      - **Left-censored**: contribution from `F(Ri)`  
      - **Right-censored**: contribution from `1 - F(Li)`  
      - **Uncensored**: contribution from `f(ti)`

  * `mle_imp()`  
    - Uses **imputation-based likelihood** for **interval- and uncensored** data  
    - For **left- and right-censored**, uses the proper **likelihood contributions without imputation**  
    - Specifically:
      - **Interval-censored**: event times imputed from `(Li, Ri)` using midpoint, random, medians, or survival-based methods  
      - **Left-censored**: contribution from `F(Ri)`  
      - **Right-censored**: contribution from `1 - F(Li)`  
      - **Uncensored**: contribution from `f(ti)`

## 📦 Installation

You can install the development version of `simIC` from GitHub:

```r
install.packages("remotes")
remotes::install_github("jayarasan/simIC")
library(simIC)

🧪 Simulate Survival Data
# Interval-censored data only (no visit window)
data <- simIC(n = 100, dist = "weibull", shape = 1.5, scale = 5, width = 2)

# Left-, right-, and uncensored data using a follow-up window and tolerance
data <- simIC(n = 100, dist = "weibull", shape = 1.5, scale = 5,
              width = 2, start_time = 0, end_time = 10, uncensored_tol = 0.1)

📈 Model Fitting Examples

# Direct MLE for interval-censored data
fit_int <- mle_int(data$left, data$right, dist = "weibull")
print(fit_int$estimates)

# Imputation-based MLE (midpoint)
fit_imp <- mle_imp(data$left, data$right, dist = "weibull", impute = "midpoint")
print(fit_imp$estimates)

