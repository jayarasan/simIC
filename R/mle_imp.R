#' Imputation-Based MLE for Censored Data
#'
#' Estimates distribution parameters using imputed event times.
#'
#' @param left Left bounds of censoring intervals
#' @param right Right bounds of censoring intervals
#' @param dist Distribution name (e.g. "weibull", "loglogistic", "EMV")
#' @param impute Imputation method: "midpoint", "random", "median",
#' "harmonic_median", "geometric_median", "random_survival"
#' @return A list containing estimates, standard errors, and log-likelihood
#' @export
mle_imp <- function(left, right, dist = "weibull",
                    impute = c("midpoint", "random", "median",
                               "harmonic_median", "geometric_median", "random_survival")) {
  impute <- match.arg(impute)
  
  prelim_times <- ifelse(is.finite(right), (left + right) / 2, left)
  prelim_times <- pmax(prelim_times, 1e-5)
  
  neg_loglik_init <- switch(dist,
                            "weibull" = function(par) -sum(dweibull(prelim_times, par[1], par[2], log = TRUE)),
                            "lognormal" = function(par) -sum(dlnorm(prelim_times, par[1], par[2], log = TRUE)),
                            "loglogistic" = function(par) {
                              shape <- par[1]; scale <- par[2]
                              logf <- log(shape) - log(scale) + (shape - 1) * log(prelim_times / scale) -
                                2 * log1p((prelim_times / scale)^shape)
                              -sum(logf)
                            },
                            "logistic" = function(par) -sum(dlogis(prelim_times, par[1], par[2], log = TRUE)),
                            "normal" = function(par) -sum(dnorm(prelim_times, par[1], par[2], log = TRUE)),
                            "gamma" = function(par) {
                              if (any(par <= 0)) return(Inf)
                              -sum(dgamma(prelim_times, par[1], par[2], log = TRUE))
                            },
                            "gompertz" = function(par) {
                              a <- par[1]; b <- par[2]
                              if (a <= 0 || b <= 0) return(Inf)
                              logf <- log(b) + a * prelim_times - (b / a) * (exp(a * prelim_times) - 1)
                              -sum(logf)
                            },
                            "EMV" = function(par) {
                              loc <- par[1]; scale <- par[2]
                              if (scale <= 0 || any(prelim_times <= 0)) return(Inf)
                              z <- (log(prelim_times) - loc) / scale
                              logf <- -log(prelim_times) - log(scale) - z - exp(-z)
                              -sum(logf)
                            },
                            "exp" = function(par) {
                              if (par[1] <= 0) return(Inf)
                              -sum(dexp(prelim_times, rate = 1 / par[1], log = TRUE))
                            },
                            stop("Unsupported distribution.")
  )
  
  init <- switch(dist,
                 "weibull" = c(1, 1),
                 "lognormal" = c(mean(log(prelim_times)), sd(log(prelim_times))),
                 "loglogistic" = c(1, 1),
                 "logistic" = c(median(prelim_times), sd(prelim_times)),
                 "normal" = c(mean(prelim_times), sd(prelim_times)),
                 "gamma" = c(1, 1),
                 "gompertz" = c(0.1, 0.1),
                 "EMV" = c(mean(log(prelim_times)), sd(log(prelim_times))),
                 "exp" = c(mean(prelim_times))
  )
  
  fit0 <- optim(par = init, fn = neg_loglik_init, method = "L-BFGS-B", lower = rep(1e-5, length(init)))
  est_par <- fit0$par
  
  # Distribution CDF and Quantile functions for survival-based imputation
  cdf_fn <- switch(dist,
                   "weibull" = function(t) pweibull(t, shape = est_par[1], scale = est_par[2]),
                   "lognormal" = function(t) plnorm(t, meanlog = est_par[1], sdlog = est_par[2]),
                   "loglogistic" = function(t) 1 / (1 + (t / est_par[2])^(-est_par[1])),
                   "logistic" = function(t) plogis(t, location = est_par[1], scale = est_par[2]),
                   "normal" = function(t) pnorm(t, mean = est_par[1], sd = est_par[2]),
                   "gamma" = function(t) pgamma(t, shape = est_par[1], scale = est_par[2]),
                   "gompertz" = function(t) 1 - exp(-(est_par[2] / est_par[1]) * (exp(est_par[1] * t) - 1)),
                   "EMV" = function(t) 1 - exp(-exp(-(log(t) - est_par[1]) / est_par[2])),
                   "exp" = function(t) pexp(t, rate = 1 / est_par[1])
  )
  
  quantile_fn <- switch(dist,
                        "weibull" = function(p) qweibull(p, shape = est_par[1], scale = est_par[2]),
                        "lognormal" = function(p) qlnorm(p, meanlog = est_par[1], sdlog = est_par[2]),
                        "loglogistic" = function(p) est_par[2] * ((1 - p) / p)^(-1 / est_par[1]),
                        "logistic" = function(p) qlogis(p, location = est_par[1], scale = est_par[2]),
                        "normal" = function(p) qnorm(p, mean = est_par[1], sd = est_par[2]),
                        "gamma" = function(p) qgamma(p, shape = est_par[1], scale = est_par[2]),
                        "gompertz" = function(p) {
                          val <- tryCatch(log(1 - log(1 - p) * (est_par[1] / est_par[2])) / est_par[1],
                                          error = function(e) NA)
                          ifelse(is.finite(val) & val > 0, val, NA)
                        },
                        "EMV" = function(p) exp(est_par[1] - est_par[2] * log(-log(1 - p)) ),
                        "exp" = function(p) qexp(p, rate = 1 / est_par[1])
  )
  
  times <- numeric(length(left))
  for (i in seq_along(left)) {
    if (!is.finite(right[i])) {
      times[i] <- left[i]
    } else if (impute == "midpoint") {
      times[i] <- (left[i] + right[i]) / 2
    } else if (impute == "random") {
      times[i] <- runif(1, left[i], right[i])
    } else if (impute == "median") {
      times[i] <- median(runif(100, left[i], right[i]))
    } else if (impute == "harmonic_median") {
      x <- runif(100, left[i], right[i]); times[i] <- 1 / median(1 / x)
    } else if (impute == "geometric_median") {
      x <- runif(100, left[i], right[i]); times[i] <- exp(median(log(x)))
    } else if (impute == "random_survival") {
      S_L <- 1 - cdf_fn(left[i])
      S_R <- 1 - cdf_fn(right[i])
      if (S_R >= S_L || S_R < 0 || S_L > 1) {
        times[i] <- (left[i] + right[i]) / 2
      } else {
        u_vals <- runif(100, min = S_R, max = S_L)
        t_vals <- suppressWarnings(quantile_fn(1 - u_vals))
        t_vals <- t_vals[is.finite(t_vals) & t_vals > 0]
        times[i] <- if (length(t_vals) == 0) (left[i] + right[i]) / 2 else sample(t_vals, 1)
      }
    }
  }
  
  times <- pmax(times, 1e-5)
  
  neg_loglik_final <- switch(dist,
                             "weibull" = function(par) -sum(dweibull(times, par[1], par[2], log = TRUE)),
                             "lognormal" = function(par) -sum(dlnorm(times, par[1], par[2], log = TRUE)),
                             "loglogistic" = function(par) {
                               shape <- par[1]; scale <- par[2]
                               logf <- log(shape) - log(scale) + (shape - 1) * log(times / scale) -
                                 2 * log1p((times / scale)^shape)
                               -sum(logf)
                             },
                             "logistic" = function(par) -sum(dlogis(times, par[1], par[2], log = TRUE)),
                             "normal" = function(par) -sum(dnorm(times, par[1], par[2], log = TRUE)),
                             "gamma" = function(par) {
                               if (any(par <= 0)) return(Inf)
                               -sum(dgamma(times, par[1], par[2], log = TRUE))
                             },
                             "gompertz" = function(par) {
                               a <- par[1]; b <- par[2]
                               if (a <= 0 || b <= 0) return(Inf)
                               logf <- log(b) + a * times - (b / a) * (exp(a * times) - 1)
                               -sum(logf)
                             },
                             "EMV" = function(par) {
                               loc <- par[1]; scale <- par[2]
                               if (scale <= 0 || any(times <= 0)) return(Inf)
                               z <- (log(times) - loc) / scale
                               logf <- -log(times) - log(scale) - z - exp(-z)
                               -sum(logf)
                             },
                             "exp" = function(par) {
                               if (par[1] <= 0) return(Inf)
                               -sum(dexp(times, rate = 1 / par[1], log = TRUE))
                             }
  )
  
  fit <- optim(par = est_par, fn = neg_loglik_final, method = "L-BFGS-B",
               lower = rep(1e-5, length(est_par)), hessian = FALSE)
  hessian <- tryCatch(optimHess(fit$par, neg_loglik_final), error = function(e) NULL)
  
  param_names <- switch(dist,
                        "weibull"     = c("shape", "scale"),
                        "loglogistic" = c("shape", "scale"),
                        "gamma"       = c("shape", "scale"),
                        "gompertz"    = c("shape", "scale"),
                        "exp"         = "scale",
                        "lognormal"   = c("meanlog", "sdlog"),
                        "logistic"    = c("location", "scale"),
                        "normal"      = c("location", "scale"),
                        "EMV"         = c("location", "scale"),
                        c("param1", "param2")
  )
  param_names <- param_names[seq_along(fit$par)]
  estimates <- setNames(fit$par, param_names)
  
  if (fit$convergence == 0 && !is.null(hessian)) {
    vcov <- tryCatch(solve(hessian), error = function(e) matrix(NA, length(estimates), length(estimates)))
    se <- sqrt(diag(vcov))
    tvals <- estimates / se
    pvals_raw <- 2 * (1 - pnorm(abs(tvals)))
    pvals <- ifelse(pvals_raw < 2e-16, "<2e-16", formatC(pvals_raw, digits = 6, format = "f"))
    results <- data.frame(Estimate = estimates, Std.Error = se, t.value = tvals, p.value = pvals)
  } else {
    results <- data.frame(Estimate = estimates, Std.Error = NA, t.value = NA, p.value = NA)
  }
  
  return(list(
    distribution = dist,
    imputation = impute,
    estimates = results,
    logLik = -fit$value,
    converged = fit$convergence == 0
  ))
}
