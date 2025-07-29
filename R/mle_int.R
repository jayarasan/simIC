#' Interval-Censored Maximum Likelihood Estimation
#'
#' Estimates distribution parameters by maximizing the interval-censored likelihood.
#'
#' @param left Left bounds of censoring intervals
#' @param right Right bounds of censoring intervals
#' @param dist Distribution name (e.g. "weibull", "loglogistic", "EMV")
#' @return A list containing estimates, standard errors, log-likelihood, and convergence status
#' @examples
#' # Simulate data from a log-logistic distribution
#' set.seed(123)
#' data <- simIC(n = 100, dist = "loglogistic", shape = 1.5, scale = 5, width = 2,
#'               study_start = 1, study_end = 8, uncensored_tol = 0.1)
#' # Fit the model
#' fit <- mle_int(left = data$left, right = data$right, dist = "loglogistic")
#' print(fit$estimates)
#' print(fit$logLik)
#' print(fit$converged)
#' @export
#' @importFrom stats dexp dgamma dlnorm dlogis dnorm dweibull median optim optimHess pexp pgamma plnorm plogis pnorm pweibull qexp qgamma qlnorm qlogis qnorm qweibull rexp rgamma rlnorm rlogis rnorm runif rweibull sd setNames

mle_int <- function(left, right, dist) {
  left <- pmax(left, 1e-5)
  right <- pmax(right, 1e-5)

  loglik <- function(params) {
    logL <- switch(
      dist,
      "weibull" = {
        shape <- params[1]; scale <- params[2]
        mapply(function(l, r) {
          if (is.infinite(r)) log(pmax(1 - pweibull(l, shape, scale), 1e-10))
          else if (l == r) dweibull(l, shape, scale, log = TRUE)
          else log(pmax(pweibull(r, shape, scale) - pweibull(l, shape, scale), 1e-10))
        }, left, right)
      },
      "exp" = {
        rate <- 1 / params[1]
        mapply(function(l, r) {
          if (is.infinite(r)) log(pmax(1 - pexp(l, rate), 1e-10))
          else if (l == r) dexp(l, rate, log = TRUE)
          else log(pmax(pexp(r, rate) - pexp(l, rate), 1e-10))
        }, left, right)
      },
      "logistic" = {
        location <- params[1]; scale <- params[2]
        mapply(function(l, r) {
          if (l <= 0 || r <= 0) return(-1e10)
          if (is.infinite(r)) log(pmax(1 - plogis(l, location, scale), 1e-10))
          else if (l == r) dlogis(l, location, scale, log = TRUE)
          else log(pmax(plogis(r, location, scale) - plogis(l, location, scale), 1e-10))
        }, left, right)
      },
      "loglogistic" = {
        shape <- params[1]; scale <- params[2]
        pllogis <- function(x) 1 / (1 + (scale / x)^shape)
        dllogis <- function(x) {
          top <- (shape / scale) * (x / scale)^(shape - 1)
          bot <- (1 + (x / scale)^shape)^2
          top / bot
        }
        mapply(function(l, r) {
          if (is.infinite(r)) log(pmax(1 - pllogis(l), 1e-10))
          else if (l == r) log(pmax(dllogis(l), 1e-10))
          else log(pmax(pllogis(r) - pllogis(l), 1e-10))
        }, left, right)
      },
      "normal" = {
        location <- params[1]; scale <- params[2]
        mapply(function(l, r) {
          if (is.infinite(r)) log(pmax(1 - pnorm(l, location, scale), 1e-10))
          else if (l == r) dnorm(l, location, scale, log = TRUE)
          else log(pmax(pnorm(r, location, scale) - pnorm(l, location, scale), 1e-10))
        }, left, right)
      },
      "EMV" = {
        location <- params[1]; scale <- params[2]
        z <- function(t) (log(t) - location) / scale
        pemv <- function(t) 1 - exp(-exp(z(t)))
        demv <- function(t) (1 / (t * scale)) * exp(z(t)) * exp(-exp(z(t)))
        mapply(function(l, r) {
          if (l <= 0 || r <= 0) return(-1e10)
          if (is.infinite(r)) log(pmax(1 - pemv(l), 1e-10))
          else if (l == r) log(pmax(demv(l), 1e-10))
          else log(pmax(pemv(r) - pemv(l), 1e-10))
        }, left, right)
      },
      "gamma" = {
        shape <- params[1]; scale <- params[2]
        mapply(function(l, r) {
          if (l <= 0 || r <= 0) return(-1e10)
          if (is.infinite(r)) log(pmax(1 - pgamma(l, shape, scale = scale), 1e-10))
          else if (l == r) dgamma(l, shape, scale = scale, log = TRUE)
          else log(pmax(pgamma(r, shape, scale = scale) - pgamma(l, shape, scale = scale), 1e-10))
        }, left, right)
      },
      "gompertz" = {
        a <- params[1]; b <- params[2]
        pgomp <- function(t) 1 - exp(-(b / a) * (exp(a * t) - 1))
        dgomp <- function(t) b * exp(a * t) * exp(-(b / a) * (exp(a * t) - 1))
        mapply(function(l, r) {
          if (a <= 0 || b <= 0 || l < 0 || r < 0) return(-1e10)
          if (is.infinite(r)) log(pmax(1 - pgomp(l), 1e-10))
          else if (l == r) log(pmax(dgomp(l), 1e-10))
          else log(pmax(pgomp(r) - pgomp(l), 1e-10))
        }, left, right)
      },
      "lognormal" = {
        meanlog <- params[1]; sdlog <- params[2]
        mapply(function(l, r) {
          if (is.infinite(r)) log(pmax(1 - plnorm(l, meanlog, sdlog), 1e-10))
          else if (l == r) dlnorm(l, meanlog, sdlog, log = TRUE)
          else log(pmax(plnorm(r, meanlog, sdlog) - plnorm(l, meanlog, sdlog), 1e-10))
        }, left, right)
      },
      stop("Unsupported distribution.")
    )

    if (any(is.nan(logL) | is.infinite(logL))) return(Inf)
    -sum(logL)
  }

  init <- switch(
    dist,
    "weibull" = c(1, 1),
    "lognormal" = c(1, 1),
    "exp" = c(1),
    c(1, 1)
  )

  lower_bounds <- rep(1e-5, length(init))
  fit <- optim(init, loglik, method = "L-BFGS-B", lower = lower_bounds, hessian = TRUE)

  estimates <- fit$par
  hessian <- fit$hessian

  if (is.null(hessian) || any(is.na(hessian)) || inherits(try(solve(hessian), silent = TRUE), "try-error")) {
    se <- rep(NA, length(estimates))
    warning("Hessian not invertible.")
  } else {
    se <- sqrt(diag(solve(hessian)))
  }

  tvals <- estimates / se
  pvals <- 2 * (1 - pnorm(abs(tvals)))

  # Define parameter names based on distribution
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
                        paste0("param", seq_along(estimates))
  )

  result <- data.frame(
    Estimate = estimates,
    SE = se,
    t = tvals,
    p = formatC(pvals, digits = 6, format = "f")
  )

  rownames(result) <- param_names

  list(
    distribution = dist,
    estimates = result,
    logLik = -fit$value,
    converged = fit$convergence == 0
  )
}

