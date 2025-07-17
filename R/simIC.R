#' Simulate Interval-Censored Survival Data
#'
#' This function simulates interval-censored survival data from various parametric distributions.
#'
#' @param n Number of samples.
#' @param dist The distribution to use. Options: "weibull", "exp", "lognormal", "loglogistic",
#'   "normal", "logistic", "EMV", "gamma", "gompertz".
#' @param shape Shape parameter (for Weibull, Log-Logistic, Gamma, Gompertz).
#' @param scale Scale parameter (for Weibull, Log-Logistic, Gamma, Gompertz, Logistic).
#' @param meanlog Mean log (for Log-Normal).
#' @param sdlog Standard deviation of log (for Log-Normal).
#' @param location Location parameter (for Logistic, Normal, EMV).
#' @param dist_params Optional list of distribution-specific parameters (not used currently).
#' @param width Width of censoring intervals.
#' @param visit_start Starting time for visit schedule (default = 0).
#'
#' @return A data frame with columns: \code{id}, \code{left}, \code{right}, \code{event}, and \code{true_time}.
#' @export
#'
#' @examples
#' simIC(n = 15, dist = "weibull", shape = 1.2, scale = 5, width = 4)
#' simIC(n = 10, dist = "lognormal", meanlog = 3, sdlog = 1, width = 5)
#'
#' @importFrom stats rexp rgamma rlnorm rlogis rnorm runif rweibull
simIC <- function(n = 100,
                  dist = "weibull",
                  shape = 2,
                  scale = 1,
                  meanlog = 0,
                  sdlog = 1,
                  location = 0,
                  dist_params = list(),
                  width = 1,
                  visit_start = 0) {

  if (width <= 0) stop("width must be a positive number")

  get_time <- switch(
    dist,
    "weibull" = function(n) rweibull(n, shape = shape, scale = scale),
    "exp" = function(n) rexp(n, rate = 1 / scale),
    "loglogistic" = function(n) {
      u <- runif(n)
      scale * (u / (1 - u))^(1 / shape)
    },
    "lognormal" = function(n) rlnorm(n, meanlog = meanlog, sdlog = sdlog),
    "logistic" = function(n) rlogis(n, location = location, scale = scale),
    "normal" = function(n) rnorm(n, mean = location, sd = scale),
    "EMV" = function(n) {
      u <- runif(n)
      exp(location + scale * log(-log(1 - u)))
    },
    "gamma" = function(n) rgamma(n, shape = shape, scale = scale),
    "gompertz" = function(n) {
      a <- shape
      b <- scale
      if (a <= 0 || b <= 0) stop("Shape and scale must be positive for Gompertz")
      u <- runif(n)
      (1 / a) * log(1 - (a / b) * log(1 - u))
    },
    stop("Unsupported distribution.")
  )

  true_times <- get_time(n)
  visit_times <- seq(visit_start, max(true_times) + width, by = width)

  left <- right <- numeric(n)

  for (i in seq_len(n)) {
    t_i <- true_times[i]
    idx <- which(visit_times >= t_i)

    if (length(idx) == 0) {
      left[i] <- max(visit_times) - width
      right[i] <- max(visit_times)
    } else if (idx[1] == 1) {
      left[i] <- visit_start
      right[i] <- visit_times[1]
    } else {
      left[i] <- visit_times[idx[1] - 1]
      right[i] <- visit_times[idx[1]]
    }
  }

  data.frame(
    id = seq_len(n),
    left = left,
    right = right,
    event = ifelse(is.infinite(right), 0, 1),
    true_time = true_times
  )
}
