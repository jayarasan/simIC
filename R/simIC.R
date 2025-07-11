#' Simulate Interval-Censored Survival Data
#'
#' This function simulates interval-censored survival data based on the specified distribution.
#' It supports Weibull, Exponential, Log-Normal, Logistic, and more.
#'
#' @param dist_params Optional list of additional distribution parameters.
#' @param visit_start Optional start time for censoring intervals. Default is 0.
#' @param n Number of samples.
#' @param dist The distribution to use. Options include: "weibull", "exp", "lognormal",
#'   "loglogistic", "logistic", "normal", "EMV", "gamma", "gompertz".
#' @param shape Shape parameter for Weibull, Log-Logistic, etc.
#' @param scale Scale parameter for Weibull, Log-Logistic, etc.
#' @param meanlog Mean log for Log-Normal distribution.
#' @param sdlog Standard deviation for Log-Normal distribution.
#' @param location Location parameter for Logistic, EMV, etc.
#' @param width Width of the interval for censoring.
#'
#' @return A data frame with the following columns: \code{id}, \code{left}, \code{right}, \code{event}, and \code{true_time}.
#'
#' @examples
#' simIC(n = 15, dist = "weibull", shape = 1.2, scale = 5, width = 4)
#' simIC(n = 10, dist = "lognormal", meanlog = 3, sdlog = 1, width = 5)
#'
#' @export
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
      a <- shape  # shape parameter > 0
      b <- scale  # scale parameter > 0

      u <- runif(n)

      # Make sure parameters are positive to avoid invalid values
      if (a <= 0 || b <= 0) stop("Shape and scale must be positive for Gompertz")

      t <- (1 / a) * log(1 - (a / b) * log(u))
      return(t)
    },
    stop("Unsupported distribution.")
  )

  true_times <- get_time(n)
  visit_times <- seq(visit_start, max(true_times + width), by = width)

  left <- right <- rep(NA, n)

  for (i in 1:n) {
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
    id = 1:n,
    left = left,
    right = right,
    event = ifelse(is.infinite(right), 0, 1),
    true_time = true_times
  )
}
