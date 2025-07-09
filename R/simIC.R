#' Simulate Interval-Censored Survival Data
#'
#' This function simulates interval-censored survival data based on the specified distribution.
#' It supports Weibull, Exponential, Log-Normal, Logistic, and more.
#'
#' @param dist_params Optional list of additional distribution parameters.
#' @param visit_start Optional start time for censoring intervals. Default is 0.
#' @param n Number of samples.
#' @param dist The distribution to use. Options include: "weibull", "exp", "lognormal",
#'   "loglogistic", "logistic", "normal", "gumbel", "gamma", "gompertz".
#' @param shape Shape parameter for Weibull, Log-Logistic, etc.
#' @param scale Scale parameter for Weibull, Log-Logistic, etc.
#' @param meanlog Mean log for Log-Normal distribution.
#' @param sdlog Standard deviation for Log-Normal distribution.
#' @param location Location parameter for Logistic, Gumbel, etc.
#' @param interval_width Width of the interval for censoring.
#'
#' @return A data frame with the following columns: \code{id}, \code{L}, \code{R}, \code{status}, and \code{true_time}.
#'
#' @examples
#' simIC(n = 15, dist = "weibull", shape = 1.2, scale = 5, interval_width = 4)
#' simIC(n = 10, dist = "lognormal", meanlog = 3, sdlog = 1, interval_width = 5)
#'
#' @export
#' @importFrom stats rexp rgamma rlnorm rlogis rnorm runif rweibull


simIC<- function(n = 100,
                 dist = "weibull",
                 shape = 2,
                 scale = 1,
                 meanlog = 0,
                 sdlog = 1,
                 location = 0,
                 dist_params = list(),
                 interval_width = 1,
                 visit_start = 0) {
  if (interval_width <= 0) stop("interval_width must be a positive number")
  # Simulate true event times based on the distribution specified
  get_time <- switch(
    dist,
    "weibull" = function(n) rweibull(n, shape = shape, scale = scale),
    "exp" = function(n) rexp(n, rate = 1/scale),
    "loglogistic" = function(n) {
      u <- runif(n)
      q <- scale * (u / (1 - u))^(1 / shape)
      q
    },
    "lognormal" = function(n) rlnorm(n, meanlog = meanlog, sdlog = sdlog),
    "logistic" = function(n) rlogis(n, location = location, scale = scale),
    "normal" = function(n) rnorm(n, mean = location, sd = scale),
    "gumbel" = function(n) {
      u <- runif(n)
      location - scale * log(-log(u))
    },
    "gamma" = function(n) rgamma(n, shape = shape, scale = scale),
    "gompertz" = function(n) {
      u <- runif(n)
      location + scale * log(-log(u))  # Gompertz approximation
    },
    stop("Unsupported distribution.")
  )

  true_times <- get_time(n)

  # Adjust visit_times based on interval_width
  visit_times <- seq(0, max(true_times + interval_width), by = interval_width)

  # Assign intervals based on the visit times and true event times
  left <- right <- rep(NA, n)

  for (i in 1:n) {
    t_i <- true_times[i]

    # Find the closest visit time intervals
    idx <- which(visit_times >= t_i)

    # If no visit time is greater than the event time, assign right-censoring
    if (length(idx) == 0) {
      left[i] <- max(visit_times) - interval_width
      right[i] <- max(visit_times)
    } else if (idx[1] == 1) { # If the first visit time is the smallest, the event happens before the first visit
      left[i] <- 0
      right[i] <- visit_times[1]
    } else {
      left[i] <- visit_times[idx[1] - 1]
      right[i] <- visit_times[idx[1]]
    }
  }

  # Return a data.frame with the simulated data
  data.frame(id = 1:n,
             left = left,
             right = right,
             event = ifelse(is.infinite(right), 0, 1),
             true_time = true_times)
}
