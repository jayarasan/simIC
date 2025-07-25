#' Simulate Interval-, Left-, Right-, and Uncensored Survival Data
#'
#' Simulates survival data with optional left-censoring, right-censoring, and uncensoring thresholds.
#'
#' @param n Number of samples.
#' @param dist Distribution name ("weibull", "exp", "lognormal", "loglogistic", "normal", "logistic", "EMV", "gamma", "gompertz").
#' @param shape,scale Distribution parameters for applicable distributions.
#' @param meanlog,sdlog For lognormal.
#' @param location For normal, logistic, and EMV.
#' @param width Visit interval width.
#' @param visit_start First visit time.
#' @param study_start Optional: left-censoring cutoff.
#' @param study_end Optional: right-censoring cutoff.
#' @param uncensored_tol Tolerance to treat (left, right) as exact event.
#' @return A data frame with columns: id, left, right, true_time, censoring
#' @export
simIC <- function(n = 100,
                  dist = "weibull",
                  shape = 2,
                  scale = 1,
                  meanlog = 0,
                  sdlog = 1,
                  location = 0,
                  width = 1,
                  visit_start = 0,
                  study_start = NULL,
                  study_end = NULL,
                  uncensored_tol = 0.1) {
  
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
    "logistic" = function(n) {
      x <- rlogis(n, location = location, scale = scale)
      while (any(x < 0)) {
        x[x < 0] <- rlogis(sum(x < 0), location = location, scale = scale)
      }
      x
    },
    "normal" = function(n) {
      x <- rnorm(n, mean = location, sd = scale)
      x[x <= 0] <- 1e-5
      x
    },
    "EMV" = function(n) {
      u <- runif(n)
      exp(location + scale * log(-log(1 - u)))
    },
    "gamma" = function(n) rgamma(n, shape = shape, scale = scale),
    "gompertz" = function(n) {
      a <- shape; b <- scale
      if (a <= 0 || b <= 0) stop("Shape and scale must be positive for Gompertz")
      u <- runif(n)
      (1 / a) * log(1 - (a / b) * log(1 - u))
    },
    stop("Unsupported distribution.")
  )
  
  true_times <- get_time(n)
  visit_times <- seq(visit_start, max(true_times) + width, by = width)
  
  left <- right <- numeric(n)
  censoring <- character(n)
  
  for (i in seq_len(n)) {
    t_i <- true_times[i]
    
    # Right censoring
    if (!is.null(study_end) && t_i > study_end) {
      left[i] <- study_end
      right[i] <- Inf
      censoring[i] <- "right"
      next
    }
    
    # Left censoring
    if (!is.null(study_start) && t_i < study_start) {
      left[i] <- 0
      right[i] <- study_start
      censoring[i] <- "left"
      next
    }
    
    # Uncensored (if close to any visit time)
    if (!is.null(uncensored_tol)) {
      time_diffs <- abs(visit_times - t_i)
      if (min(time_diffs) <= uncensored_tol) {
        left[i] <- right[i] <- t_i
        censoring[i] <- "uncensored"
        next
      }
    }
    
    # Interval censoring
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
    
    censoring[i] <- "interval"
  }
  
  data.frame(
    id = seq_len(n),
    left = left,
    right = right,
    true_time = true_times,
    censoring = censoring
  )
}
