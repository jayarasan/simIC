# simIC

`simIC` is an R package for simulating **interval-censored survival data** from various distributions, such as Weibull, exponential, log-normal, gamma, logistic, and more.

## 🔧 Installation

```r
# Install from GitHub
install.packages("remotes")  # if not already installed
remotes::install_github("jayarasan/simIC")


# Example

library(simIC)

# Simulate data from a Weibull distribution
data <- simIC(n = 15, dist = "weibull", shape = 1.2, scale = 5, interval_width = 4)

head(data)

#Features
Simulates interval-censored survival data

Supports common distributions: Weibull, Exponential, Log-normal, Logistic, Gamma, etc.

Easy-to-use interface for teaching, testing, or model evaluation

#License


---

### ✅ After saving:

1. **Save the file** (`Ctrl + S` or click Save)
2. In your console, commit the change:

```r
usethis::use_git_commit("Add custom README")

