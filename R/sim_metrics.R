## Metrics (% bias etc) go here
# See https://cran.r-project.org/web/packages/rsimsum/vignettes/A-introduction.html

#' Compute the parameter bias of a vector estimated parameters versus the ground-truth value of that parameter
#'
#' @param true_param_value numeric. Value of the ground-truth parameter.
#' @param simulated_param_values k-length numeric vector. k >= 1 and holds the parameter values of the estimated parameters.
#'
#' @return numeric vector. Contains two elements: (1) average bias of simulated values versus the ground-truth value and (2) MCMC SE of the bias value
#'
#' @details This function computes the percentage bias by using the signed mean difference.
#'
#' @export
bias <- function(true_param_value, simulated_param_values) {
  bias <- mean(simulated_param_values - true_param_value)
  bmcse <- bias_MCMC_SE(simulated_param_values)
  ret <- c(bias, bmcse)
  names(ret) <- c("bias", "MCMC_SE")
  return(ret)
}

#' MCMC Standard Error of the bias value
#'
#' @param x k-length numeric vector. k >= 1 and holds the parameter values of the estimated parameters.
#'
#' @return Bias MCMC SE
bias_MCMC_SE <- function(x) {
  nsim <- length(x)
  a <- (1 / (nsim - 1))
  b <- sum((x - mean(x))^2)
  return(
    sqrt((a * b)/(nsim))
  )
}

#' Compute the emperical SE
#'
#' @param x numeric vector. Simulated parameter estimates
#'
#' @return numeric vector with two values: (1) Emperical standard error and (2) MCMC SE of this value
#'
#' @export
emperical_SE <- function(x) {
  ESE <- sqrt((1/(length(x) - 1)) * sum((x - mean(x))^2))
  MCMCSE <- emperical_MCMC_SE(ESE, length(x))
  ret <- c(ESE, MCMCSE)
  names(ret) <- c("emperical_se", "MCMC_SE")
  return(ret)
}

#' Compute the MCMC SE
#'
#' @param x numeric vector. Simulated parameter estimates.
#'
#' @return numeric scalar. MCMC standard error.
emperical_MCMC_SE <- function(emp_se, nsim) {
  emp_se / sqrt(2*nsim - 1)
}

#' Compute coverage
#'
#' @param CI list. Each element contains the upper and lower values of the 95\% CI of an interation.
#' @param true_param_value numeric scalar. True value of the parameter.
#'
#' @return numeric vector with two values: (1) Probability that the 95\% CI contains the true value and (2) MCMC Standard Error
#'
#' @export
coverage <- function(CI, true_param_value) {
  cvr <- mean(vapply(CI, function(x) (x[1] <= true_param_value) & (true_param_value <= x[2]), 0))
  cvrSE <- coverage_MCMC_SE(cvr, length(CI))
  ret <- c(cvr, cvrSE)
  names(ret) <- c("coverage", "MCMC_SE")
  return(ret)
}

#' Coverage MCMC Standard Error
#'
#' @param coverage XX
#' @param nsim XX
#'
#' @return XX
coverage_MCMC_SE <- function(coverage, nsim) {
  sqrt((coverage * (1-coverage))/(nsim))
}

#' Mean Squared Error (MSE)
#'
#' @param values XX
#' @param ref XX
#'
#' @return XX
#'
#' @export
MSE <- function(values, ref) {
  return(
    mean((values - ref)^2)
  )
}

#' MSE MCMC Standard Error
#'
#' @param MSE
#' @param estimates
#' @param true_value
#' @param nsim
#'
#' @return YY
MSE_MCMC_SE <- function(MSE, estimates, true_value, nsim) {
  a <- sum(((estimates - true_value)^2 - MSE)^2)
  b <- nsim * (nsim - 1)
  return(sqrt(a/b))
}
