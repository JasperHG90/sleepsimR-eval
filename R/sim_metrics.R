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
#' @seealso Morris, Tim P., Ian R. White, and Michael J. Crowther. "Using simulation studies to evaluate statistical methods." Statistics in medicine 38.11 (2019): 2074-2102.
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
#' @seealso Morris, Tim P., Ian R. White, and Michael J. Crowther. "Using simulation studies to evaluate statistical methods." Statistics in medicine 38.11 (2019): 2074-2102.
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
#' @seealso Morris, Tim P., Ian R. White, and Michael J. Crowther. "Using simulation studies to evaluate statistical methods." Statistics in medicine 38.11 (2019): 2074-2102.
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
#' @seealso Morris, Tim P., Ian R. White, and Michael J. Crowther. "Using simulation studies to evaluate statistical methods." Statistics in medicine 38.11 (2019): 2074-2102.
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
#' @seealso Morris, Tim P., Ian R. White, and Michael J. Crowther. "Using simulation studies to evaluate statistical methods." Statistics in medicine 38.11 (2019): 2074-2102.
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
#' @param coverage scalar. Coverage of the scenario under investigation.
#' @param nsim scalar. Number of iterations used in the scenario.
#'
#' @seealso Morris, Tim P., Ian R. White, and Michael J. Crowther. "Using simulation studies to evaluate statistical methods." Statistics in medicine 38.11 (2019): 2074-2102.
#'
#' @return MCMC standard error of coverage.
coverage_MCMC_SE <- function(coverage, nsim) {
  sqrt((coverage * (1-coverage))/(nsim))
}

#' Mean Squared Error (MSE)
#'
#' @param simulated_param_values k-length numeric vector. k >= 1 and holds the parameter values of the estimated parameters.
#' @param true_param_value numeric. Value of the ground-truth parameter.
#'
#' @return Mean-Squared Error
#'
#' @seealso Morris, Tim P., Ian R. White, and Michael J. Crowther. "Using simulation studies to evaluate statistical methods." Statistics in medicine 38.11 (2019): 2074-2102.
#'
#' @export
MSE <- function(simulated_param_values, true_param_value) {
  return(
    mean((simulated_param_values - true_param_value)^2)
  )
}

#' MSE MCMC Standard Error
#'
#' @param MSE Mean-Squared Error. \link[sleepsimReval]{MSE}.
#' @param estimates k-length numeric vector. k >= 1 and holds the parameter values of the estimated parameters.
#' @param true_value numeric. Value of the ground-truth parameter.
#' @param nsim scalar. Number of iterations used in the scenario.
#'
#' @seealso Morris, Tim P., Ian R. White, and Michael J. Crowther. "Using simulation studies to evaluate statistical methods." Statistics in medicine 38.11 (2019): 2074-2102.
#'
#' @return MCMC Standard Error of MSE.
MSE_MCMC_SE <- function(MSE, estimates, true_value, nsim) {
  a <- sum(((estimates - true_value)^2 - MSE)^2)
  b <- nsim * (nsim - 1)
  return(sqrt(a/b))
}

#' Wrapper function that computes all simulation metrics
#'
#' This wrapper function is compatible with dplyr and can be called in a dplyr chain
#'  of arguments.
#'
#' @param true_values vector of length n. Contains true population values.
#' @param simulated_values vector of length n. Contains simulated values.
#' @param lower_cci vector of length n. Lower 95\% CCI.
#' @param lower_cci vector of length n. Upper 95\% CCI.
#' @param compute_multimodal boolean. If TRUE, this function will compute a test on the simulated parameter values to check if the distribution of parameter estimates is multimodal. See \link[multimode]{modetest}.
#'
#' @importFrom multimode modetest
#'
#' @return data frame containing:
#' \describe{
#'   \item{bias}{percent bias of the simulation scenario, computed as a percentage relative to the true value.}
#'   \item{bias_mcmc_se}{MCMC standard error of bias estimate, computed as a percentage relative to the bias estimate.}
#'   \item{empirical_se}{empirical standard error computed from the simulated values.}
#'   \item{empirical_se_mcmc_se}{MCMC standard error of the empirical SE.}
#'   \item{MSE}{Mean-Squared Error of the simulated values.}
#'   \item{MSE_mcmc_se}{MCMC standard error of the simulated values.}
#'   \item{coverage}{Coverage given the 95\% CCI.}
#'   \item{coverage_mcmc_se}{MCMC standard error of coverage.}
#'   \item{bias_corr_coverage}{Bias-adjusted coverage. Instead of using the true population value, use the mean of the simulated values. Useful to check whether poor coverage is the result of bias. See 'See also' for reference.}
#'   \item{bias_corr_coverage_mcmc_se}{MCMC standard error of bias-adjusted coverage.}
#'   \item{multimodal}{p-value of the \link[multimode]{modetest} used to check for multimodal distributions.}
#' }
#'
#' @seealso Morris, Tim P., Ian R. White, and Michael J. Crowther. "Using simulation studies to evaluate statistical methods." Statistics in medicine 38.11 (2019): 2074-2102.
#'
#' @examples
#' \dontrun{
#' tst <- data %>%
#'    group_by(scenario_id) %>%
#'    # This is how you use the function
#'    do(summarize_simulation_metrics(.$true_values, .$simulated_values
#'                                    .$lower_cci, .$upper_cci, FALSE))
#' }
#'
#' @export
summarize_simulation_scenario <- function(true_values, simulated_values,
                                          lower_cci, upper_cci,
                                          compute_multimodal = FALSE) {
  out_bias <- unname(bias(true_values[1], simulated_values))
  out_MSE <- unname(MSE(simulated_values, true_values[1]))
  out_ESE <- unname(emperical_SE(simulated_values))
  # Compute coverage
  cci <- map2(lower_cci, upper_cci, function(x,y) c(x, y))
  out_coverage <- unname(coverage(cci, true_values[1]))
  # Compute bias-corrected coverage
  out_coverage_bc <- unname(coverage(cci, mean(simulated_values)))
  df <- data.frame(
    "bias" = (out_bias[1] / true_values[1]) * 100,
    "bias_mcmc_se" = (out_bias[2] / out_bias[1]) * 100,
    "empirical_se" = out_ESE[1],
    "empirical_se_mcmc_se" = out_ESE[2],
    "MSE" = out_MSE,
    "MSE_mcmc_se" = MSE_MCMC_SE(out_MSE,
                                simulated_values,
                                true_values[1],
                                length(simulated_values)),
    "coverage" = out_coverage[1],
    "coverage_mcmc_se" = out_coverage[2],
    "bias_corr_coverage" = out_coverage_bc[1],
    "bias_corr_coverage_mcmc_se" = out_coverage_bc[2]
  )
  if(compute_multimodal) {
    df$multimodal <- modetest(simulated_values, mod0 = 1)$p.value
  }
  return(df)
}

