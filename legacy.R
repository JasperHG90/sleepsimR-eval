# Legacy

#' Parse multiple sleepsimR result files
#'
#' @param folder_path
#'
#' @return
#'
#' @importFrom jsonlite read_json
#' @importFrom future.apply future_lapply
#' @importFrom future plan
#'
#' @export
parse_sleepsimR_results <- function(folder_path) {
  # Check folder exists
  if(!dir.exists(folder_path)) {
    stop(paste0("Folder '", folder_path, "' does not exist."))
  }
  # Get list of files
  f <- list.files(folder_path)
  # Make filepaths
  fp <- file.path(folder_path, f)
  m <- future::availableCores()
  # Split into chunks of size fp // m
  cs <- ceiling(length(fp) / unname(m))
  csout <- split(fp, ceiling(seq_along(fp)/cs))
  # Read data
  t1 <- Sys.time()
  res <- future_lapply(csout,
                       FUN = function(x) {
                         # For each x, parse result
                         lapply(x, function(y) {
                           parse_sleepsimR_result(y)
                         })
                       },
                       future.packages = c("jsonlite"),
                       future.globals = c("parse_sleepsimR_result"))
  Sys.time() - t1
  # Unlist
  return(res)
  # Set names on list
  names(res) <- vapply(res, function(x) x$uid, "char")
  # Add structure
  class(res) <- "sleepsimR_results"
  # Return
  return(res)
}

#' @export
#'
#' @importFrom future.apply future_lapply
get.sleepsimR_results <- function(x, var = c('uid','scenario_uid','iteration_uid', 'PD_subj','emiss_mu_bar','gamma_prob_bar',
                                             'emiss_var_bar','emiss_varmu_bar','credible_intervals','label_switch',
                                             'state_order'),
                                  type = c("list", "data_frame")) {
  # Match arg
  var <- match.arg(var)
  type <- match.arg(type)
  # Subset
  out <- future_lapply(x, FUN = function(y) get(y, var=var),
                       future.globals = c(
                         "postprocess_subject_specific",
                         "postprocess_param_est",
                         "postprocess_gamma_int",
                         "postprocess_ci",
                         "postprocess_order"
                       ),
                       future.packages = c("sleepsimReval"))
  if(type == "list") {
    return(out)
  } else {
    return(do.call(rbind.data.frame, out))
  }
}

#' Get method
#'
#' @param x an object of type sleepsimR_result
#'
#' @return
#'
#' @export
get <- function(x, ...) {
  UseMethod("get", x)
}
#' @export
get.sleepsimR_result <- function(x, var = c('uid','scenario_uid','iteration_uid', 'PD_subj', 'emiss_mu_bar','gamma_prob_bar',
                                            'emiss_var_bar','emiss_varmu_bar','credible_intervals','state_order')) {
  # Match arg
  var <- match.arg(var)
  # Get number of states
  m <- length(x$emiss_mu_bar[[1]][[1]])
  # Switch. Depends on output type how it should be postprocessed
  out <- switch(var,
                uid = x[[var]],
                scenario_uid = x[[var]],
                iteration_uid = x[[var]],
                PD_subj = sleepsimReval:::postprocess_subject_specific(x[[var]], m),
                emiss_mu_bar = sleepsimReval:::postprocess_param_est(x[[var]],m),
                gamma_prob_bar = sleepsimReval:::postprocess_gamma_int(x[[var]],m),
                emiss_var_bar = sleepsimReval:::postprocess_param_est(x[[var]],m),
                emiss_varmu_bar = sleepsimReval:::postprocess_param_est(x[[var]],m),
                credible_intervals = sleepsimReval:::postprocess_ci(x[[var]],m),
                state_order = sleepsimReval:::postprocess_order(x[[var]],m))
  return(out)
}
