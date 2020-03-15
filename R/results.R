# Classes and method for sleepsimR results

#' Parse a sleepsimR result file
#'
#' @param path
#'
#' @return
#'
#' @importFrom jsonlite read_json
#'
#' @export
parse_sleepsimR_result <- function(path) {
  # Check file
  if(!file.exists(path)) {
    stop(paste0("File '", path, "' does not exist."))
  }
  # Load json
  rf <- jsonlite::read_json(f)
  # Add structure
  class(rf) <- "sleepsimR_result"
  # Return
  return(rf)
}

#' Parse multiple sleepsimR result files
#'
#' @param folder_path
#'
#' @return
#'
#' @importFrom jsonlite read_json
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
  # Read data
  res <- lapply(fp, parse_sleepsimR_result)
  # Set names on list
  names(res) <- vapply(res, function(x) x$uid, "char")
  # Add structure
  class(res) <- "sleepsimR_results"
  # Return
  return(res)
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
get.sleepsimR_result <- function(x, var = c('uid','scenario_uid','iteration_uid','emiss_mu_bar','gamma_int_bar',
                                            'emiss_var_bar','emiss_varmu_bar','credible_intervals','label_switch',
                                            'state_order')) {
  # Match arg
  var <- match.arg(var)
  # Subset
  out <- x[[var]]
  # Switch. Depends on output type how it should be postprocessed
  # ...
  return(out)
}
#' @export
get.sleepsimR_results <- function(x, var = c('uid','scenario_uid','iteration_uid','emiss_mu_bar','gamma_int_bar',
                                             'emiss_var_bar','emiss_varmu_bar','credible_intervals','label_switch',
                                             'state_order'),
                                  type = c("list", "data_frame")) {
  # Match arg
  var <- match.arg(var)
  # Subset
  out <- lapply(x, function(y) get(y, var=var))
  return(out)
}

# Postprocessing utility function for parameter estimates
