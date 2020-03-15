# Classes and method for sleepsimR results

#' Parse a sleepsimR result file
#'
#' @param path
#'
#' @return
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
