# Functions that deal with allocations file go here

#' Read a sleepsimR allocations file
#'
#' @param path path to the allocations file
#'
#' @return S3 object of type 'sleepsimR_allocations'
#'
#' @importFrom jsonlite read_json
#' @importFrom lubridate as_datetime
#'
#' @details A sleepsimR allocations file is a json file that has been created
#'     by the sleepsimR-api program <https://github.com/JasperHG90/sleepsimR-api>.
#'     It contains two elements: (1) a list of dicts where the key is the container-id
#'     (i.e. the container that processed an iteration of the study) and the values
#'     contain (a) the unique id of the iteration on which the program worked, (b)
#'     the status of the program ('working', 'error' or 'completed') and (c) two timestamp
#'     values (if completed). The first time stamp value indicates when the iteration
#'     began, and the second timestamp value indicates when it was finished. This function
#'     only returns containers which completed their task.
#'
#' @export
read_allocations_file <- function(path) {
  # Check file
  if(!file.exists(path)) {
    stop(paste0("File '", path, "' does not exist."))
  }
  # Read data
  allocs <- jsonlite::read_json(path)
  # Select first element and convert timestamps
  allocsfe <- lapply(allocs$allocations, function(x) {
    if(x$status != "completed") {
      return(NULL)
    }
    x$ts_request <- lubridate::as_datetime(x$ts_request)
    x$ts_finished <- lubridate::as_datetime(x$ts_finished)
    x$status <- NULL
    x$runtime_minutes <- as.numeric(difftime(x$ts_finished, x$ts_request))
    return(x)
  })
  # Filter null
  allocsfe <- allocsfe[vapply(allocsfe, function(x) !is.null(x), TRUE)]
  # Add structure
  class(allocsfe) <- "sleepsimR_allocations"
  # Return
  return(allocsfe)
}

#' Convert allocations to data frame
#'
#' @param x object of type sleepsimR_allocations
#'
#' @return a data frame with 5 columns
#' \describe{
#'     \item{container_id}{Unique id of the container that ran the iteration.}
#'     \item{iteration_id}{Unique id of the iteration. See also the function 'generate_scenarios' in the sleepsimR library <https://github.com/JasperHG90/sleepsimR>}
#'     \item{ts_request}{POSIXct. Date/time stamp indicating when the iteration started.}
#'     \item{ts_finished}{POSIXct. Date/time stamp indicating when the iteration finished.}
#'     \item{runtime_minutes}{Numeric. Time, in minutes, between ts_request and ts_finished.}
#' }
#'
#' @export
as.data.frame.sleepsimR_allocations <- function(x) {
  # Convert list to df
  cont_ids <- names(x)
  x_df <- do.call(rbind.data.frame, unname(x))
  # Add container id
  x_df$container_id <- cont_ids
  # Reshuffle and return
  return(x_df[,c(5,1,2,3,4)])
}

#' Compute the average runtime of the allocations
#'
#' @param x object of type sleepsimR_allocations
#'
#' @return numeric scalar. Average runtime across all iterations.
#'
#' @export
average_runtime <- function(x, ...) {
  UseMethod("average_runtime", x)
}
#' @export
average_runtime.sleepsimR_allocations <- function(x) {
  # Return mean of runtime data
  return(mean(vapply(x, function(x) x$runtime_minutes, 0.0)))
}
