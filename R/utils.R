## Generic utility functions

#' Set session-wise plan for the library future to use
#'
#' This library uses the 'future' library <https://github.com/HenrikBengtsson/future>. The future library was designed to easily work with asynchronous processes in R. Users can specify a 'plan' (e.g. sequential or multisession) that tells R how to process functions.
#'
#' @param plan use either 'sequential' or 'multisession' plan. See <https://github.com/HenrikBengtsson/future> for additional information.
#'
#' @return TRUE.
#'
#' @export
set_future_plan <- function(plan = c("sequential", "multisession")) {
  # Match arg
  plan <- match.arg(plan)
  # Set plan
  options(
    sleepsimReval = list(
      # Future library plan
      "plan" = plan
    )
  )
  return(TRUE)
}
