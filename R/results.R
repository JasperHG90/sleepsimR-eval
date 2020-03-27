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
  rf <- jsonlite::read_json(path)
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
# TODO: add label switch postprocess
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
         PD_subj = postprocess_subject_specific(x[[var]], m),
         emiss_mu_bar = postprocess_param_est(x[[var]],m),
         gamma_prob_bar = postprocess_gamma_int(x[[var]],m),
         emiss_var_bar = postprocess_param_est(x[[var]],m),
         emiss_varmu_bar = postprocess_param_est(x[[var]],m),
         credible_intervals = postprocess_ci(x[[var]],m),
         state_order = postprocess_order(x[[var]],m))
  return(out)
}
#' @export
get.sleepsimR_results <- function(x, var = c('uid','scenario_uid','iteration_uid', 'PD_subj','emiss_mu_bar','gamma_prob_bar',
                                             'emiss_var_bar','emiss_varmu_bar','credible_intervals','label_switch',
                                             'state_order'),
                                  type = c("list", "data_frame")) {
  # Match arg
  var <- match.arg(var)
  type <- match.arg(type)
  # Subset
  out <- lapply(x, function(y) get(y, var=var))
  if(type == "list") {
    return(out)
  } else {
    return(do.call(rbind.data.frame, out))
  }
}

#' Postprocessing utility function for parameter estimates
#'
#' @param z
#' @param m integer. Number of hidden states
#'
#' @return
postprocess_param_est <- function(z, m) {
  # Create names
  nams <- paste0("state", 1:m)
  # n dep
  n_dep <- length(z)
  # Add to each
  for(idx in seq_along(z)) {
    names(z[[idx]]$mean) <- nams
    names(z[[idx]]$median) <- nams
    names(z[[idx]]$SE) <- nams
  }
  # To data frame
  df <- as.data.frame(z)
  # Replace periods by underscores
  colnames(df) <- gsub("\\.", "_", colnames(df))
  # Return
  return(df)
}

#' Postprocess utility function for gamma_prob_bar
#'
#' @param z
#' @param m
#'
#' @return
postprocess_gamma_int <- function(z, m) {
  # Number of values is equal to m x (m-1)
  # Make names
  nams <- c()
  for(idx_col in 1:m) {
    for(idx_row in 1:m) {
      nams <- c(nams, paste0("int_S",idx_col, "toS", idx_row))
    }
  }
  # Subset mean
  smean <- z$mean
  # Ignore SE --> names
  names(smean) <- nams
  # data frame and return
  return(data.frame(smean))
}

#' Postprocess credible intervals
#'
#' @param z
#' @param m
#'
#' @return
postprocess_ci <- function(z, m) {
  # Create names
  out_mp <- vector("list", length(z))
  for(lst_idx in seq_along(z)) {
    tmp <- z[[lst_idx]]
    nm <- names(z)[lst_idx]
    if(nm == "gamma_prob_bar") {
      # Make names
      nams <- c()
      for(idx_col in 1:m) {
        for(idx_row in 2:m) {
          for(rngnm in c("lower", "upper")) {
            nams <- c(nams, paste0("gamma_prob_bar_S",idx_col, "toS", idx_row, "_", rngnm))
          }
        }
      }
      # Add names
      names(tmp) <- nams
      # To data frame
      out_mp[[lst_idx]] <- as.data.frame(tmp)
    } else {
      for(var_idx in seq_along(tmp)) {
        nms <- c()
        for(ele_idx in 1:(length(tmp[[var_idx]])/2)) {
          nms <- c(nms, c(paste0("state", ele_idx,"_lower"), paste0("state", ele_idx,"_upper")))
        }
        names(tmp[[var_idx]]) <- nms
      }
      tmpdf <- as.data.frame(tmp)
      colnames(tmpdf) <- gsub("\\.", "_", colnames(tmpdf))
      out_mp[[lst_idx]] <-tmpdf
    }
  }
  # Cbind
  return(
    do.call(cbind.data.frame, out_mp)
  )
}

#' Postprocess state order
#'
#' @param z
#' @param m
#'
#' @return
postprocess_order <- function(z, m) {
  tmp_out <- vector("list", length(z))
  for(idx in seq_along(z)) {
    names(z[[idx]]) <- paste0("state", 1:m)
  }
  # Return data frame
  tmpdf <- as.data.frame(z)
  colnames(tmpdf) <- gsub("\\.", "_", colnames(tmpdf))
  return(
    tmpdf
  )
}

#' Postprocess subject-specific metrics
#'
#' @param z
#' @param m
#'
#' @return
postprocess_subject_specific <- function(z, m) {
  # Number of subjects
  subjs <- paste0("subject_", 1:length(z))
  # Depvar/state names
  depvars <- c("EEG_mean_beta", "EOG_median_theta", "EOG_min_beta")
  # States
  states <- c("state_1", "state_2", "state_3")
  # Expand grid
  grid <- expand.grid(states, depvars)
  # Paste
  statenames <- paste0(grid$Var2, "_", grid$Var1)
  # For each subject, make data frame
  subj_out <- vector("list", length(z))
  names(subj_out) <- subjs
  for(idx in seq_along(subjs)) {
    subj_est <- z[[idx]]
    for(est_idx in seq_along(subj_est)) {
      names(subj_est[[est_idx]]) <- statenames
    }
    # Bind
    subj_est <- do.call(cbind.data.frame, subj_est)
    # Replace periods
    colnames(subj_est) <- gsub("\\.", "_", colnames(subj_est))
    # Add to results
    subj_out[[idx]] <- subj_est
  }
  # Bind
  b <- do.call(rbind.data.frame, subj_out)
  # Unname rows
  row.names(b) <- 1:nrow(b)
  # Return
  return(b)
}
