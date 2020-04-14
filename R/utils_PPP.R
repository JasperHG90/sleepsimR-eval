## Utility functions for PPP

## Functions used to simulate data are placed here

#' Simulate a sleep dataset used in posterior predictive checks.
#'
#' @param n int. number of subjects.
#' @param n_t int. number of observations for each subject.
#' @param gamma vector of length m times m containing transition probabilities.
#' @param mu vector of length m times m with emission distribution means.
#' @param eps vector of length m times m with residual variances.
#' @param zeta vector of floats. between-subject variance of the emission distributions.
#' @param Q float. between-subject variance of the transition probability matrix.
#' @param seed int. random seed used to generate the dataset.
#'
#' @return list containing states and emission distribution observed data. See \link[mHMMbayes]{sim_mHMM}.
#'
#' @importFrom mHMMbayes sim_mHMM
#'
#' @export
simulate_dataset_PPP <- function(n, n_t, gamma, mu, eps, zeta, Q, seed) {
  # Retrieve parameter values from options
  # (set in zzz.R)
  m <- sqrt(length(gamma))
  n_dep <- sqrt(length(mu))
  gamma <- matrix(gamma, nrow = m, ncol= m, byrow=TRUE)
  emiss <- getOption("sleepsimR_simulate")[["emission_bar"]]
  # Set means/variances
  emiss[[1]][,1] <- mu[1:3]
  emiss[[1]][,2] <- eps[1:3]
  emiss[[2]][,1] <- mu[4:6]
  emiss[[2]][,2] <- eps[4:6]
  emiss[[3]][,1] <- mu[7:9]
  emiss[[3]][,2] <- eps[7:9]
  # Set seed
  set.seed(seed)
  # Simulate dataset
  data_simulated <- mHMMbayes::sim_mHMM(
    # Number of observations for each person
    n_t = n_t,
    # Number of persons
    n = n,
    # Type of emission distributions
    data_distr = "continuous",
    # Number of states
    m = m,
    # Number of emission distributions
    n_dep = n_dep,
    # Start state (Awake)
    start_state = 1,
    # Transition probabilities
    gamma = gamma,
    # Emission distribution means + var
    emiss_distr = emiss,
    # Between-subject variance for TPM
    var_gamma = Q,
    # Between-subject variance for emission distributions
    var_emiss = zeta
  )
  # Return
  return(data_simulated)
}

# Sample from posterior
library(sleepsimR)
# Function to draw from posterior and simulate new data
simulate_new_data <- function(posterior) {
  # Fixed quantities
  n <- 41
  n_t <- 1440
  seed <- sample.int(999999999, 1)
  # Draw of transition probabilities
  gammas <- unname(unlist(posterior$gamma_prob_bar[sample.int(nrow(posterior$gamma_prob_bar), 1),]))
  # Draw of emission means and residual variances
  emiss_var <- unname(unlist(posterior$emiss_var_bar[sample.int(nrow(posterior$emiss_var_bar), 1),]))
  emiss_varmu <- unname(unlist(posterior$emiss_varmu_bar[sample.int(nrow(posterior$emiss_varmu_bar), 1),]))
  emiss_varmu <- c(mean(emiss_varmu[1:3]),mean(emiss_varmu[4:6]),mean(emiss_varmu[7:9]))
  emiss_mu <- unname(unlist(posterior$emiss_mu_bar[sample.int(nrow(posterior$emiss_mu_bar), 1),]))
  # New data
  data_new <- simulate_dataset(n, n_t, gammas, emiss_mu, emiss_var, emiss_varmu, 0.1, seed)
  # Combine datasets
  data_out <- as.data.frame(data_new$obs)
  data_out$states <- data_new$states[,2]
  data_return <- list(
    "data" = data_out,
    "seed" = seed
  )
  # Return
  return(data_return)
}

# Posterior predictive checks
PPP_mean_var <- function(newdata) {
  sbv <- newdata %>%
    gather(., var, val, -subj, -states) %>%
    group_by(var, states) %>%
    summarize(avg = mean(val),
              vvar = var(val))
  return(list(sbv$avg, sbv$vvar))
}

# Compute TPM from data
compute_tpm <- function(input_data) {
  m <- matrix(0L, nrow=3, ncol=3)
  state_prev <- 0
  for(idx in seq_along(1:length(input_data))) {
    if(idx == 1) {
      state_prev <- input_data[idx]
      next
    }
    state_now <- input_data[idx]
    m[state_prev, state_now] = m[state_prev, state_now] + 1
    state_prev <- state_now
  }
  # Return
  return(m)
}

# Compute average state length
compute_state_length <- function(input_data) {
  state_lengths <- list(
    "1" = c(),
    "2" = c(),
    "3" = c()
  )
  prev_state <- 0
  current_seq <- c()
  for(idx in seq_along(1:length(input_data))) {
    if(idx == 1) {
      current_seq <- c(current_seq, input_data[idx])
      prev_state <- input_data[idx]
    } else {
      if(input_data[idx] != prev_state) {
        state_lengths[[as.character(prev_state)]] <- c(state_lengths[[as.character(prev_state)]], length(current_seq))
        current_seq <- c(input_data[idx])
        prev_state <- input_data[idx]
      } else {
        current_seq <- c(current_seq, input_data[idx])
      }
    }
  }
  return(state_lengths)
}

# Check on epoch length and transition probabilities
PPP_tpm <- function(newdata) {
  states <- newdata[,c("subj", "states")]
  # Split into subsets of length 720
  from_idx <- seq.int(1, 1440, 480)
  to_idx <- (1440 - from_idx[length(from_idx):1]) + 1
  # Subset
  n_subj <- length(unique(newdata$subj))
  subj_tpms_out <- subj_state_lengths_out <- vector("list", n_subj)
  for(subj_idx in seq_along(1:n_subj)) {
    subj_data <- states %>% filter(subj == subj_idx)
    states_data_split <- states_data_lengths <- vector("list", length(from_idx))
    for(idx in seq_along(states_data_split)) {
      subj_data_curidx <- subj_data[from_idx[idx]:to_idx[idx], "states"]
      subj_data_curidx <- subj_data_curidx[!is.na(subj_data_curidx)]
      states_data_split[[idx]] <- as.vector(compute_tpm(subj_data_curidx))
      states_data_lengths[[idx]] <- compute_state_length(subj_data_curidx)
    }
    subj_tpms_out[[subj_idx]] <- states_data_split
    subj_state_lengths_out[[subj_idx]] <- states_data_lengths
  }
  names(subj_tpms_out) <- names(subj_state_lengths_out) <- 1:n_subj
  return(list(
    "tpms" = subj_tpms_out,
    "lengths" = subj_state_lengths_out
  ))
}
