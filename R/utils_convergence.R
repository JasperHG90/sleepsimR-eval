## Utility functions to check for convergence between two MCMC chains

# Trace plots
library(ggplot2)
library(coda)

#' Make trace plots of the convergence of two chains
#'
#' @seealso Explanation in \link[coda]{traceplot}.
#'
#' @param mod1 \link[mHMMbayes]{mHMM_cont} model.
#' @param mod2 \link[mHMMbayes]{mHMM_cont} model.
#' @param param character. name of the variable to plot.
#' @param var numeric. index of the emission distribution.
#' @param thin numeric. thinning interval for the posterior distribution. Defaults to 5.
#'
#' @importFrom coda traceplot
#'
#' @export
tpp <- function(mod1, mod2, param, var, thin = 5) {
  # Combine chains
  chains_comb <- combine_mcmc_chains(mod1, mod2, param = param,
                                     var= var, thin = thin)
  # Plot
  ## Max of three plots in rows
  row_out <- mod1$input$m
  col_out <- ceiling(ncol(chains_comb[[1]]) / mod1$input$m)
  par(mfrow = c(row_out, col_out),     # 2x2 layout
      mar = c(4, 4, 2, 2), # space for one row of text at ticks and to separate plots
      mgp = c(3, 1, 0),    # axis label at 2 rows distance, tick labels at 1 row
      xpd = NA)
  coda::traceplot(chains_comb)
}

#' Make density plots of the convergence of two chains
#'
#' @seealso Explanation in \link[mcmcplots]{denoverplot}.
#'
#' @param mod1 \link[mHMMbayes]{mHMM_cont} model.
#' @param mod2 \link[mHMMbayes]{mHMM_cont} model.
#' @param param character. name of the variable to plot.
#' @param var numeric. index of the emission distribution.
#' @param thin numeric. thinning interval for the posterior distribution. Defaults to 5.
#'
#' @importFrom mcmcplots denoverplot
#'
#' @export
dens_plot <- function(mod1, mod2, param, var, thin = 5) {
  # Combine chains
  chains_comb <- combine_mcmc_chains(mod1, mod2, param = param,
                                     var= var, thin = thin)
  # Plot
  ## Max of three plots in rows
  row_out <- mod1$input$m
  col_out <- ceiling(ncol(chains_comb[[1]]) / mod1$input$m)
  #par(mfrow = c(row_out, col_out),     # 2x2 layout
  #    mar = c(4, 4, 2, 2), # space for one row of text at ticks and to separate plots
  #    mgp = c(3, 1, 0),    # axis label at 2 rows distance, tick labels at 1 row
  #    xpd = NA)
  mcmcplots::denoverplot(chains_comb[[1]], chains_comb[[2]], style="plain")
}

#' Compute the Gelman-Rubin potential scale reduction factor.
#'
#' @seealso Explanation of GRS in \link[coda]{gelman.diag}.
#'
#' @param mod1 \link[mHMMbayes]{mHMM_cont} model.
#' @param mod2 \link[mHMMbayes]{mHMM_cont} model.
#' @param param character. name of the variable to plot.
#' @param var numeric. index of the emission distribution.
#' @param thin numeric. thinning interval for the posterior distribution. Defaults to 5.
#'
#' @importFrom coda gelman.diag
#'
#' @export
compute_grs <- function(mod1, mod2, param, var, thin = 5) {
  # Combine chains
  chains_comb <- combine_mcmc_chains(mod1, mod2, param, var, thin = thin)
  # Compute diag
  coda::gelman.diag(chains_comb, autoburnin = FALSE)
}

#' Make autocorrelation plots of the convergence of two chains
#'
#' @seealso Explanation in \link[coda]{autocorr.plot}.
#'
#' @param mod1 \link[mHMMbayes]{mHMM_cont} model.
#' @param mod2 \link[mHMMbayes]{mHMM_cont} model.
#' @param param character. name of the variable to plot.
#' @param var numeric. index of the emission distribution.
#' @param thin numeric. thinning interval for the posterior distribution. Defaults to 5.
#'
#' @importFrom coda autocorr.plot
#' @importFrom coda as.mcmc
#' @importFrom sleepsimR burn
#'
#' @export
autocorr_plot <- function(mod, param, var, thin = 5) {
  # Subset chain and burn
  # Retrieve posteriors
  modb <- sleepsimR::burn(mod)[[param]]
  if(param != "gamma_int_bar") {
    modb <- modb[[var]]
  }
  # Convert to mcmc
  modmcmc <- as.mcmc(modb)
  # Thin
  modmcmc <- thin_mcmc(modmcmc, thin=thin)
  # Plot
  coda::autocorr.plot(modmcmc)
}

# Helper functions

#' Thin the mcmc chain
#'
#' @importFrom coda mcmc
thin_mcmc <- function(x, thin = 5) {
  dim_chain <- dim(x)
  keep_its <- seq.int(from = 1, dim_chain[1], by = 5)
  # Remove any where index of keep_its > number of rows in mcmc chain
  keep_its <- keep_its[keep_its <= dim_chain[1]]
  # Subset
  return(coda::mcmc(x[keep_its,], thin = thin))
}

#' Combine two separate chains into one chain
#'
#' @importFrom coda mcmc
combine_mcmc_chains <- function(mod1, mod2, param, var, thin = 5) {
  # Retrieve posteriors
  mod1b <- sleepsimR::burn(mod1)[[param]]
  mod2b <- sleepsimR::burn(mod2)[[param]]
  if(!param %in% c("gamma_int_bar", "gamma_prob_bar", "gamma_int_subj", "PD_subj")) {
    mod1b <- mod1b[[var]]
    mod2b <- mod2b[[var]]
  }
  if(param %in% c("gamma_int_subj")) {
    for(idx in seq_along(mod1b)) {
      md1bt <- as.data.frame(mod1b[[idx]])
      md1bt$subj_idx <- idx
      md2bt <- as.data.frame(mod2b[[idx]])
      md2bt$subj_idx <- idx
      # To chain and thin
      md1bt <- mcmc(md1bt, thin=1)
      md2bt <- mcmc(md2bt, thin=1)
      mod1bmcmc <- thin_mcmc(md1bt, thin=thin)
      mod2bmcmc <- thin_mcmc(md2bt, thin=thin)
      # Replace in list
      mod1b[[idx]] <- mod1bmcmc
      mod2b[[idx]] <- mod2bmcmc
    }
    # Bind
    mod1b <- do.call(rbind.data.frame, mod1b)
    mod2b <- do.call(rbind.data.frame, mod2b)
    # To mcmc and return
    return(as.mcmc.list(list(mcmc(mod1b, thin=thin), mcmc(mod2b, thin=thin))))
  } else if(param %in% c("PD_subj")) {
    for(idx in seq_along(mod1b)) {
      md1bt <- as.data.frame(mod1b[[idx]])[,1:9]
      md1bt$subj_idx <- idx
      md2bt <- as.data.frame(mod2b[[idx]])[,1:9]
      md2bt$subj_idx <- idx
      # To chain and thin
      md1bt <- mcmc(md1bt, thin=1)
      md2bt <- mcmc(md2bt, thin=1)
      mod1bmcmc <- thin_mcmc(md1bt, thin=thin)
      mod2bmcmc <- thin_mcmc(md2bt, thin=thin)
      # Replace in list
      mod1b[[idx]] <- mod1bmcmc
      mod2b[[idx]] <- mod2bmcmc
    }
    # Bind
    mod1b <- do.call(rbind.data.frame, mod1b)
    mod2b <- do.call(rbind.data.frame, mod2b)
    # To mcmc and return
    return(as.mcmc.list(list(mcmc(mod1b, thin=thin), mcmc(mod2b, thin=thin))))
  }
  # Convert to mcmc
  mod1bmcmc <- mcmc(mod1b, thin=1)
  mod2bmcmc <- mcmc(mod2b, thin=1)
  # Thin
  mod1bmcmc <- thin_mcmc(mod1bmcmc, thin=thin)
  mod2bmcmc <- thin_mcmc(mod2bmcmc, thin=thin)
  # Return
  return(as.mcmc.list(list(mod1bmcmc, mod2bmcmc)))
}
