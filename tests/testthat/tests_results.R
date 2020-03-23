context("Test result methods and class")
f <- "tests/testthat/results"

test_that("Can parse result files", {
  # Parse results
  r <- parse_sleepsimR_results(f)
  # Check class
  expect_s3_class(r, "sleepsimR_results")
  # Length
  expect_length(r, 2)
  expect_length(r[[1]], 11)
  expect_named(r[[1]], c('uid','scenario_uid','iteration_uid','PD_subj','emiss_mu_bar',
                         'gamma_int_bar','emiss_var_bar','emiss_varmu_bar','credible_intervals',
                         'label_switch','state_order'))
})

test_that("Can get specific result fields", {
  # Parse results
  r <- parse_sleepsimR_results(f)
  # Get fields
  for(v in c('uid','scenario_uid','iteration_uid','PD_subj','emiss_mu_bar',
               'gamma_int_bar','emiss_var_bar','emiss_varmu_bar','credible_intervals', 'state_order')) {
    g <- get(r, v)
    if(v %in% c("uid", "scenario_uid", "iteration_uid")) {
      expect_length(g, 2)
      expect_type(g[[1]], "character")
      expect_named(g)
    } else if(v %in% c("PD_subj", "emiss_mu_bar", "gamma_int_bar", "emiss_var_bar",
                       "emiss_varmu_bar", "credible_intervals", "state_order")) {
      expect_length(g, 2)
      expect_type(g[[1]], "list")
      gc <- do.call(rbind.data.frame, g)
    }
  }
})
