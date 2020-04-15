context("Test result methods and class")
f <- "results/0a2bb17076ba5cec2336e4f6b8d07ad8.json"

test_that("Can parse result files", {
  # Parse results
  r <- parse_sleepsimR_result(f)
  # Check class
  expect_s3_class(r, "sleepsimR_result")
  # Length
  expect_length(r, 11)
  expect_named(r, c('uid','scenario_uid','iteration_uid','PD_subj','emiss_mu_bar',
                     'gamma_prob_bar','emiss_var_bar','emiss_varmu_bar','credible_intervals',
                     'label_switch','state_order'))
})
