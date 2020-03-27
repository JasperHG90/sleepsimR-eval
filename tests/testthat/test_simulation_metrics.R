context("Test simulation metrics (parameter bias etc)")

test_that("Can compute bias and MCMC error for bias", {
  # True value
  tv <- 5
  # Simulated values
  set.seed(4445896)
  sv <- tv + rnorm(500)
  # Compute
  b <- bias(tv, sv)
  expect_equal(unname(round(b, 2)), c(-0.01, 0.05))
  expect_named(b, c("bias", "MCMC_SE"))
})

test_that("Can compute emperical SE and MCMC error", {
  # True value
  tv <- 5
  # Simulated values
  set.seed(4445896)
  sv <- tv + rnorm(500)
  # Compute
  ese <- emperical_SE(sv)
  expect_equal(unname(round(ese, 2)), c(1.02, 0.03))
  expect_named(ese, c("emperical_se", "MCMC_SE"))
})

test_that("Can compute coverage and MCMC error", {
  # True value
  tv <- 5.3
  # Simulated values
  set.seed(45896)
  # Create mock CI values
  ci <- lapply(
    1:500,
    function(x) c(6 - runif(1), 6 + runif(1))
  )
  # Compute
  covered1 <- coverage(ci, tv)
  expect_equal(unname(round(covered1, 2)), c(0.3, 0.02))
  expect_named(covered1, c("coverage", "MCMC_SE"))
  # For second truth value
  # Coverage should be much higher
  tv2 <- 5.9
  covered2 <- coverage(ci, tv2)
  expect_equal(unname(round(covered2, 2)), c(0.91, 0.01))
  expect_named(covered2, c("coverage", "MCMC_SE"))
  # Expect coverage of scenario 2 to be higher than scenario 1
  expect_gt(tv2[1], tv[1])
})

test_that("Can compute MSE", {
  # True mean
  tm <- 10.8
  # Simulated values
  set.seed(367255)
  # Create mock values
  val1 <- rnorm(100, tm, 1)
  val2 <- rnorm(1000, tm, 1)
  # Compute RMSE
  out_mse1 <- round(MSE(val1, tm), 3)
  out_mse2 <- round(MSE(val2, tm), 3)
  expect_equal(out_mse1, 1.31)
  expect_equal(out_mse2, 1.009)
  expect_gt(out_mse1, out_mse2)
})
