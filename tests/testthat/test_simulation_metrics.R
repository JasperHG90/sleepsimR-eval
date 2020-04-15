context("Test simulation metrics (parameter bias etc)")

test_that("Can compute bias and MCMC error for bias", {
  # True value
  tv <- 5
  # Simulated values
  set.seed(4445896)
  noise <- rnorm(500)
  sv <- tv + noise
  # Compute
  b <- bias(tv, sv)
  # Bias should be equal to mean noise value
  expect_equal(unname(round(b, 5)), c(round(mean(noise), 5), 0.04546))
  expect_named(b, c("bias", "MCMC_SE"))
  # Check if MCMC SE goes down as number of iterations increases
  noise <- rnorm(1500)
  sv <- tv + noise
  b2 <- bias(tv, sv)
  expect_lt(b2[2], b[2])
})

test_that("Can compute empirical SE and MCMC error", {
  # True value
  tv <- 5
  # Simulated values
  set.seed(4445896)
  sv <- tv + rnorm(500)
  # Compute
  ese <- empirical_SE(sv)
  expect_equal(unname(round(ese, 2)), c(1.02, 0.03))
  expect_named(ese, c("empirical_se", "MCMC_SE"))
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
  out_mse1 <- unname(round(MSE(val1, tm), 3))
  out_mse2 <- unname(round(MSE(val2, tm), 3))
  expect_equal(out_mse1[1], 1.31)
  expect_equal(out_mse2[1], 1.009)
  expect_gt(out_mse1[1], out_mse2[2])
  # MCMC error should be larger
  expect_gt(out_mse1[2], out_mse2[2])
})

test_that("Can compute modSE", {
  # True mean
  tm <- 10.8
  # Simulated values
  set.seed(367255)
  # Create mock values

  # Scenario 1
  # 100 samples of 100 observations each
  # SE should be larg here (smallest sample)
  # MCMC SE should be largest here (smallest number of iterations)
  mv1 <- unlist(lapply(1:100, function(x) {
    val1 <- rnorm(100, tm, 1)
    # Compute SE
    return(sd(val1) / sqrt(length(val1)))
  }))

  # Scenario 2
  # 1000 samples of 100 observations each
  # SE should be ~same as in scenario 1 (sample size does not change).
  # MCMC SE should be smaller than scenario 1 (many more samples of same sample size)
  mv2 <- unlist(lapply(1:1000, function(x) {
    val1 <- rnorm(100, tm, 1)
    # Compute SE
    return(sd(val1) / sqrt(length(val1)))
  }))

  # Scenario 3
  # 100 samples of 1000 observations each
  # SE should be smaller than scenarios 1 & 2 (sample size is larger).
  # MCMC SE should be about the same relative to the SE value (same number of iterations)
  mv3 <- unlist(lapply(1:100, function(x) {
    val1 <- rnorm(1000, tm, 1)
    # Compute SE
    return(sd(val1) / sqrt(length(val1)))
  }))

  # Scenario 4
  # 1000 samples of 1000 observations each
  # SE should be ~same as in scenario 3 (sample size is the same).
  # MCMC SE should be smaller than scenario 3 (larger number of iterations)
  mv4 <- unlist(lapply(1:1000, function(x) {
    val1 <- rnorm(1000, tm, 1)
    # Compute SE
    return(sd(val1) / sqrt(length(val1)))
  }))

  # Compute modSE
  out_modse1 <- round(unname(average_model_SE(mv1)), 5)
  out_modse2 <- round(unname(average_model_SE(mv2)), 5)
  out_modse3 <- round(unname(average_model_SE(mv3)), 5)
  out_modse4 <- round(unname(average_model_SE(mv4)), 5)

  # Expectations
  expect_equal(round(out_modse1[1],2), round(out_modse2[1], 2))
  expect_gt(out_modse1[2], out_modse2[2])
  expect_lt(out_modse3[1], out_modse2[1])
  expect_equal(round(out_modse3[2] / out_modse3[1], 3),
               round(out_modse2[2] / out_modse2[1], 3))
  expect_equal(round(out_modse4[1],3), round(out_modse3[1],3))
  expect_lt(out_modse4[2], out_modse3[2])
})
