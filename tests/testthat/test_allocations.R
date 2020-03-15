context("Test allocation methods and class")

test_that("Can load allocations file from disk and convert to data frame",{
  # Load
  allocs <- read_allocations_file("allocations.json")
  expect_s3_class(allocs, "sleepsimR_allocations")
  expect_named(allocs, "abcd")
  expect_length(allocs, 1)
  expect_length(allocs[[1]], 4)
  expect_named(allocs[[1]], c("iteration_id", "ts_request", "ts_finished", "runtime_minutes"))
})
