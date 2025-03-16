library(testthat)
library(SuperLearner)

# Helper function to create test data
create_test_data <- function(n = 100, binary = TRUE) {
  W <- matrix(rnorm(n * 2), ncol = 2)
  A <- rbinom(n, 1, 0.5)
  if (binary) {
    Y <- rbinom(n, 1, plogis(-1 + 0.5*A + 0.3*W[,1] + 0.2*W[,2]))
  } else {
    Y <- -1 + 0.5*A + 0.3*W[,1] + 0.2*W[,2] + rnorm(n, 0, 0.1)
  }
  return(list(Y = Y, A = A, W = W))
}

# Test constructor
test_that("Repeated constructor works correctly", {
    set.seed(123)
  # Create test data
  test_data <- create_test_data()

  # Create AIPW object
  aipw_obj <- AIPW$new(Y = test_data$Y,
                       A = test_data$A,
                       W = test_data$W,
                       Q.SL.library = "SL.mean",
                       g.SL.library = "SL.mean",
                       k_split = 2,
                       verbose = FALSE)

  # Test valid initialization
  repeated_obj <- Repeated$new(aipw_obj)
  expect_s3_class(repeated_obj, "R6")
  expect_true("RepeatedFit" %in% class(repeated_obj))
  expect_identical(repeated_obj$aipw_obj, aipw_obj)

  # Test invalid input type
  expect_error(Repeated$new("not an AIPW object"))
})

# Test repfit method
test_that("repfit method works correctly", {
  set.seed(123)
  test_data <- create_test_data()
  aipw_obj <- AIPW$new(Y = test_data$Y,
                       A = test_data$A,
                       W = test_data$W,
                       Q.SL.library = "SL.mean",
                       g.SL.library = "SL.mean",
                       k_split = 2,
                       verbose = FALSE)
  repeated_obj <- Repeated$new(aipw_obj)

  # Test with different num_reps
  repeated_obj$repfit(num_reps = 3, stratified = FALSE)
  expect_equal(nrow(repeated_obj$repeated_estimates) %% 3, 0)
  expect_true(all(c("Estimate", "SE", "Estimand") %in% colnames(repeated_obj$repeated_estimates)))

  # Test with stratified = TRUE
  repeated_obj$repfit(num_reps = 2, stratified = TRUE)
  expect_true(repeated_obj$stratified_fitted)
  expect_equal(nrow(repeated_obj$repeated_estimates) %% 2, 0)

  # Test with continuous outcome
  cont_data <- create_test_data(binary = FALSE)
  aipw_cont <- AIPW$new(Y = cont_data$Y,
                        A = cont_data$A,
                        W = cont_data$W,
                        Q.SL.library = "SL.mean",
                        g.SL.library = "SL.mean",
                        k_split = 2,
                        verbose = FALSE)
  repeated_cont <- Repeated$new(aipw_cont)
  repeated_cont$repfit(num_reps = 2, stratified = FALSE)
  expect_true(all(c("Mean of Exposure", "Mean of Control", "Mean Difference") %in%
                   unique(repeated_cont$repeated_estimates$Estimand)))
})

# Test summary_median method
test_that("summary_median method works correctly", {
  set.seed(123)
  test_data <- create_test_data()
  aipw_obj <- AIPW$new(Y = test_data$Y,
                       A = test_data$A,
                       W = test_data$W,
                       Q.SL.library = "SL.mean",
                       g.SL.library = "SL.mean",
                       k_split = 2,
                       verbose = FALSE)
  repeated_obj <- Repeated$new(aipw_obj)
  repeated_obj$repfit(num_reps = 3, stratified = FALSE)
  repeated_obj$summary_median()

  # Check result structure
  expect_true(!is.null(repeated_obj$result))
  expect_true(!is.null(repeated_obj$repeated_results))

  # Check column names
  expected_cols <- c("Median Estimate", "Median SE", "SE of Median Estimate",
                    "95% LCL Median Estimate", "95% UCL Median Estimate")
  expect_true(all(expected_cols %in% colnames(repeated_obj$result)))

  # Check CIs contain point estimates
  expect_true(all(repeated_obj$result[, "Median Estimate"] >= repeated_obj$result[, "95% LCL Median Estimate"]))
  expect_true(all(repeated_obj$result[, "Median Estimate"] <= repeated_obj$result[, "95% UCL Median Estimate"]))
})

# Test full workflow and consistency
test_that("Full workflow produces consistent results", {
  set.seed(123)
  test_data <- create_test_data()
  aipw_obj <- AIPW$new(Y = test_data$Y,
                       A = test_data$A,
                       W = test_data$W,
                       Q.SL.library = "SL.mean",
                       g.SL.library = "SL.mean",
                       k_split = 2,
                       verbose = FALSE)

  # Run twice with same seed
  set.seed(456)
  repeated_obj1 <- Repeated$new(aipw_obj)
  repeated_obj1$repfit(num_reps = 3, stratified = FALSE)
  repeated_obj1$summary_median()
  result1 <- repeated_obj1$result

  set.seed(456)
  repeated_obj2 <- Repeated$new(aipw_obj)
  repeated_obj2$repfit(num_reps = 3, stratified = FALSE)
  repeated_obj2$summary_median()
  result2 <- repeated_obj2$result

  # Results should be identical with same seed
  expect_equal(result1, result2)

  # Test with different SuperLearner libraries
  aipw_obj_sl <- AIPW$new(Y = test_data$Y,
                          A = test_data$A,
                          W = test_data$W,
                          Q.SL.library = c("SL.mean", "SL.glm"),
                          g.SL.library = c("SL.mean", "SL.glm"),
                          k_split = 2,
                          verbose = FALSE)

  repeated_obj_sl <- Repeated$new(aipw_obj_sl)
  expect_error(repeated_obj_sl$repfit(num_reps = 2, stratified = FALSE), NA)
  expect_error(repeated_obj_sl$summary_median(), NA)
})
