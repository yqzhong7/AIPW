test_that("aipw input", {
  expect_error(aipw(exposure = rep(1,100),outcome = rep(1,100),tmle_fit = rep(1,100)))
})
