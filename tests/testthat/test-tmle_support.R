#' @title Testing AIPW constructor: input data dimension
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2021/01/26
test_that("AIPW_tmle class: tmle", {
  require(tmle)
  require(SuperLearner)
  vec <- function() sample(0:1,100,replace = T)
  df <- data.frame(replicate(4,vec()))
  names(df) <- c("A","Y","W1","W2")
  tmle_fit <- tmle(Y=df$Y,
              A=df$A,
              W=df[,3:4],
              Q.SL.library="SL.glm",
              g.SL.library="SL.glm",
              family="binomial")
  #test constructor
  expect_error(AIPW_tmle$new(A=df$A,Y=df$Y,tmle_fit = "tmle_fit",verbose = T),
               info = "The tmle_fit is neither a `tmle` or `tmle3_Fit` object")
  expect_message(aipw_tmle <- AIPW_tmle$new(A=df$A,Y=df$Y,tmle_fit = tmle_fit,verbose = T),
                info = "Cross-fitting is supported only within the outcome model from a fitted tmle object (with cvQinit = TRUE)")
  #correctly print output
  expect_output(aipw_tmle$summary(), regexp = "Estimate")
  #check any null values after calculating results
  expect_false(any(sapply(aipw_tmle$estimates, is.null)))
  expect_false(any(sapply(aipw_tmle$libs, is.null)))
  expect_false(any(sapply(aipw_tmle$obs_est, is.null)))
  expect_false(is.null(aipw_tmle$result))
})
