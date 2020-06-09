#' @title Testing AIPW constructor: input data dimension
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2020/06/08
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
                info = "Sample splitting was not supported with a fitted tmle object")
  #correctly print output
  expect_output(aipw_tmle$summary(), regexp = "Estimate")
  #check any null values after calculating results
  expect_false(any(sapply(aipw_tmle$estimates, is.null)))
  expect_false(any(sapply(aipw_tmle$libs, is.null)))
  expect_false(any(sapply(aipw_tmle$obs_est, is.null)))
  expect_false(is.null(aipw_tmle$result))
})

test_that("AIPW_tmle class: tmle3", {
  require(tmle3)
  require(sl3)
  vec <- function() sample(0:1,100,replace = T)
  df <- data.frame(replicate(4,vec()))
  names(df) <- c("A","Y","W1","W2")
  #tmle3 setup
  node_list <- list(A = "A",Y = "Y",W = c("W1","W2"))
  or_spec <- tmle_OR(
    baseline_level = "0",
    contrast_level = "1"
  )
  tmle_task <- or_spec$make_tmle_task(df,node_list)
  lrnr_glm <- make_learner(Lrnr_glm)

  sl_Y <- Lrnr_sl$new(learners = list(lrnr_glm))
  sl_A <- Lrnr_sl$new(learners = list(lrnr_glm))

  learner_list <- list(A = sl_A, Y = sl_Y)
  tmle_fit <- tmle3(or_spec, data=df, node_list, learner_list)

  #test constructor
  expect_message(aipw_tmle <- AIPW_tmle$new(A=df$A,Y=df$Y,tmle_fit = tmle_fit,verbose = T),
                 info = "Propensity scores from fitted tmle3 object are by default truncated (0.025)")
  #correctly print output
  expect_output(aipw_tmle$summary(), regexp = "Estimate")
  #check any null values after calculating results
  expect_false(any(sapply(aipw_tmle$estimates, is.null)))
  expect_false(any(sapply(aipw_tmle$libs, is.null)))
  expect_false(any(sapply(aipw_tmle$obs_est, is.null)))
  expect_false(is.null(aipw_tmle$result))
})
