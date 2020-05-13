#' @title Tesing calculate_result: SuperLeaner & k_split
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2020/05/13
test_that("AIPW calculate_result: SuperLeaner & k_split", {
  require(SuperLearner)
  ##k_split == 1: no sample splitting
  vec <- function() sample(0:1,100,replace = T)
  sl.lib <- c("SL.mean","SL.glm")
  aipw <-  AIPW$new(Y=vec(),
                    A=vec(),
                    W.Q =vec(),
                    W.g =vec(),
                    Q.SL.library=sl.lib,
                    g.SL.library=sl.lib,
                    k_split = 1,verbose = FALSE)
  aipw$fit()
  invisible(capture.output(est_mat <- aipw$calculate_result()))
  #correctly print output
  expect_output(aipw$calculate_result(), regexp = "Estimate")
  #whether output is a matrix
  expect_true(is.matrix(est_mat))
  #check any null values after calculating results
  expect_false(any(sapply(aipw$estimates, is.null)))
  expect_false(any(sapply(aipw$libs, is.null)))
  expect_false(any(sapply(aipw$obs_est, is.null)))
  expect_false(is.null(aipw$result))

  ##k_split >0: sample splitting == k_split
  aipw <-  AIPW$new(Y=vec(),
                    A=vec(),
                    W.Q =vec(),
                    W.g =vec(),
                    Q.SL.library=sl.lib,
                    g.SL.library=sl.lib,
                    k_split = 5,verbose = FALSE)
  aipw$fit()
  invisible(capture.output(est_mat <- aipw$calculate_result()))
  #correctly print output
  expect_output(aipw$calculate_result(), regexp = "Estimate")
  #whether output is a matrix
  expect_true(is.matrix(est_mat))
  #check any null values after calculating results
  expect_false(any(sapply(aipw$estimates, is.null)))
  expect_false(any(sapply(aipw$libs, is.null)))
  expect_false(any(sapply(aipw$obs_est, is.null)))
  expect_false(is.null(aipw$result))
})

#' @title Tesing calculate_result: sl3 & k_split
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2020/05/13
test_that("AIPW calculate_result: sl3 & k_split", {
  require(sl3)
  ##k_split == 1: no sample splitting
  vec <- function() sample(0:1,100,replace = T)
  lrnr_glm <- sl3::Lrnr_glm$new()
  lrnr_mean <- sl3::Lrnr_mean$new()
  stacklearner <- sl3::Stack$new(lrnr_glm, lrnr_mean)
  metalearner <- sl3::Lrnr_nnls$new()
  sl3.lib <- sl3::Lrnr_sl$new(learners = stacklearner,
                              metalearner = metalearner)
  aipw <-  AIPW$new(Y=vec(),
                    A=vec(),
                    W.Q =vec(),
                    W.g =vec(),
                    Q.SL.library=sl3.lib,
                    g.SL.library=sl3.lib,
                    k_split = 1,verbose = FALSE)
  aipw$fit()
  invisible(capture.output(est_mat <- aipw$calculate_result()))
  #correctly print output
  expect_output(aipw$calculate_result(), regexp = "Estimate")
  #whether output is a matrix
  expect_true(is.matrix(est_mat))
  #check any null values after calculating results
  expect_false(any(sapply(aipw$estimates, is.null)))
  expect_false(any(sapply(aipw$libs, is.null)))
  expect_false(any(sapply(aipw$obs_est, is.null)))
  expect_false(is.null(aipw$result))

  ##k_split >0: sample splitting == k_split
  aipw <-  AIPW$new(Y=vec(),
                    A=vec(),
                    W.Q =vec(),
                    W.g =vec(),
                    Q.SL.library=sl3.lib,
                    g.SL.library=sl3.lib,
                    k_split = 2,verbose = FALSE)
  aipw$fit()
  invisible(capture.output(est_mat <- aipw$calculate_result()))
  #correctly print output
  expect_output(aipw$calculate_result(), regexp = "Estimate")
  #whether output is a matrix
  expect_true(is.matrix(est_mat))
  #check any null values after calculating results
  expect_false(any(sapply(aipw$estimates, is.null)))
  expect_false(any(sapply(aipw$libs, is.null)))
  expect_false(any(sapply(aipw$obs_est, is.null)))
  expect_false(is.null(aipw$result))
})

#' @title Tesing calculate_result: verbose
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2020/05/13
test_that("AIPW calculate_result: verbose", {
  #verbose == TRUE: w/ progression bar & "Done!"
  vec <- function() sample(0:1,100,replace = T)
  sl.lib <- c("SL.mean","SL.glm")
  aipw <-  AIPW$new(Y=vec(),
                    A=vec(),
                    W.Q =vec(),
                    W.g =vec(),
                    Q.SL.library=sl.lib,
                    g.SL.library=sl.lib,
                    k_split = 1,verbose = T)
  expect_output(aipw$fit(),regexp = "Done")
})

