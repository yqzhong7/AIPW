#' @title Testing fit: SuperLeaner & k_split
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2020/05/13
test_that("AIPW fit: SuperLeaner & k_split", {
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
  #check any null values after calculating results
  expect_false(any(sapply(aipw$libs, is.null)))
  expect_false(any(sapply(aipw$obs_est[1:4], is.na))) #mu - raw_p_score
  expect_true(any(sapply(aipw$obs_est[5:7], is.null)))
  expect_true(is.null(aipw$result))
  expect_true(is.null(aipw$estimate))

  ##k_split >0: sample splitting == k_split
  aipw <-  AIPW$new(Y=vec(),
                    A=vec(),
                    W.Q =vec(),
                    W.g =vec(),
                    Q.SL.library=sl.lib,
                    g.SL.library=sl.lib,
                    k_split = 2,verbose = FALSE)
  aipw$fit()
  expect_false(any(sapply(aipw$libs, is.null)))
  expect_false(any(sapply(aipw$obs_est[1:4], is.na))) #mu - raw_p_score
  expect_true(any(sapply(aipw$obs_est[5:7], is.null)))
  expect_true(is.null(aipw$result))
  expect_true(is.null(aipw$estimate))
})

#' @title Testing calculate_result: sl3 & k_split
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2020/05/09
test_that("AIPW fit: sl3 & k_split", {
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
  expect_false(any(sapply(aipw$libs, is.null)))
  expect_false(any(sapply(aipw$obs_est[1:4], is.na))) #mu - raw_p_score
  expect_true(any(sapply(aipw$obs_est[5:7], is.null)))
  expect_true(is.null(aipw$result))
  expect_true(is.null(aipw$estimate))

  ##k_split >0: sample splitting == k_split
  aipw <-  AIPW$new(Y=vec(),
                    A=vec(),
                    W.Q =vec(),
                    W.g =vec(),
                    Q.SL.library=sl3.lib,
                    g.SL.library=sl3.lib,
                    k_split = 2,verbose = FALSE)
  aipw$fit()
  expect_false(any(sapply(aipw$libs, is.null)))
  expect_false(any(sapply(aipw$obs_est[1:4], is.na))) #mu - raw_p_score
  expect_true(any(sapply(aipw$obs_est[5:7], is.null)))
  expect_true(is.null(aipw$result))
  expect_true(is.null(aipw$estimate))
})


#' @title Testing fit: verbose
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2020/05/22
test_that("AIPW fit: verbose", {
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
  expect_output(aipw$fit(),regexp = "=====")
})

