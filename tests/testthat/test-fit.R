#' @title Testing fit: SuperLeaner & k_split
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2020/08/09
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
                    k_split = 1,verbose = FALSE,
                    save.sl.fit = TRUE)
  aipw$fit()
  #check any null values after calculating results
  expect_false(any(sapply(aipw$libs, is.null)))
  expect_false(any(sapply(aipw$obs_est[1:4], is.na))) #mu - raw_p_score
  expect_true(any(sapply(aipw$obs_est[5:7], is.null)))
  expect_true(is.null(aipw$result))
  expect_true(is.null(aipw$estimate))

  expect_error(aipw <-  AIPW$new(Y=vec(),
                                 A=vec(),
                                 W.Q =vec(),
                                 W.g =vec(),
                                 Q.SL.library=sl.lib,
                                 g.SL.library=sl.lib,
                                 k_split = 2,verbose = FALSE))

  ##k_split >=3: sample splitting == k_split
  aipw <-  AIPW$new(Y=vec(),
                    A=vec(),
                    W.Q =vec(),
                    W.g =vec(),
                    Q.SL.library=sl.lib,
                    g.SL.library=sl.lib,
                    k_split = 3,verbose = FALSE)


  aipw <-  AIPW$new(Y=vec(),
                    A=vec(),
                    W.Q =vec(),
                    W.g =vec(),
                    Q.SL.library=sl.lib,
                    g.SL.library=sl.lib,
                    k_split = 3,verbose = FALSE,
                    save.sl.fit = TRUE)
  aipw$fit()
  expect_false(any(sapply(aipw$libs, is.null)))
  expect_false(any(sapply(aipw$obs_est[1:4], is.na))) #mu - raw_p_score
  expect_true(any(sapply(aipw$obs_est[5:7], is.null)))
  expect_true(is.null(aipw$result))
  expect_true(is.null(aipw$estimate))
})

#' @title Testing summary: sl3 & k_split
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2020/08/09
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
                    k_split = 1,verbose = FALSE,
                    save.sl.fit = TRUE)
  aipw$fit()
  expect_false(any(sapply(aipw$libs, is.null)))
  expect_false(any(sapply(aipw$obs_est[1:4], is.na))) #mu - raw_p_score
  expect_true(any(sapply(aipw$obs_est[5:7], is.null)))
  expect_true(is.null(aipw$result))
  expect_true(is.null(aipw$estimate))

  ##k_split >=3: sample splitting == k_split
  aipw <-  AIPW$new(Y=vec(),
                    A=vec(),
                    W.Q =vec(),
                    W.g =vec(),
                    Q.SL.library=sl3.lib,
                    g.SL.library=sl3.lib,
                    k_split = 3,verbose = FALSE,
                    save.sl.fit = TRUE)
  aipw$fit()
  expect_false(any(sapply(aipw$libs, is.null)))
  expect_false(any(sapply(aipw$obs_est[1:4], is.na))) #mu - raw_p_score
  expect_true(any(sapply(aipw$obs_est[5:7], is.null)))
  expect_true(is.null(aipw$result))
  expect_true(is.null(aipw$estimate))
})


#' @title Testing fit: verbose and progressr
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2021/01/20
test_that("AIPW fit: verbose", {
  #verbose == TRUE: "Done!"
  library(SuperLearner)
  vec <- function() sample(0:1,100,replace = T)
  sl.lib <- c("SL.mean","SL.glm")
  aipw <-  AIPW$new(Y=vec(),
                    A=vec(),
                    W.Q =vec(),
                    W.g =vec(),
                    Q.SL.library=sl.lib,
                    g.SL.library=sl.lib,
                    k_split = 1,verbose = T)
  expect_message(aipw$fit(),regexp = "Done!")

  ##progressr
  #when progressr not loaded
  expect_false(aipw$.__enclos_env__$private$isLoaded_progressr)
  #when loaded
  library(progressr)
  aipw$fit()
  expect_true(aipw$.__enclos_env__$private$isLoaded_progressr)
})


