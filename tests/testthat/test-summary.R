#' @title Testing summary: SuperLeaner & k_split
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2020/08/09
test_that("AIPW summary: SuperLeaner & k_split", {
  require(SuperLearner)
  ##k_split == 1: no cross-fitting
  vec <- function() sample(0:1,100,replace = T)
  sl.lib <- c("SL.mean","SL.glm")
  aipw <-  AIPW$new(Y=vec(),
                    A=vec(),
                    W.Q =vec(),
                    W.g =vec(),
                    Q.SL.library=sl.lib,
                    g.SL.library=sl.lib,
                    k_split = 1,verbose = TRUE,
                    save.sl.fit = TRUE)
  #correctly print output
  expect_output(aipw$fit()$summary(), regexp = "Estimate")
  #check any null values after calculating results
  expect_false(any(sapply(aipw$estimates, is.null)))
  expect_false(any(sapply(aipw$libs, is.null)))
  expect_false(any(sapply(aipw$obs_est, is.null)))
  expect_false(is.null(aipw$result))

  ##k_split >0: cross-fitting == k_split
  aipw <-  AIPW$new(Y=vec(),
                    A=vec(),
                    W.Q =vec(),
                    W.g =vec(),
                    Q.SL.library=sl.lib,
                    g.SL.library=sl.lib,
                    k_split = 3,verbose = TRUE,
                    save.sl.fit = TRUE)
  #correctly print output
  expect_output(aipw$fit()$summary(), regexp = "Estimate")
  #check any null values after calculating results
  expect_false(any(sapply(aipw$estimates, is.null)))
  expect_false(any(sapply(aipw$libs, is.null)))
  expect_false(any(sapply(aipw$obs_est, is.null)))
  expect_false(is.null(aipw$result))
})



#' @title Testing summary: verbose
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2020/05/22
test_that("AIPW summary: verbose", {
  require(SuperLearner)
  vec <- function() sample(0:1,100,replace = T)
  sl.lib <- c("SL.mean","SL.glm")
  ##verbose==False
  aipw <-  AIPW$new(Y=vec(),
                    A=vec(),
                    W.Q =vec(),
                    W.g =vec(),
                    Q.SL.library=sl.lib,
                    g.SL.library=sl.lib,
                    k_split = 1,verbose = FALSE)
  expect_equal(capture.output(aipw$fit()$summary()),character(0))
})

#' @title Testing summary: g.bound input
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2021/01/26
test_that("Testing summary: g.bound input", {
  require(SuperLearner)
  vec <- function() sample(0:1,100,replace = T)
  sl.lib <- c("SL.mean","SL.glm")
  aipw <-  AIPW$new(Y=vec(),
                    A=vec(),
                    W.Q =vec(),
                    W.g =vec(),
                    Q.SL.library=sl.lib,
                    g.SL.library=sl.lib,
                    k_split = 1,verbose = FALSE)$fit()
  #Check error messages
  expect_warning(aipw$summary(c(0.1,0.2,0.3)),
                 info = "More than two g.bound are provided. Only the first two will be used.")
  expect_equal(aipw$.__enclos_env__$private$g.bound, c(0.1,0.2))
  expect_error(aipw$summary(0.5),
               info = "g.bound >= 0.5 is not allowed when only one g.bound value is provided")
  expect_error(aipw$summary("C"),
                 info = "g.bound must be numeric")
  expect_error(aipw$summary(-1),
               info = "g.bound must between 0 and 1")
  expect_error(aipw$summary(1),
               info = "g.bound must between 0 and 1")
})


#' @title Testing summary: continuous outcome reporting
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2021/03/24
test_that("AIPW summary: continuous outcome", {
  require(SuperLearner)
  ##k_split == 1: no cross-fitting
  vec <- function() sample(0:1,100,replace = T)
  sl.lib <- c("SL.mean","SL.glm")
  ##verbose==False
  aipw <-  AIPW$new(Y=rnorm(100,10),
                    A=vec(),
                    W =vec(),
                    Q.SL.library=sl.lib,
                    g.SL.library=sl.lib,
                    k_split = 1,verbose = FALSE)$fit()$summary()
  #check results with RR and OR
  expect_equal(nrow(aipw$result),3)
})


#' @title Testing summary: missing outcome reporting (N)
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2021/03/24
test_that("AIPW summary: missing outcome", {
  require(SuperLearner)
  ##k_split == 1: no cross-fitting
  vec <- function() sample(0:1,100,replace = T)
  sl.lib <- c("SL.mean","SL.glm")
  ##verbose==False
  expect_warning(aipw <-  AIPW$new(Y=c(NA,vec()[2:100]),
                    A=c(1,vec()[2:100]),
                    W =vec(),
                    Q.SL.library=sl.lib,
                    g.SL.library=sl.lib,
                    k_split = 1,verbose = FALSE)$fit()$summary())
  #Check N reporting
  expect_equal(aipw$result[3,5], 100)
})

#' @title Testing summary: ATT reporting
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2021/03/17
test_that("AIPW summary: continuous outcome", {
  require(SuperLearner)
  set.seed(123)
  vec <- function() sample(0:1,100,replace = T)
  sl.lib <- c("SL.mean","SL.glm")
  aipw <-  AIPW$new(Y=vec(),
                    A=vec(),
                    W.Q =vec(),
                    W.g =vec(),
                    Q.SL.library=sl.lib,
                    g.SL.library=sl.lib,
                    k_split = 1,verbose = FALSE)
  expect_warning(aipw$stratified_fit()$summary())
  #check results with RR and OR
  expect_equal(nrow(aipw$result),7)
})

