#' @title Tesing AIPW constructor: input data dimension
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2020/05/09
test_that("AIPW constructor: input data dimension", {
  ##correct dimension
  #single column W.g
  vec <- rep(1,100)
  sl.lib <- c("SL.mean","SL.glm")
  aipw <-  AIPW$new(Y=vec,
             A=vec,
             W.Q =vec,
             W.g =vec,
             Q.SL.library=sl.lib,
             g.SL.library=sl.lib,
             k_split = 1,verbose = FALSE)
  expect_equal(aipw$n,100)
  expect_equal(length(aipw$.__enclos_env__$private$A),100)
  expect_equal(length(aipw$.__enclos_env__$private$Y),100)
  expect_equal(dim(aipw$.__enclos_env__$private$g.set)[1],100)
  expect_equal(dim(aipw$.__enclos_env__$private$Q.set),c(100,2))
  #multiple columns W.g
  mat <- matrix(rep(1,200),ncol = 2)
  aipw <-  AIPW$new(Y=vec,
                    A=vec,
                    W.Q =vec,
                    W.g =mat,
                    Q.SL.library=sl.lib,
                    g.SL.library=sl.lib,
                    k_split = 1,verbose = FALSE)
  expect_equal(dim(aipw$.__enclos_env__$private$g.set),c(100,2))
  ##wrong dimension
  #single column W.g (A and W.Q have the same nrow)
  expect_error(
    AIPW$new(Y=vec,
             A=rep(1,80),
             W.Q =rep(1,80),
             W.g =vec,
             Q.SL.library=sl.lib,
             g.SL.library=sl.lib,
             k_split = 1,verbose = FALSE),
    regexp = "Please check the dimension of the data"
  )
  #A and W.Q have different nrow
  expect_error(
    AIPW$new(Y=vec,
             A=rep(1,80),
             W.Q =vec,
             W.g =vec,
             Q.SL.library=sl.lib,
             g.SL.library=sl.lib,
             k_split = 1,verbose = FALSE),
    regexp = "arguments imply differing number of rows"
  )
  #multiple columns W.g
  expect_error(
    AIPW$new(Y=vec,
             A=vec,
             W.Q =vec,
             W.g =matrix(rep(1,200),ncol = 1),
             Q.SL.library=sl.lib,
             g.SL.library=sl.lib,
             k_split = 1,verbose = FALSE),
    regexp = "Please check the dimension of the data"
  )
})


#' @title Tesing AIPW constructor: SL libraries
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2020/05/09
test_that("AIPW constructor: SL libraries", {
  vec <- rep(1,100)
  sl.lib <- c("SL.mean","SL.glm")
  ##SuperLearner
  aipw <-  AIPW$new(Y=vec,
                    A=vec,
                    W.Q =vec,
                    W.g =vec,
                    Q.SL.library=sl.lib,
                    g.SL.library=sl.lib,
                    k_split = 1,verbose = FALSE)
  #sl.library writing
  expect_identical(aipw$libs$Q.SL.library,sl.lib)
  expect_identical(aipw$libs$g.SL.library,sl.lib)
  expect_false(is.null(aipw$sl.fit))
  expect_false(is.null(aipw$sl.predict))
  #wrong SL library
  expect_error(
    AIPW$new(Y=vec,
             A=vec,
             W.Q =vec,
             W.g =vec,
             Q.SL.library=c("screen.randomForest"),
             g.SL.library=sl.lib,
             k_split = 1,verbose = FALSE),
    regexp = "Input Q.SL.library and/or g.SL.library is not a valid SuperLearner library"
  )

  #sl3
  lrnr_glm <- sl3::Lrnr_glm$new()
  lrnr_mean <- sl3::Lrnr_mean$new()
  stacklearner <- sl3::Stack$new(lrnr_glm, lrnr_mean)
  metalearner <- sl3::Lrnr_nnls$new()
  sl3.lib <- sl3::Lrnr_sl$new(learners = stacklearner,
                              metalearner = metalearner)
  aipw <- AIPW$new(Y=vec,
             A=vec,
             W.Q =vec,
             W.g =vec,
             Q.SL.library=sl3.lib,
             g.SL.library=sl3.lib,
             k_split = 1,verbose = FALSE)
  #sl3 lib writing
  expect_identical(aipw$libs$Q.SL.library,sl3.lib)
  expect_identical(aipw$libs$g.SL.library,sl3.lib)
  expect_false(is.null(aipw$sl.fit))
  expect_false(is.null(aipw$sl.predict))
  #warning for using stack learners only
  expect_warning(AIPW$new(Y=vec,
                          A=vec,
                          W.Q =vec,
                          W.g =vec,
                          Q.SL.library=stacklearner,
                          g.SL.library=stacklearner,
                          k_split = 1,verbose = FALSE),
                 regexp = "sl3::Stack")
  # wrong libs
  expect_error(AIPW$new(Y=vec,
                        A=vec,
                        W.Q =vec,
                        W.g =vec,
                        Q.SL.library=sl.lib,
                        g.SL.library=sl3.lib,
                        k_split = 1,verbose = FALSE),
               regexp = "Input Q.SL.library and/or g.SL.library is not a valid SuperLearner/sl3 library")
})


#' @title Tesing AIPW constructor: k_split
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2020/05/09
test_that("AIPW constructor: k_split", {
  #sample splitting
  vec <- rep(1,100)
  sl.lib <- c("SL.mean","SL.glm")
  aipw <-  AIPW$new(Y=vec,
                    A=vec,
                    W.Q =vec,
                    W.g =vec,
                    Q.SL.library=sl.lib,
                    g.SL.library=sl.lib,
                    k_split = 5,verbose = FALSE)
  expect_identical(aipw$.__enclos_env__$private$k_split,5)
  #k_split out of range
  expect_error(
    AIPW$new(Y=rep(1,100),
             A=rep(1,100),
             W.Q =rep(1,100),
             W.g =rep(1,100),
             Q.SL.library=c("SL.mean","SL.glm"),
             g.SL.library=c("SL.mean","SL.glm"),
             k_split = -1,verbose = FALSE),
    regexp = "k_split is not valid"
  )
  expect_error(
    AIPW$new(Y=rep(1,100),
             A=rep(1,100),
             W.Q =rep(1,100),
             W.g =rep(1,100),
             Q.SL.library=c("SL.mean","SL.glm"),
             g.SL.library=c("SL.mean","SL.glm"),
             k_split = 100,verbose = FALSE),
    regexp = "k_split is not valid"
  )
})


#' @title Tesing AIPW constructor: verbose
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2020/05/09
test_that("AIPW constructor: verbose", {
  vec <- rep(1,100)
  sl.lib <- c("SL.mean","SL.glm")
  aipw <-  AIPW$new(Y=vec,
                    A=vec,
                    W.Q =vec,
                    W.g =vec,
                    Q.SL.library=sl.lib,
                    g.SL.library=sl.lib,
                    k_split = 5,verbose = TRUE)
  expect_true(aipw$.__enclos_env__$private$verbose)
  #wrong verbose value
  expect_error(
    AIPW$new(Y=rep(1,100),
             A=rep(1,100),
             W.Q =rep(1,100),
             W.g =rep(1,100),
             Q.SL.library=c("SL.mean","SL.glm"),
             g.SL.library=c("SL.mean","SL.glm"),
             k_split = 5,verbose = -1),
    regexp = "verbose is not valid"
  )
})
