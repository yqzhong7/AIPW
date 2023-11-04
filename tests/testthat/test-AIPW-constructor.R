#' @title Testing AIPW constructor: input data dimension
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2020/05/28
test_that("AIPW constructor: input data dimension", {
  ##correct dimension
  #single column W.g
  vec <- rep(1,100)
  sl.lib <- c("SL.mean","SL.glm")
  expect_warning(aipw <-  AIPW$new(Y=vec,
                                   A=vec,
                                   W.Q =vec,
                                   W.g =vec,
                                   Q.SL.library=sl.lib,
                                   g.SL.library=sl.lib,
                                   k_split = 1,verbose = FALSE),
                 info = "Either `SuperLearner` or `sl3` package is not loaded.")
  expect_equal(aipw$n,100)
  expect_equal(length(aipw$.__enclos_env__$private$A),100)
  expect_equal(length(aipw$.__enclos_env__$private$Y),100)
  expect_equal(dim(aipw$.__enclos_env__$private$g.set)[1],100)
  expect_equal(dim(aipw$.__enclos_env__$private$Q.set),c(100,2))
  #multiple columns W.g
  mat <- matrix(rep(1,200),ncol = 2)
  expect_warning(aipw <-  AIPW$new(Y=vec,
                                   A=vec,
                                   W.Q =vec,
                                   W.g =mat,
                                   Q.SL.library=sl.lib,
                                   g.SL.library=sl.lib,
                                   k_split = 1,verbose = FALSE),
                 info = "Either `SuperLearner` or `sl3` package is not loaded.")
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
    regexp = "Please check the dimension of the data"
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
    regexp = "Please check the dimension of the covariates"
  )
})


#' @title Testing AIPW constructor: W input logic
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2020/05/28
test_that("AIPW constructor: W input logic", {
  ##correct dimension
  #single column W.g
  vec <- rep(1,100)
  sl.lib <- c("SL.mean","SL.glm")
  expect_error(aipw <-  AIPW$new(Y=vec,
                                 A=vec,
                                 Q.SL.library=sl.lib,
                                 g.SL.library=sl.lib,
                                 k_split = 1,verbose = FALSE),
               info="No sufficient covariates were provided.")
  expect_warning(aipw <-  AIPW$new(Y=vec,
                                   A=vec,
                                   W=vec,
                                   Q.SL.library=sl.lib,
                                   g.SL.library=sl.lib,
                                   k_split = 1,verbose = FALSE),
                 info = "Either `SuperLearner` or `sl3` package is not loaded.")
  expect_equal(aipw$n,100)
  expect_equal(length(aipw$.__enclos_env__$private$A),100)
  expect_equal(length(aipw$.__enclos_env__$private$Y),100)
  expect_equal(dim(aipw$.__enclos_env__$private$g.set)[1],100)
  expect_equal(dim(aipw$.__enclos_env__$private$Q.set),c(100,2))
  #multicolumns
  expect_warning(aipw <-  AIPW$new(Y=vec,
                                   A=vec,
                                   W=cbind(vec,vec),
                                   Q.SL.library=sl.lib,
                                   g.SL.library=sl.lib,
                                   k_split = 1,verbose = FALSE),
                 info = "Either `SuperLearner` or `sl3` package is not loaded.")
  expect_equal(dim(aipw$.__enclos_env__$private$g.set),c(100,2))
  expect_equal(dim(aipw$.__enclos_env__$private$Q.set),c(100,3))

})

#' @title Testing AIPW constructor: SL libraries
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2020/05/15
test_that("AIPW constructor: SL libraries", {
  vec <- rep(1,100)
  sl.lib <- c("SL.mean","SL.glm")
  ##SuperLearner
  #if SuperLearner package is not loaded
  expect_warning(aipw <-  AIPW$new(Y=vec,
                    A=vec,
                    W.Q =vec,
                    W.g =vec,
                    Q.SL.library=sl.lib,
                    g.SL.library=sl.lib,
                    k_split = 1,verbose = FALSE),
  regexp = "Either `SuperLearner` or `sl3` package is not loaded.")
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
})


#' @title Testing AIPW constructor: k_split
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2021/01/26
test_that("AIPW constructor: k_split", {
  require("SuperLearner")
  #cross-fitting
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
    regexp = "`k_split` < 1 is not allowed."
  )
  expect_error(
    AIPW$new(Y=rep(1,100),
             A=rep(1,100),
             W.Q =rep(1,100),
             W.g =rep(1,100),
             Q.SL.library=c("SL.mean","SL.glm"),
             g.SL.library=c("SL.mean","SL.glm"),
             k_split = 100,verbose = FALSE),
    regexp = "`k_split` >= number of observations is not allowed."
  )
})


#' @title Testing AIPW constructor: verbose
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2020/05/09
test_that("AIPW constructor: verbose", {
  require("SuperLearner")
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
    regexp = "`verbose` is not valid"
  )
})



#' @title Testing AIPW constructor: missing exposure and outcome
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2021/03/24
test_that("AIPW constructor: input data dimension", {
  vec <- rep(1,100)
  vec_na <- c(NA,vec[2:100])
  sl.lib <- c("SL.mean","SL.glm")
  #missing exposure
  expect_error(aipw <-  AIPW$new(Y=vec,
                                   A=vec_na,
                                   W.Q =vec,
                                   W.g =vec,
                                   Q.SL.library=sl.lib,
                                   g.SL.library=sl.lib,
                                   k_split = 1,verbose = FALSE),
                 info = "Missing exposure is not allowed.")
  #missing outcome w/ W.Q and W.g specified
  expect_error(suppressWarnings({aipw <-  AIPW$new(Y=vec_na,
                                   A=vec,
                                   W.Q =vec,
                                   W.g =vec,
                                   Q.SL.library=sl.lib,
                                   g.SL.library=sl.lib,
                                   k_split = 1,verbose = FALSE)}),
                 info = "`W.Q` and `W.g` are disabled when missing outcome is detected. Please provide covariates in `W`")
  #missing outcome working
  expect_warning(aipw <-  AIPW$new(Y=vec_na,
                                 A=vec,
                                 W =vec,
                                 Q.SL.library=sl.lib,
                                 g.SL.library=sl.lib,
                                 k_split = 1,verbose = FALSE),
               info = "Missing outcome is detected. Analysis assumes missing at random (MAR).")
  expect_equal(aipw$n,100)
  expect_equal(length(aipw$.__enclos_env__$private$A),100)
  expect_equal(length(aipw$.__enclos_env__$private$Y),100)
  expect_equal(dim(aipw$.__enclos_env__$private$Q.set)[1],100)
  expect_equal(dim(aipw$.__enclos_env__$private$g.set)[1],100)
  expect_equal(sum(aipw$.__enclos_env__$private$observed),99)
  expect_equal(length(aipw$.__enclos_env__$private$observed),100)
})
