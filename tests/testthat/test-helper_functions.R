#' @title Testing helper functions
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2020/05/26
test_that("AIPW helper functions", {
  #get CI
  expect_equal(AIPW:::get_ci(est=0,se=1,ratio = F),c(-1.96,1.96))
  expect_equal(AIPW:::get_ci(est=1,se=1,ratio = T),exp(c(-1.96,1.96)))
  expect_warning(AIPW:::get_ci(est=-1,se=1,ratio = T))

  #private methods: get_RD, get_RR, get_OR, get_sigma_covar
  vec <- rep(1,100)
  root_n = sqrt(length(vec))
  sl.lib <- c("SL.mean","SL.glm")
  aipw <-  AIPW$new(Y=vec,
                    A=vec,
                    W.Q =vec,
                    W.g =vec,
                    Q.SL.library=sl.lib,
                    g.SL.library=sl.lib,
                    k_split = 1,verbose = FALSE)
  #get_RD
  expect_identical(as.numeric(aipw$.__enclos_env__$private$get_RD(aipw_eif0 = vec,
                                                                  aipw_eif1 = vec,
                                                                  root_n=root_n)),
               rep(0,4))
  #get_sigma_covar
  expect_identical(aipw$.__enclos_env__$private$get_sigma_covar(aipw_eif0 = vec,
                                                                aipw_eif1 = vec),
                   matrix(rep(0,4),ncol=2))
  #get_RR
  # sigma_covar <- matrix(c(0,0,0,0),ncol=2)
  # expect_identical(as.numeric(aipw$.__enclos_env__$private$get_RR(aipw_eif0 = vec,
  #                                                                 aipw_eif1 = vec,
  #                                                                 sigma_covar = sigma_covar,
  #                                                                 Z_norm=Z_norm)),
  #                  rep(0,4))
  #get_OR
  # expect_identical(as.numeric(aipw$.__enclos_env__$private$get_OR(aipw_eif0 = vec,
  #                                                                 aipw_eif1 = vec,
  #                                                                 sigma_covar = sigma_covar,
  #                                                                 Z_norm=Z_norm)),
  #                  rep(0,4))
  #.bound


  }
  )


#' @title Testing make new training index function
#'
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2021/01/27
test_that("make new training index: .new_cv_index",{
  vec <- rep(1,100)
  root_n = sqrt(length(vec))
  sl.lib <- c("SL.mean","SL.glm")
  aipw <-  AIPW$new(Y=vec,
                    A=vec,
                    W.Q =vec,
                    W.g =vec,
                    Q.SL.library=sl.lib,
                    g.SL.library=sl.lib,
                    k_split = 3,verbose = FALSE)
  suppressWarnings(aipw$fit())
  k_index = aipw$.__enclos_env__$private$cv$k_index
  fold_index = aipw$.__enclos_env__$private$cv$fold_index
  fold_length = aipw$.__enclos_env__$private$cv$fold_length
  expect_equal(sum(fold_length),100)
  k_split = aipw$.__enclos_env__$private$k_split
  for (i in 1:k_split){
    new_index = aipw$.__enclos_env__$private$.new_cv_index(val_fold=i,fold_length = fold_length, k_split=k_split)
    expect_equal(length(new_index),k_split-1) #same length of list
    expect_equal(as.numeric(unlist(new_index)), 1:length(as.numeric(unlist(new_index)))) #consecutive index
  }
})



#' @title Testing propensity score truncation function
#'
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2021/01/26
test_that("propensity score truncation: .bound",{
  vec <- rep(1,100)
  sl.lib <- c("SL.mean","SL.glm")
  aipw <-  AIPW$new(Y=vec,
                    A=vec,
                    W.Q =vec,
                    W.g =vec,
                    Q.SL.library=sl.lib,
                    g.SL.library=sl.lib,
                    k_split = 3,verbose = FALSE)

  ps <- seq(0, 1, length.out = 100)
  expect_equal(range(aipw$.__enclos_env__$private$.bound(ps, 0.1)), c(0.1, 0.9))
  expect_equal(range(aipw$.__enclos_env__$private$.bound(ps, 0.49)), c(0.49, 0.51))
  expect_equal(range(aipw$.__enclos_env__$private$.bound(ps, c(0.1,0.8))), c(0.1,0.8))
  expect_equal(range(aipw$.__enclos_env__$private$.bound(ps, c(0.8,0.1))), c(0.1,0.8))
  expect_equal(range(aipw$.__enclos_env__$private$.bound(ps, c(0.7,0.6))), c(0.6,0.7))
})
