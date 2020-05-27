#' @title Augmented Inverse Probablity Weighting (AIPW)
#'
#' @description Define an R6Class aipw object
#'
#' @docType class
#'
#' @importFrom R6 R6Class
#'
#' @export
#'
#' @details create an AIPW object
#'
#' @return \code{AIPW} object
#'
#' @format \code{\link{R6Class}} object.
AIPW <- R6::R6Class(
  "AIPW",
  portable = TRUE,
  inherit = AIPW_base,
  public = list(
    #' @field libs SuperLearner or sl3 libraries and their fitted objects
    libs =list(Q.SL.library=NULL,
               Q.fit = NULL,
               g.SL.library=NULL,
               g.fit = NULL,
               num_val_index = NULL),
    #' @field sl.fit a wrapper for fitting SuperLearner or sl3
    sl.fit = NULL,
    #' @field sl.predict a wrapper using \code{sl.fit} to predict
    sl.predict = NULL,

    #' @description
    #' Create a new `AIPW` object.
    #'
    #' @param Y outcome (binary integer: 0 or 1)
    #' @param A exposure (binary integer: 0 or 1)
    #' @param W.Q covariates for outcome model (vector, matrix or data.frame)
    #' @param W.g covariates for exposure model (vector, matrix or data.frame)
    #' @param Q.SL.library SuperLearner libraries or sl3 learner object (Lrnr_base) for outcome model
    #' @param g.SL.library SuperLearner libraries or sl3 learner object (Lrnr_base) for exposure model
    #' @param k_split number of splitting (integer; range: from 1 to number of observation-1):
    #'   if k_split=1, no sample splitting;
    #'   if k_split>1, use similar technique of cross-validation
    #'   (e.g., k_split=5, use 4/5 of the data to estimate and the remaining 1/5 leftover to predict)
    #' @param verbose whether to show progression bar and print the result (logical; Default = FALSE)
    #'
    #' @return A new `AIPW` obejct
    #'
    #' @examples
    #' library(SuperLearner)
    #' aipw_sl <- AIPW$new(Y=rbinom(100,1,0.5), A=rbinom(100,1,0.5),
    #'                     W.Q=rbinom(100,1,0.5), W.g=rbinom(100,1,0.5),
    #'                     Q.SL.library="SL.mean",g.SL.library="SL.mean",
    #'                     k_split=1,verbose=FALSE)
    initialize = function(Y=NULL, A=NULL,W.Q=NULL, W.g=NULL,
                          Q.SL.library=NULL,g.SL.library=NULL,
                          k_split=1,verbose=FALSE){
      #initialize from AIPW_base class
      super$initialize(Y=Y,A=A)
      #save input into private fields
      private$Q.set=cbind(A, as.data.frame(W.Q))
      private$g.set=as.data.frame(W.g)
      private$k_split=k_split
      private$verbose=verbose
      #check data length
      if (!(length(private$Y)==dim(private$Q.set)[1] & length(private$A)==dim(private$g.set)[1])){
        stop("Please check the dimension of the data")
      }
      #determine SuperLearner or sl3 and change accordingly
      if (is.character(Q.SL.library) & is.character(g.SL.library)) {
        if (any(grepl("SL.",Q.SL.library)) & any(grepl("SL.",g.SL.library))){
          #change future package loading
          private$sl.pkg <- "SuperLearner"
          #change wrapper functions
          self$sl.fit = function(Y, X, SL.library){
            fit <- SuperLearner::SuperLearner(Y = Y, X = X, SL.library = SL.library, family="binomial")
            return(fit)
          }
          self$sl.predict = function(fit, newdata){
            pred <- as.numeric(predict(fit,newdata = newdata)$pred)
            return(pred)
          }
        } else{
          stop("Input Q.SL.library and/or g.SL.library is not a valid SuperLearner library")
        }
      } else if (any(class(Q.SL.library) == "Lrnr_base") & any(class(g.SL.library) == "Lrnr_base")) {
        #only using Stack in sl3 will return estimates of each library separately
        if (any(class(Q.SL.library) == "Stack") & any(class(g.SL.library) == "Stack")){
          warning("Only using sl3::Stack may cause problem. Please consider using metalearners for the stacked libraries!")
        } else {
          #change wrapper functions
          self$sl.fit = function(X, Y, SL.library){
            dat <- data.frame(cbind(Y,X))
            dat_colnames <- colnames(dat)
            task <- sl3::sl3_Task$new(dat, covariates = colnames(dat)[-1],
                                      outcome = colnames(dat)[1], outcome_type = "binomial"
            )
            fit <- SL.library$train(task)
            return(fit)
          }
          self$sl.predict = function(fit, newdata){
            new_task <- sl3::sl3_Task$new(newdata, covariates = colnames(newdata))
            pred <- fit$predict(new_task)
            return(pred)
          }
        }
      } else {
        stop("Input Q.SL.library and/or g.SL.library is not a valid SuperLearner/sl3 library")
      }

      #input sl libraries
      self$libs$Q.SL.library=Q.SL.library
      self$libs$g.SL.library=g.SL.library
      #validation set index
      self$libs$num_val_index <- rep(NA,self$n)
      #check k_split value
      if (private$k_split<1 | private$k_split>=self$n){
        stop("k_split is not valid")
      }
      #check verbose value
      if (!is.logical(private$verbose)){
        stop("verbose is not valid")
      }
      #check if SuperLearner and/or sl3 library is loaded
      if (!any(names(sessionInfo()$otherPkgs) %in% c("SuperLearner","sl3"))){
        warning("Either `SuperLearner` or `sl3` package is not loaded.")
      }
      #check if future.apply is loaded otherwise lapply would be used.
      if (any(names(sessionInfo()$otherPkgs) %in% c("future.apply"))){
        private$.f_lapply = function(iter,func) future.apply::future_lapply(iter,func,future.seed = T,future.packages = private$sl.pkg)
      }else{
        private$.f_lapply = function(iter,func) lapply(iter,func)
        }
    },
    #' @description
    #' Fitting the data into the `AIPW` object with/without sample splitting to estimate the influence functions
    #'
    #' @return A fitted `AIPW` obejct
    #'
    #' @examples
    #' library(SuperLearner)
    #' aipw_sl <- AIPW$new(Y=rbinom(100,1,0.5), A=rbinom(100,1,0.5),
    #'                     W.Q=rbinom(100,1,0.5), W.g=rbinom(100,1,0.5),
    #'                     Q.SL.library="SL.mean",g.SL.library="SL.mean",
    #'                     k_split=1,verbose=FALSE)
    #' aipw_sl$fit()
    fit = function(){
      #create index for sample splitting
      k_index <- sample(rep(1:private$k_split,ceiling(self$n/private$k_split))[1:self$n],replace = F)
      #progress bar setup
      progressr::handlers("progress")
      iter <- 1:private$k_split
      progressr::with_progress(enable = private$verbose,{
        #progress bar
        pb <- progressr::progressor(along = iter)
        #parallelization with future.apply
        fitted <- private$.f_lapply(
          iter=iter,
          func=function(i,...){
            #check whether to split samples
            if (private$k_split==1){
              train_index <- validation_index <- k_index==i
            } else{
              train_index <- k_index!=i
              validation_index <- k_index==i
            }

            num_val_index <- which(validation_index)
            names(num_val_index) <- rep(i,length(num_val_index))
            #split the sample based on the index
            #Q outcome set
            train_set.Q <- private$Q.set[train_index,]
            validation_set.Q <- private$Q.set[validation_index,]
            #g exposure set
            train_set.g <- data.frame(private$g.set[train_index,])
            validation_set.g <- data.frame(private$g.set[validation_index,])
            colnames(train_set.g)=colnames(validation_set.g)=colnames(private$g.set) #make to g df colnames consistent

            #Q model(outcome model: g-comp)
            #fit with train set
            Q.fit <- self$sl.fit(Y = private$Y[train_index],
                                           X = train_set.Q,
                                           SL.library = self$libs$Q.SL.library)
            # predict on validation set
            mu0 <- self$sl.predict(Q.fit,newdata=transform(validation_set.Q, A = 0)) #Q0_pred
            mu1  <- self$sl.predict(Q.fit,newdata=transform(validation_set.Q, A = 1)) #Q1_pred

            #g model(exposure model: propensity score)
            # fit with train set
            g.fit <- self$sl.fit(Y=private$A[train_index],
                                           X=train_set.g,
                                           SL.library = self$libs$g.SL.library)
            # predict on validation set
            raw_p_score  <- self$sl.predict(g.fit,newdata = validation_set.g)  #g_pred

            pb(sprintf("Iteration=%g/%g", i,private$k_split))
            output <- list(num_val_index,Q.fit,mu0,mu1,g.fit,raw_p_score)
            names(output) <- c("num_val_index","Q.fit","mu0","mu1","g.fit","raw_p_score")
            return(output)
          })
        })

      #store fitted values from future to member variables
      self$libs$num_val_index <- unlist(lapply(fitted,function(x) x$num_val_index))
      self$libs$Q.fit <- lapply(fitted,function(x) x$Q.fit)
      self$libs$g.fit <- lapply(fitted,function(x) x$g.fit)
      self$obs_est$mu0[self$libs$num_val_index] <- unlist(lapply(fitted,function(x) x$mu0))
      self$obs_est$mu1[self$libs$num_val_index] <- unlist(lapply(fitted,function(x) x$mu1))
      self$obs_est$raw_p_score[self$libs$num_val_index] <- unlist(lapply(fitted,function(x) x$raw_p_score))

      self$obs_est$mu  <- (self$obs_est$mu0*(1-private$A) + self$obs_est$mu1*(private$A)) #Q_pred

      if (private$verbose){
        cat("Done!\n")
      }

      invisible(self)
    }
  ),
  private = list(
    #input
    Q.set=NULL,
    g.set=NULL,
    k_split=NULL,
    verbose=NULL,
    g.bound=NULL,
    sl.pkg =NULL,
    #private methods
    #lapply or future_lapply
    .f_lapply =NULL,
  )
)
