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
               validation_index = NULL),
    #' @field sl.fit A wrapper function for fitting SuperLearner or sl3
    sl.fit = NULL,
    #' @field sl.predict A wrapper using \code{sl.fit} to predict
    sl.predict = NULL,

    #' @description
    #' Create a new `AIPW` object.
    #'
    #' @param Y Outcome (binary integer: 0 or 1)
    #' @param A Exposure (binary integer: 0 or 1)
    #' @param verbose Whether to print the result (logical; Default = FALSE)
    #' @param W covariates for both exposure and outcome models  (vector, matrix or data.frame). If null, this function will seek for
    #' inputs from `W.Q` and `W.g`.
    #' @param W.Q Only valid when `W` is null, otherwise it would be replaced by `W`.
    #' Covariates for outcome model (vector, matrix or data.frame).
    #' @param W.g Only valid when `W` is null, otherwise it would be replaced by `W`.
    #' Covariates for exposure model (vector, matrix or data.frame)
    #' @param Q.SL.library SuperLearner libraries or sl3 learner object (Lrnr_base) for outcome model
    #' @param g.SL.library SuperLearner libraries or sl3 learner object (Lrnr_base) for exposure model
    #' @param k_split Number of splitting (integer; range: from 1 to number of observation-1):
    #'   if k_split=1, no sample splitting;
    #'   if k_split>1, use similar technique of cross-validation
    #'   (e.g., k_split=10, use 9/10 of the data to estimate and the remaining 1/10 leftover to predict)
    #'   NOTE: it's recommended to use sample splitting.
    #'
    #' @return A new `AIPW` object
    #'
    #' @examples
    #' library(SuperLearner)
    #' aipw_sl <- AIPW$new(Y=rbinom(100,1,0.5), A=rbinom(100,1,0.5),
    #'                     W.Q=rbinom(100,1,0.5), W.g=rbinom(100,1,0.5),
    #'                     Q.SL.library="SL.mean",g.SL.library="SL.mean",
    #'                     k_split=1,verbose=FALSE)
    initialize = function(Y=NULL, A=NULL, verbose=FALSE,
                          W=NULL, W.Q=NULL, W.g=NULL,
                          Q.SL.library=NULL, g.SL.library=NULL,
                          k_split=10){
      #-----initialize from AIPW_base class-----#
      super$initialize(Y=Y,A=A,verbose=verbose)
      #decide covariate set(s): W.Q and W.g only works when W is null.
      if (is.null(W)){
        if (any(is.null(W.Q),is.null(W.g))) {
          stop("No sufficient covariates were provided.")
        } else{
          private$Q.set=cbind(A, as.data.frame(W.Q))
          private$g.set=as.data.frame(W.g)
        }
      } else{
        private$Q.set=cbind(A, as.data.frame(W))
        private$g.set=as.data.frame(W)
      }
      #save input into private fields
      private$k_split=k_split
      #check data length
      if (!(length(private$Y)==dim(private$Q.set)[1] & length(private$A)==dim(private$g.set)[1])){
        stop("Please check the dimension of the data")
      }
      #-----determine SuperLearner or sl3 and change accordingly-----#
      if (is.character(Q.SL.library) & is.character(g.SL.library)) {
        if (any(grepl("SL.",Q.SL.library)) & any(grepl("SL.",g.SL.library))){
          #change future package loading
          private$sl.pkg <- "SuperLearner"
          #create a new local env for superlearner
          private$sl.env = new.env()
          #find the learners in global env and assign them into sl.env
          private$sl.learners = grep("SL.",lsf.str(globalenv()),value = T)
          lapply(private$sl.learners, function(x) assign(x=x,value=get(x,globalenv()),envir=private$sl.env))
          #change wrapper functions
          self$sl.fit = function(Y, X, SL.library, CV){
            suppressMessages({
              fit <- SuperLearner::SuperLearner(Y = Y, X = X, SL.library = SL.library, family="binomial",
                                                env=private$sl.env, cvControl = CV)
            })
            return(fit)
          }
          self$sl.predict = function(fit, newdata){
            suppressMessages({
              pred <- as.numeric(predict(fit,newdata = newdata)$pred)
            })
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
          self$sl.fit = function(X, Y, SL.library, CV){
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
      #------input checking-----#
      #check k_split value
      if (private$k_split<1 | private$k_split>=self$n){
        stop("`k_split` is not valid")
      } else if (private$k_split %in% 2){
        warning("One fold cross-validation will be used.")
      }
      #check verbose value
      if (!is.logical(private$verbose)){
        stop("`verbose` is not valid")
      }
      #check if SuperLearner and/or sl3 library is loaded
      if (!any(names(sessionInfo()$otherPkgs) %in% c("SuperLearner","sl3"))){
        warning("Either `SuperLearner` or `sl3` package is not loaded.")
      }
      #-------check if future.apply is loaded otherwise lapply would be used.------#
      if (any(names(sessionInfo()$otherPkgs) %in% c("future.apply"))){
        private$.f_lapply = function(iter,func) {
          future.apply::future_lapply(iter,func,future.seed = T,future.packages = private$sl.pkg,future.globals = TRUE)
        }
      }else{
        private$.f_lapply = function(iter,func) lapply(iter,func)
        }
    },
    #' @description
    #' Fitting the data into the `AIPW` object with/without sample splitting to estimate the influence functions
    #'
    #' @return A fitted `AIPW` object
    #'
    #' @examples
    #' library(SuperLearner)
    #' aipw_sl <- AIPW$new(Y=rbinom(100,1,0.5), A=rbinom(100,1,0.5),
    #'                     W.Q=rbinom(100,1,0.5), W.g=rbinom(100,1,0.5),
    #'                     Q.SL.library="SL.mean",g.SL.library="SL.mean",
    #'                     k_split=1,verbose=FALSE)
    #' aipw_sl$fit()
    fit = function(){
      #----------create index for sample splitting---------#
      private$cv$k_index <- sample(rep(1:private$k_split,ceiling(self$n/private$k_split))[1:self$n],replace = F)
      private$cv$fold_index = split(1:self$n, private$cv$k_index)
      private$cv$fold_length = sapply(private$cv$fold_index,length)
      iter <- 1:private$k_split

      #----------------progress bar setup----------#
      #check if progressr is loaded
      if (!any(names(sessionInfo()$otherPkgs) %in% c("progressr"))){
        private$isLoaded_progressr = TRUE
        pb <- progressr::progressor(along = iter)
      }

      #---------parallelization with future.apply------#
      fitted <- private$.f_lapply(
        iter=iter,
        func=function(i,...){
          #when k_split in 1:2, no cvControl will be used (same cv for k_split)
          if (private$k_split==1){
            train_index <- validation_index <- as.numeric(unlist(private$cv$fold_index))
            cv_param = list()
          } else if (private$k_split==2){
            train_index <- as.numeric(unlist(private$cv$fold_index[-i]))
            validation_index <- as.numeric(unlist(private$cv$fold_index[i]))
            cv_param = list()
          } else{
            train_index <- as.numeric(unlist(private$cv$fold_index[-i]))
            validation_index <- as.numeric(unlist(private$cv$fold_index[i]))
            cv_param = list(V=private$k_split-1,
                            validRows= private$.new_cv_index(val_fold=i))
          }

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
                               SL.library = self$libs$Q.SL.library,
                               CV= cv_param)
          # predict on validation set
          mu0 <- self$sl.predict(Q.fit,newdata=transform(validation_set.Q, A = 0)) #Q0_pred
          mu1  <- self$sl.predict(Q.fit,newdata=transform(validation_set.Q, A = 1)) #Q1_pred

          #g model(exposure model: propensity score)
          # fit with train set
          g.fit <- self$sl.fit(Y=private$A[train_index],
                               X=train_set.g,
                               SL.library = self$libs$g.SL.library,
                               CV= cv_param)
          # predict on validation set
          raw_p_score  <- self$sl.predict(g.fit,newdata = validation_set.g)  #g_pred

          #add metadata
          names(validation_index) <- rep(i,length(validation_index))

          if (private$isLoaded_progressr){
            pb(sprintf("No.%g iteration", i,private$k_split))
          }
          output <- list(validation_index,Q.fit,mu0,mu1,g.fit,raw_p_score)
          names(output) <- c("validation_index","Q.fit","mu0","mu1","g.fit","raw_p_score")
          return(output)
        })

      #store fitted values from future to member variables
      for (i in fitted){
        #add estimates based on the val index
        self$obs_est$mu0[i$validation_index] <- i$mu0
        self$obs_est$mu1[i$validation_index] <- i$mu1
        self$obs_est$raw_p_score[i$validation_index] <- i$raw_p_score
        #append fitted objects
        self$libs$Q.fit = append(self$libs$Q.fit, i$Q.fit)
        self$libs$g.fit = append(self$libs$g.fit, i$g.fit)
        self$libs$validation_index = append(self$libs$validation_index, i$validation_index)
      }
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
    cv = list(
      #a vector stores the groups for splitting
      k_index= NULL,
      #a list of indices for each fold
      fold_index= NULL,
      #a vector of length(fold_index[[i]])
      fold_length = NULL
    ),
    fitted=NULL,
    sl.pkg =NULL,
    sl.env=NULL,
    sl.learners = NULL,
    isLoaded_progressr = FALSE,
    #private methods
    #lapply or future_lapply
    .f_lapply =NULL,
    #create new index for training set
    .new_cv_index = function(val_fold,fold_length=private$cv$fold_length, k_split=private$k_split){
      train_fold_length = c(0,fold_length[-val_fold])
      train_fold_cumsum = cumsum(train_fold_length)
      new_train_index= lapply(1:(k_split-1),
                              function(x) {
                                (1:train_fold_length[[x+1]])+ train_fold_cumsum[[x]]
                              }
      )
      names(new_train_index) = names(train_fold_length[-1])
      return(new_train_index)
    }
  )
)
