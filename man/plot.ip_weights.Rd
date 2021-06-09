% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AIPW_base.R
\name{plot.ip_weights}
\alias{plot.ip_weights}
\title{Plot the inverse probability weights using truncated propensity scores by exposure status}
\value{
\code{ip_weights.plot} (public variable): A box plot of inverse probability weights using truncated propensity scores by exposure status (\code{ggplot2::geom_boxplot})
}
\description{
Plot and check the balance of propensity scores by exposure status
}
\section{R6 Usage}{

\code{$plot.ip_weights()}
}

\examples{
library(SuperLearner)
library(ggplot2)
aipw_sl <- AIPW$new(Y=rbinom(100,1,0.5), A=rbinom(100,1,0.5),
                    W.Q=rbinom(100,1,0.5), W.g=rbinom(100,1,0.5),
                    Q.SL.library="SL.mean",g.SL.library="SL.mean",
                    k_split=1,verbose=FALSE)$fit()
#after average treatment effect calculation calculation
aipw_sl$summary(g.bound=0.025)$plot.ip_weights()
}
\seealso{
\link{AIPW} and \link{AIPW_tmle}
}