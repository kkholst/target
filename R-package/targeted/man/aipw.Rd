% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aipw.R
\name{aipw}
\alias{aipw}
\title{AIPW estimator}
\usage{
aipw(response_model, data, formula = ~1, missing_model, ...)
}
\arguments{
\item{response_model}{Model for the response given covariates (ml_model or
formula)}

\item{data}{data.frame}

\item{formula}{design specifying the OLS estimator with outcome given by the
EIF}

\item{missing_model}{Optional missing_model (ml_model or formula). By
default will use the same design as the response_model.}

\item{...}{arguments to cate}
}
\description{
AIPW for the mean (and linear projections of the EIF) with missing
observations
}
\examples{
m <- lvm(y ~ x+z, r ~ x)
distribution(m,~ r) <- binomial.lvm()
transform(m, y0~r+y) <- function(x) { x[x[,1]==0,2] <- NA; x[,2] }
d <- sim(m,1e3,seed=1)

aipw(y0 ~ x, data=d)
}
