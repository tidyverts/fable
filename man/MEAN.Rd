% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mean.R
\name{MEAN}
\alias{MEAN}
\alias{report.model_mean}
\title{Mean models}
\usage{
MEAN(formula, ...)
}
\arguments{
\item{formula}{Model specification.}

\item{...}{Not used.}
}
\description{
\code{MEAN()} returns an iid model applied to the formula's response variable.
}
\details{
The model does not use any specials, and so everything on the formula's
right-hand-side will be ignored.
}
\section{Specials}{


This model does not support usage of any specials. It only computes the mean!
}

\examples{
library(tsibbledata)
vic_elec \%>\% 
  model(rw = MEAN(Demand))

}
\seealso{
\href{https://otexts.com/fpp3/simple-methods.html}{Forecasting: Principles and Practices, Some simple forecasting methods (section 3.2)}
}