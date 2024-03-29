% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/croston.R
\name{fitted.croston}
\alias{fitted.croston}
\title{Extract fitted values from a fable model}
\usage{
\method{fitted}{croston}(object, ...)
}
\arguments{
\item{object}{A model for which forecasts are required.}

\item{...}{Other arguments passed to methods}
}
\value{
A vector of fitted values.
}
\description{
Extracts the fitted values.
}
\examples{
library(tsibble)
sim_poisson <- tsibble(
  time = yearmonth("2012 Dec") + seq_len(24),
  count = rpois(24, lambda = 0.3),
  index = time
)

sim_poisson \%>\%
  model(CROSTON(count)) \%>\%
  tidy()
}
