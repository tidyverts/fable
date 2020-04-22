% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/croston.R
\name{forecast.croston}
\alias{forecast.croston}
\title{Forecast a model from the fable package}
\usage{
\method{forecast}{croston}(object, new_data, specials = NULL, ...)
}
\arguments{
\item{object}{The time series model used to produce the forecasts}

\item{new_data}{A \code{tsibble} containing future information used to forecast.}

\item{specials}{(passed by \code{\link[fabletools:forecast.mdl_df]{fabletools::forecast.mdl_df()}}).}

\item{...}{Additional arguments for forecast model methods.}
}
\value{
A list of forecasts.
}
\description{
Produces forecasts from a trained model.
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
  forecast()
}