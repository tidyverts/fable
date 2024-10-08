% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/var.R
\name{generate.VAR}
\alias{generate.VAR}
\title{Generate new data from a fable model}
\usage{
\method{generate}{VAR}(x, new_data, specials, ...)
}
\arguments{
\item{x}{A fitted model.}

\item{new_data}{A tsibble containing the time points and exogenous regressors to produce forecasts for.}

\item{specials}{(passed by \code{\link[fabletools:forecast]{fabletools::forecast.mdl_df()}}).}

\item{...}{Other arguments passed to methods}
}
\description{
Simulates future paths from a dataset using a fitted model. Innovations are
sampled by the model's assumed error distribution. If \code{bootstrap} is \code{TRUE},
innovations will be sampled from the model's residuals. If \code{new_data}
contains the \code{.innov} column, those values will be treated as innovations.
}
\examples{
as_tsibble(USAccDeaths) \%>\%
  model(ETS(log(value) ~ season("A"))) \%>\%
  generate(times = 100)
}
\seealso{
\code{\link[fabletools:generate.mdl_df]{fabletools::generate.mdl_df}}
}
