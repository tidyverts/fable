#' Forecasting time series
#'
#' \code{forecast} is a generic function for forecasting from time series or
#' time series models. The function invokes particular \emph{methods} which
#' depend on the class of the first argument.
#'
#' @param object a time series or time series model for which forecasts are
#' required
#' @param h Number of periods for forecasting
#' @param ... Additional arguments affecting the forecasts produced. 
#' @return An object of class "\code{forecast}".
#'
#' @author Rob J Hyndman
#' @keywords ts
#' @importFrom forecast forecast
#' @export
# forecast <- function(object, h, ...) UseMethod("forecast")

