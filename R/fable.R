#' Forecasting Functions for Tidy Time Series
#' 
#' Methods and tools for displaying and analysing univariate time series
#' forecasts in a tidy format including exponential smoothing via state space models and
#' automatic ARIMA modelling.
#'
#' \tabular{ll}{ Package: \tab fable\cr Type: \tab Package\cr License: \tab
#' GPL3\cr LazyLoad: \tab yes\cr }
#'
#' @docType package
#' @name fable-package
#' @author Rob J Hyndman
#'
#' Maintainer: Rob.Hyndman@monash.edu
#' @keywords package
NULL # Instead of "_PACKAGE" to remove inclusion of \alias{forecast}
# "_PACKAGE"

globalVariables(".")

#' @import rlang
#' @import tsibble
#' @import tsibblestats
#' @importFrom purrr map map2 map_lgl map_chr map_dbl
NULL

#' Create a new fable
#' 
#' @inheritParams mable
#' @param forecast A list of tsibble forecasts
#' @export
fable <- function(key_vals, data, model, forecast){
  new_tibble(tibble(!!!key_vals,
                    data=data,
                    model=enclass(model, "lst_mdl"),
                    forecast=enclass(forecast, "lst_fc")),
             subclass = c("fable", "lst_ts"))
}

#' Coerce a dataset to a fable
#' 
#' @inheritParams as_mable
#' @param forecast A bare input containing the forecast column's name
#' 
#' @export
as_fable <- function(data, model, forecast){
  model <- enexpr(model)
  forecast <- enexpr(model)
  data %>%
    mutate(!!!list(model = expr(enclass(!!model, "lst_mdl")),
                   forecast = expr(enclass(!!forecast, "lst_mdl")))) %>%
    enclass("fable")
}

#' @importFrom tibble tbl_sum
#' @export
tbl_sum.fable <- function(x){
  intervals <- x %>%
    pull(!!sym("data")) %>%
    map(interval) %>%
    unique
  if(length(intervals)==1){
    int_disp <- format(intervals[[1]])
  }
  else{
    int_disp <- "MIXED"
  }
  
  out <- c(`A fable` = sprintf("%s forecasts [%s]", big_mark(NROW(x)), int_disp))
  
  if(!is_empty(key_vars(x))){
    nk <- big_mark(n_keys(x))
    out <- c(out, Keys = sprintf("%s [%s]", paste0(key_vars(x), collapse = ", "), nk))
  }
  
  out
}

#' @importFrom dplyr mutate_if
#' @export
summary.fable <- function(object, level=c(80,95), ...){
  suppressWarnings(
    object %>% 
      select(!!!syms(key_vars(object)), "forecast") %>%
      mutate(
        forecast = map(forecast, 
                       function(fc){
                         fc %>%
                           mutate(!!!set_names(map(level, ~ expr(hilo(!!sym("distribution"), !!.x))), paste0(level, "%"))) %>%
                           select(exclude("distribution"))
                       }
        )
      ) %>%
      unnest(forecast, key = syms(key_vars(object))) %>%
      mutate_if(is.list, enclass, "hilo")
    )
}

key_vars.fable <- function(x){
  setdiff(colnames(x), c("data", "model", "forecast"))
}

#' @export
n_keys.fable <- function (x){
  key <- key_vars(x)
  if (is_empty(key)) {
    return(1L)
  }
  NROW(distinct(ungroup(as_tibble(x)), !!!syms(key)))
}