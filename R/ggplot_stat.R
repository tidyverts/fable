#' @rdname geom_forecast
#' @export
StatForecast <- ggplot2::ggproto(
  "StatForecast", ggplot2::Stat,
  required_aes = c("x", "y"),
  
  compute_group = function(data, scales, params, showgap=TRUE, series=NULL,
                           model=ETS(y), fc.args = list(), level = c(80, 95), ...) {
    model <- enexpr(model)
    if(inherits(scales$x, "ScaleContinuousDatetime")){
      index <- as.POSIXct(data$x, origin = "1970-01-01")
    }
    else if(inherits(scales$x, "ScaleContinuousDate")){
      index <- as.Date(data$x, origin = "1970-01-01")
    }
    else if(inherits(scales$x, "ScaleContinuous")){
      index <- data$x
    }
    else{
      abort("Cannot determine time scale on x-axis. Perhaps produce your forecasts separately, then use `autolayer()`")
    }
    if(inherits(scales$x, "ScaleContinuousYearmonth")){
      index <- tsibble::yearmonth(index)
    }
    else if(inherits(scales$x, "ScaleContinuousYearquarter")){
      index <- tsibble::yearquarter(index)
    }
    else if(inherits(scales$x, "ScaleContinuousYearweek")){
      index <- tsibble::yearweek(index)
    }

    plot_data <- tsibble(x = index, y = data$y, index=x)

    fit <- eval_tidy(quo(plot_data %>% !!model))
    fcast <- do.call("forecast", append(list(fit), fc.args))
    fcast <- fortify(fcast, level = level, showgap = showgap) %>%
      rename(!!!exprs(y = !!sym("mean"), ymin = !!sym("lower"), ymax = !!sym("upper"))) %>%
      mutate(x := as.numeric(!!sym("x")))
    if (!is.null(series)) {
      if (data$group[1] > length(series)) {
        message("Recycling series argument, please provide a series name for each time series")
      }
      fcast <- transform(fcast, series = series[(abs(data$group[1]) - 1) %% length(series) + 1])
    }
    fcast
  }
)
