#' @importFrom ggplot2 ggplot aes geom_line guides guide_legend xlab
#' @export
autoplot.tbl_ts <- function(object, var = NULL, ...){
  if(quo_is_null(enquo(var))){
    inform(sprintf(
      "Plot variable not specified, automatically selected `var = %s`",
      measured_vars(object)[1]
    ))
    var <- sym(measured_vars(object)[1])
  }
  else{
    var <- enexpr(var)
  }
  
  aes_spec <- list(x = index(object), y = var)
  if(n_keys(object) > 1){
    aes_spec["colour"] <- list(expr(interaction(!!!syms(key_vars(object)), sep = "/")))
  }
  ggplot(object, eval_tidy(expr(aes(!!!aes_spec)))) + 
    geom_line() +
    guides(colour = guide_legend(paste0(map(syms(key_vars(object)), expr_text), collapse = "/"))) + 
    xlab(paste0(expr_text(index(object)), " [", format(interval(object)), "]"))
}

#' @export
autoplot.mable <- function(object, ...){
  if(length(object$model) > 1){
    inform("Only univariate models are supported at the moment, plotting the first model.")
  }
  autoplot(object$model[[1]])
}

# Add multiple via facets by rows
#' @export
autoplot.fable <- function(object, level = c(80, 95), ...){
  if(NROW(object)>1){
    warn("Only univariate forecast plots are currently supported. Plotting the first forecast.")
  }
  suppressMessages(autoplot(object$data[[1]], !!(object$model[[1]]%@%"response"))) +
    autolayer(object, level = level, ...)
}

#' @export
autolayer.fable <- function(object, level = c(80, 95), ...){
  data <- fortify(object, level = level, ...)
  mapping <- eval_tidy(quo(aes(x = !!index(data), y = !!sym("mean"))))
  if(!is.null(level)){
    mapping$level <- sym("level")
    mapping$ymin <- sym("lower")
    mapping$ymax <- sym("upper")
  }
  geom_forecast(mapping = mapping, stat = "identity", data = data)
}

#' @importFrom ggplot2 fortify
#' @export
fortify.fable <- function(object, level = c(80, 95), showgap = TRUE){
  # Tidy format with repeated predicted values
  if(!showgap){
    extract_last_obs <- function(data, model) {
      data %>%
        filter(!!index(.) == last(!!index(.))) %>%
        transmute(!!!syms(key_vars(.)),
                  !!index(.),
                  mean = !!attr(!!sym("model"), "response"),
                  !!!set_names(map(level, ~ expr(new_hilo(mean, mean, !!.x))), level)
        )
    }
    gap <- suppressWarnings(object %>% 
      transmute(
        !!!syms(key_vars(.)),
        gap = map2(!!sym("data"), !!sym("model"), extract_last_obs)
      ) %>%
      unnest(gap, key = id(!!!syms(key_vars(object))))
    )
  }
  else{
    gap <- NULL
  }
  tsbl <- suppressWarnings(
    object %>% 
      select(!!!syms(key_vars(object)), forecast) %>%
      mutate(
        forecast = map(forecast, 
                       function(fc){
                         fc %>%
                           mutate(!!!set_names(map(level, ~ expr(hilo(!!sym("distribution"), !!.x))), level)) %>%
                           select(exclude("distribution")) %>%
                           rbind(gap)
                       })
      ) %>%
      unnest(forecast, key = syms(key_vars(object)))
  )
  if(!is.null(level)){
    tsbl <- tsbl %>%
      gather(level, hilo, -(!!index(.)), -mean) %>%
      mutate(hilo = enclass(hilo, "hilo"),
             level = level(hilo),
             lower = lower(hilo),
             upper = upper(hilo)) %>%
      select(exclude("hilo"))
  }
  tsbl
}

