#' @importFrom ggplot2 ggplot aes geom_line guides guide_legend
#' @export
autoplot.tbl_ts <- function(object, var = sym(measured_vars(object)[1]), ...){
  if(!missing(var)){
    var <- enexpr(var)
  }
  else if(length(measured_vars(object)) > 1){
    inform(sprintf(
      "Plot variable not specified, automatically selected `var = %s`",
      measured_vars(object)[1]
    ))
  }
  
  aes_spec <- list(x = index(object), y = var)
  if(n_keys(object) > 1){
    aes_spec["colour"] <- list(expr(interaction(!!!key(object), sep = "/")))
  }
  
  ggplot(object, eval_tidy(expr(aes(!!!aes_spec)))) + 
    geom_line() + 
    guides(colour = guide_legend(paste0(map(key(object), expr_text), collapse = "/")))
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
  autoplot(suppressMessages(as_tsibble(unnest(object, !!sym("data"))))) +
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
                  mean = !!attr(model, "response"),
                  !!!set_names(map(level, ~ expr(new_hilo(mean, mean, !!.x))), level)
        )
    }
    gap <- suppressWarnings(object %>% 
      transmute(
        !!!syms(key_vars(.)),
        gap = map2(data,model, extract_last_obs)
      ) %>%
      unnest(gap, key = id(!!!key_vars(object)))
    )
  }
  else{
    gap <- NULL
  }
  tsbl <- suppressWarnings(
    object %>% 
      select(!!!key_vars(object), forecast) %>%
      mutate(
        forecast = map(forecast, 
                       function(fc){
                         fc %>%
                           mutate(!!!set_names(map(level, ~ expr(hilo(!!sym("quantile"), !!.x))), level)) %>%
                           select(exclude("quantile")) %>%
                           mutate(mean = as.numeric(mean)) %>%
                           rbind(gap)
                       })
      ) %>%
      unnest(forecast, key = id(!!!key_vars(object))) %>%
      gather(level, hilo, -(!!index(.)), -mean) %>%
      mutate(hilo = enclass(hilo, "hilo"),
             level = level(hilo),
             lower = lower(hilo),
             upper = upper(hilo)) %>%
      select(exclude("hilo"))
  )
}

