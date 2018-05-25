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
#' @importFrom ggplot2 fortify
#' @export
fortify.fable <- function(object, level = c(80, 95)){
  # Tidy format with repeated predicted values
  suppressWarnings(
    object %>% 
      select(!!!key(object), forecast) %>%
      mutate(
        forecast = map(forecast, 
                       function(fc){
                         fc %>%
                           mutate(!!!set_names(map(level, ~ expr(hilo(!!sym("quantile"), !!.x))), level)) %>%
                           select(exclude("quantile"))
                       })
      ) %>%
      enclass("lst_ts") %>%
      unnest(forecast) %>% #, .with = id(!!!key(object))) %>%
      gather(level, hilo, -(!!index(.)), -mean) %>%
      mutate(hilo = enclass(hilo, "hilo"),
             lower = lower(hilo),
             upper = upper(hilo)) %>%
      select(exclude("hilo"))
  )
}

