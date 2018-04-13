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
