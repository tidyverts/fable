#' @importFrom stats residuals
#' @export
residuals.mable <- function(object, ...){
  object %>%
    transmute(!!!key_vars(.),
              residuals = map2(!!sym("data"), !!sym("model"),
                               function(data, model) {
                                 data %>% transmute(residuals = data[[expr_text(attr(model, "response"))]] - as.numeric(fitted(model)))
                               }
              )
    )%>%
    unnest(key = syms(key_vars(object)))
}

#' @importFrom stats fitted
#' @export
fitted.mable <- function(object, ...){
  object %>%
    transmute(!!!key_vars(.),
              fitted = map2(!!sym("data"), !!sym("model"),
                             function(data, model) {
                               data %>% transmute(fitted = as.numeric(fitted(model)))
                             }
              )
    ) %>%
    unnest(key = syms(key_vars(object)))
}

#' @export
summary.mable <- function(object, ...){
  map(object$model, ~capture.output(summary(.x))) %>%
    invoke(cat, ., sep="\n")
  invisible(object)
}

key_vars.mable <- function(x){
  syms(setdiff(colnames(x), c("data", "model")))
}