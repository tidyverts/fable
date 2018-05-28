#' @importFrom stats residuals
#' @export
residuals.mable <- function(object, ...){
  object %>%
    transmute(!!!key(.),
              residuals = map2(!!sym("data"), !!sym("model"),
                               function(data, model) {
                                 data %>% transmute(residuals = data[[expr_text(attr(model, "response"))]] - as.numeric(fitted(model)))
                               }
              )
    )%>%
    unnest(key = id(!!!key(object)))
}

#' @importFrom stats fitted
#' @export
fitted.mable <- function(object, ...){
  object %>%
    transmute(!!!key(.),
              fitted = map2(!!sym("data"), !!sym("model"),
                             function(data, model) {
                               data %>% transmute(fitted = as.numeric(fitted(model)))
                             }
              )
    ) %>%
    unnest(key = id(!!!key(object)))
}

key.mable <- function(x){
  syms(setdiff(colnames(x), c("data", "model")))
}