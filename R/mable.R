#' @importFrom stats residuals
#' @export
residuals.mable <- function(object, ...){
  object %>%
    transmute(!!!key(.),
              residuals = map2(data, model,
                               function(data, model) {
                                 data %>% transmute(residuals = data[[expr_text(attr(model, "response"))]] - as.numeric(fitted(model)))
                               }
              )
    )%>%
    unnest(.with = id(!!!key(object)))
}

#' @importFrom stats fitted
#' @export
fitted.mable <- function(object, ...){
  object %>%
    transmute(!!!key(.),
              fitted = map2(data, model,
                             function(data, model) {
                               data %>% transmute(fitted = as.numeric(fitted(model)))
                             }
              )
    ) %>%
    unnest(.with = id(!!!key(object)))
}

key.mable <- function(x){
  syms(setdiff(colnames(x), c("data", "model")))
}