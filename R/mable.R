#' Create a new mable
#' 
#' @param key_vals A list of values for keys
#' @param data The nested data (one row per model)
#' @param model A list of models
#'
#' @export
mable <- function(key_vals, data, model){
  new_tibble(tibble(!!!key_vals, data=data, model=enclass(model, "lst_mdl")), subclass = c("mable", "lst_ts"))
}

#' Coerce a dataset to a mable
#' 
#' @param data A dataset containing a list model column
#' @param model A bare input containing the model column's name
#' 
#' @export
as_mable <- function(data, model){
  model <- enexpr(model)
  data %>%
    mutate(!!!list(model = expr(enclass(!!model, "lst_mdl")))) %>%
    enclass("mable")
}

#' @importFrom tibble tbl_sum
#' @importFrom tsibble key_sum
#' @export
tbl_sum.mable <- function(x){
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
  
  out <- c(`A mable` = sprintf("%s model%s [%s]", big_mark(NROW(x)), ifelse(NROW(x)==1, "", "s"), int_disp))
  
  if(!is_empty(key_vars(x))){
    out <- c(out, key_sum(x))
  }
  
  out
}

#' Provide a succinct summary of a model
#' 
#' Similarly to pillar's type_sum and obj_sum, model_sum is used to provide brief model summaries.
#' 
#' @param x The model to summarise
#' 
#' @export
model_sum <- function(x){
  UseMethod("model_sum")
}

#' @export
#' @importFrom forecast forecast
#' @importFrom dplyr mutate
forecast.mable <- function(object, ...){
  fable(
    key_vals=as.list(object)[key_vars(object)],
    data=object$data, 
    model=object$model, 
    forecast=map2(object$model, object$data, forecast, ...)
  )
}

#' @importFrom stats residuals
#' @export
residuals.mable <- function(object, ...){
  object %>%
    transmute(!!!syms(key_vars(.)),
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
    transmute(!!!syms(key_vars(.)),
              fitted = map2(!!sym("data"), !!sym("model"),
                             function(data, model) {
                               data %>% transmute(fitted = as.numeric(fitted(model)))
                             }
              )
    ) %>%
    unnest(key = syms(key_vars(object)))
}

#' @importFrom forecast getResponse
#' @export
getResponse.mable <- function(object, ...){
  object %>%
    transmute(!!!syms(key_vars(.)),
              fitted = map2(!!sym("data"), !!sym("model"),
                            function(data, model) {
                              response <- model%@%"response"
                              data %>% transmute(response = data[[expr_text(attr(model, "response"))]])
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

#' @export
components.mable <- function(object, ...){
  object %>%
    transmute(
      !!!syms(key_vars(.)),
      components = map(object$model, components)
    ) %>%
    unnest(key = syms(key_vars(object)))
}

#' @export
key.mable <- key.dable

#' @export
key_vars.mable <- function(x){
  setdiff(colnames(x), c("data", "model"))
}

#' @export
n_keys.mable <- function (x){
  key <- key_vars(x)
  if (is_empty(key)) {
    return(1L)
  }
  NROW(distinct(ungroup(as_tibble(x)), !!!syms(key)))
}
