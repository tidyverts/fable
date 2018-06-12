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

#' @importFrom tibble tbl_sum
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
  
  out <- c(`A mable` = sprintf("%s models [%s]", big_mark(NROW(x)), int_disp))
  
  if(!is_empty(key_vars(x))){
    nk <- big_mark(n_keys(x))
    out <- c(out, Keys = sprintf("%s [%s]", paste0(key_vars(x), collapse = ", "), nk))
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

#' @export
summary.mable <- function(object, ...){
  map(object$model, ~capture.output(summary(.x))) %>%
    invoke(cat, ., sep="\n")
  invisible(object)
}

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
