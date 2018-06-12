#' Create a new mable
#' 
#' @param data The nested data (one row per model)
#' @param model A list of models
#'
#' @export
mable <- function(data, model){
  new_tibble(list(data=data, model=enclass(model, "lst_mdl")), subclass = "mable")
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
