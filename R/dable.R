#' Create a new dable
#' 
#' @param key_vals A list of values for keys
#' @param data The nested data (one row per decomposition)
#' @param decomposition A list of decompositions
#'
#' @export
dable <- function(key_vals, data, decomposition){
  new_tibble(tibble(!!!key_vals, data=data, decomposition=enclass(decomposition, "lst_dcmp")), subclass = c("dable", "lst_ts"))
}

#' Coerce a dataset to a dable
#' 
#' @param data A dataset containing a list decomposition column
#' @param decomposition A bare input containing the decomposition column's name
#' 
#' @export
as_dable <- function(data, decomposition){
  decomposition <- enexpr(decomposition)
  data %>%
    mutate(!!!list(decomposition = expr(enclass(!!decomposition, "lst_mdl")))) %>%
    enclass("mable")
}

#' @importFrom tibble tbl_sum
#' @export
tbl_sum.dable <- function(x){
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
  
  out <- c(`A dable` = sprintf("[%s]", int_disp))
  
  if(!is_empty(key_vars(x))){
    nk <- big_mark(n_keys(x))
    out <- c(out, Keys = sprintf("%s [%s]", paste0(key_vars(x), collapse = ", "), nk))
  }
  
  out
}

#' @export
components.dable <- function(object, ...){
  object %>%
    transmute(
      !!!syms(key_vars(.)),
      components = map(object$decomposition, components)
    ) %>%
    unnest(key = syms(key_vars(object)))
}

#' @export
key_vars.dable <- function(x){
  setdiff(colnames(x), c("data", "decomposition"))
}

#' @export
n_keys.dable <- function (x){
  key <- key_vars(x)
  if (is_empty(key)) {
    return(1L)
  }
  NROW(distinct(ungroup(as_tibble(x)), !!!syms(key)))
}