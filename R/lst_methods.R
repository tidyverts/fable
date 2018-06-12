#' @export
type_sum.lst_mdl <- function(x){
  "model"
}

#' @export
pillar_shaft.lst_mdl <- function(x, ...){
  pillar::new_pillar_shaft_simple(format(x))
}

#' @export
format.lst_mdl <- function(x, ...){
  x %>% map_chr(model_sum)
}

#' @export
c.lst_mdl <- function(x, ...){
  enclass(NextMethod(), "lst_mdl")
}

#' @export
`[.lst_mdl` <- c.lst_mdl

#' @export
print.lst_mdl <- function(x, ...){
  class(x) <- "list"
  print(x)
}

