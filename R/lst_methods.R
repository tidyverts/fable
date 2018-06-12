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

#' @export
type_sum.lst_fc <- function(x){
  "fc"
}

#' @export
pillar_shaft.lst_fc <- function(x, ...){
  pillar::new_pillar_shaft_simple(format(x))
}

#' @export
format.lst_fc <- function(x, ...){
  x %>% map_chr(function(x){
    dist <- attr(x$distribution, "qname")
    h <- NROW(x)
    if(attr(x$distribution, "trans")){
      dist <- sprintf("t(%s)", dist)
    }
    sprintf("~%s [h=%i]", dist, h)
  })
}

#' @export
c.lst_fc <- function(x, ...){
  enclass(NextMethod(), "lst_fc")
}

#' @export
`[.lst_fc` <- c.lst_fc

#' @export
print.lst_fc <- function(x, ...){
  class(x) <- "list"
  print(x)
}
