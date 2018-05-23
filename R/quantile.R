#' Create a quantile object
#'  
#' @param f A quantile function
#' @param ... Arguments for the quantile function
#' @param transformation Transformation to be applied to quantiles
#' @param abbr Abbreviation for display purposes, defaults to the name of the quantile function
#' 
#' @examples 
#' qt <- new_quantile(qnorm, mean = rep(3, 10), sd = seq(0, 1, length.out=10),
#'  transformation = exp, abbr = "N")
#' qt
#' hilo(qt, 95)
#' @export
new_quantile <- function(f, ..., transformation = ~ .x, abbr = NULL){
  f_quo <- enquo(f)
  t_fn <- as_mapper(transformation)
  pmap(dots_list(...), list) %>%
    enclass("quantile",
            f = f,
            t = t_fn,
            qname = abbr%||%quo_text(f_quo),
            trans = !is.name(body(t_fn)))
}

#' @export
type_sum.quantile <- function(x){
  "qt"
}

#' @export
obj_sum.quantile <- function(x) {
  format(x)
}

#' @export
print.quantile <- function(x, ...) {
  print(format(x, ...), quote = FALSE)
  invisible(x)
}

#' @export
format.quantile <- function(x, ...){
  x %>%
    map_chr(function(qt){
      args <- qt %>%
        imap(~ paste0(ifelse(nchar(.y)>0, paste0(.y, " = "), ""), format(.x, trim = TRUE, digits = 2))) %>%
        invoke("paste", ., sep = ", ")
      out <- paste0(
        attr(x, "qname"),
        "(", args, ")"
      )
      if(attr(x, "trans")){
        paste0("t(", out, ")")
      }
      else{
        out
      }
    })
}

#' @export
`[.quantile` <- function(x, ...){
  enclass(NextMethod(), "quantile", 
          !!!attributes(x))
}

#' @export
c.quantile <- function(...){
  sameAttr <- dots_list(...) %>%
    map(~ if(!inherits(.x, "quantile")) {abort("Only combining quantiles is supported")} else {attributes(.x)}) %>%
    duplicated %>%
    .[-1]
  if(any(!sameAttr)){
    abort("Cannot combine quantiles of different types.")
  }
    
  enclass(NextMethod(), "quantile", 
          !!!attributes(..1))
}

#' @export
length.quantile <- function(x){
  NextMethod()
}

# #' @importFrom ggplot2 aes_
# autoplot.quantile <- function(q_fn, q_range = c(0.0001, 0.9999), precision = 0.01){
#   tibble(x = seq(q_range[1], q_range[2], by = precision)) %>%
#     mutate(!!"density":=q_fn(!!sym("x"))) %>%
#     ggplot(aes_(x=~x, y=~density)) + 
#     geom_line()
# }
#' @export
hilo.quantile <- function(x, level = 95, ...){
  if(length(level)!=1){
    abort("Only one value of 'level' is supported.")
  }
  if (level < 0 || level > 100) {
    abort("'level' can't be negative or greater than 100.")
  }
  list(lower = 50-level/2, upper = 50+level/2) %>%
    map(function(level){
      eval_tidy(quo(attr(x, "t")(attr(x, "f")(level/100, !!!merge_named_list(!!!x))))) %>%
        unlist(recursive = FALSE, use.names = FALSE)
    }) %>%
    append(list(level = level)) %>%
    invoke("new_hilo", .)
}
