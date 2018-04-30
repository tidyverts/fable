#' Create a quantile object
#'  
#' @param f A quantile function
#' @param ... Arguments for the quantile function
#' @param transformation Transformation to be applied to quantiles
#' @param abbr Abbreviation for display purposes, defaults to the name of the quantile function
#' 
#' @examples 
#' qt <- new_quantile(qnorm, mean = rep(100, 100), sd = 1:100, transformation = log, abbr = "N")
#' qt
#' qt(0.5)
#' @export
new_quantile <- function(f, ..., transformation = ~ .x, abbr = NULL){
  f_quo <- enquo(f)
  t_fn <- as_mapper(transformation)
  new_function(alist(level = c(80, 95)), expr({
    if(any(level > 1)){
      level <- level/100
    }
    eval_quantile(params, level)
  })) %>% 
    set_env(., child_env(environment(.),
                         params = list(list(
                           f = f,
                           t = t_fn,
                           args = dots_list(...),
                           qname = abbr%||%quo_text(f_quo),
                           trans = !is.name(body(t_fn))
                         )))) %>%
    enclass("quantile")
}

eval_quantile <- function(params, level){
  params %>% map(function(param){
    eval_tidy(quo(param$t(param$f(level, !!!param$args))))
  }) %>% 
    unlist
}

#' @export
type_sum.quantile <- function(x){
  format(x)
}

#' @export
print.quantile <- function(x, ...) {
  print(format(x, ...), quote = FALSE)
  invisible(x)
}

#' @export
format.quantile <- function(x, ...){
  environment(x)$params %>%
    map(function(param){
      args <- param$args %>%
        imap(~ paste0(ifelse(nchar(.y)>0, paste0(.y, " = "), ""), format(.x, trim = TRUE, digits = 2))) %>%
        invoke("paste", ., sep = ", ")
      out <- paste0(
        param$qname,
        "(", args, ")"
      )
      if(param$trans){
        paste0("t(", out, ")")
      }
      else{
        out
      }
    }) %>%
    invoke("c", .)
}

#' @export
`[.quantile` <- function(x, ...){
  environment(x)$args <- environment(x)$args %>%
    map(function(x) x[...])
  x
}

#' @export
c.quantile <- function(...){
  x <- dots_list(...)[[1]]
  environment(x)$params <- dots_list(...) %>% map(~ environment(.x)$params) %>% invoke("c", .)
  x
}

#' @export
length.quantile <- function(x){
  environment(x)$params %>%
    map(function(param){
      param %>%
        map(length) %>%
        invoke("max", .)
    }) %>%
    invoke("sum", .)
}

#' @importFrom ggplot2 aes_
autoplot.quantile <- function(q_fn, q_range = c(0.0001, 0.9999), precision = 0.01){
  tibble(x = seq(q_range[1], q_range[2], by = precision)) %>%
    mutate(!!"density":=q_fn(!!sym("x"))) %>%
    ggplot(aes_(x=~x, y=~density)) + 
    geom_line()
}
