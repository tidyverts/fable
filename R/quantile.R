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
  
  list(list(
    f = f,
    t = t_fn,
    args = dots_list(...),
    qname = abbr%||%quo_text(f_quo),
    trans = !is.name(body(t_fn))
  )) %>%
    enclass("quantile")
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
  x %>%
    map(function(qt){
      args <- qt$args %>%
        imap(~ paste0(ifelse(nchar(.y)>0, paste0(.y, " = "), ""), format(.x, trim = TRUE, digits = 2))) %>%
        invoke("paste", ., sep = ", ")
      out <- paste0(
        qt$qname,
        "(", args, ")"
      )
      if(qt$trans){
        paste0("t(", out, ")")
      }
      else{
        out
      }
    }) %>%
    invoke("c", .)
}

#' @export
`[.quantile` <- function(x, i, ...){
  if(is_logical(i)){
    i <- which(i)
  }
  abs_i <- abs(i)
  sign_i <- sign(i)
  
  is_neg <- any(sign_i<=0)
  if(is_neg && any(sign_i>0)){
    abort("only 0's may be mixed with negative subscripts")
  }
  if(is_neg){
    sign_i <- -1
  }
  else{
    sign_i <- 1
  }
  
  qt_lens <- qt_lengths(x)
  offset <- cumsum(qt_lens)-qt_lens[1]
  x <- map2(qt_lens, offset,
       ~ sign_i*(intersect(abs_i,seq_len(.x) + .y) - .y)) %>%
    map2(x, function(i, qt){
      if(!is_neg || length(i)>0){
        qt$args <- qt$args %>%
          map(function(arg){arg[i, ...]})
      }
      qt
    })
  x[qt_lengths(x)>0] %>%
    enclass("quantile")
}

#' @export
c.quantile <- function(...){
  dots_list(...) %>%
    map(~ .x[[1]]) %>%
    enclass("quantile")
}

#' @export
length.quantile <- function(x){
  sum(qt_lengths(x))
}

qt_lengths <- function(x){
  x %>%
    map_dbl(function(qt){
      len <- qt$args %>% 
        map(length)
      if(any(len %>% map_lgl(~length(.x)>0))){
        len %>% invoke("max", .)
      }
      else{
        0
      }
    })
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
      x %>%
        map(~ eval_tidy(quo(.x$t(.x$f(level/100, !!!(.x$args)))))) %>%
        unlist(recursive = FALSE, use.names = FALSE)
    }) %>%
    append(list(level = level)) %>%
    invoke("new_hilo", .)
}
