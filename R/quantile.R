#' @export
new_quantile <- function(f, ..., transformation = ~ .x, abbr = NULL){
  f_quo <- enquo(f)
  t_fn <- as_mapper(transformation)
  new_function(alist(level = c(80, 95)), expr({
    if(any(level > 1)){
      level <- level/100
    }
    eval_quantile(f, t, level, args)
  })) %>% 
    set_env(., child_env(environment(.),
                         f = f,
                         t = t_fn,
                         args = dots_list(...))) %>%
    enclass("quantile",
            qname = abbr%||%quo_text(f_quo),
            trans = !is.name(body(t_fn)))
}

eval_quantile <- function(f, t, level, args){
  eval_tidy(quo(t(f(level, !!!args))))
}

#' @export
type_sum.quantile <- function(x){
  paste0(
    ifelse(x%@%"trans", "", "Transformed "),
    x%@%"qname"
  )
}

#' @export
print.quantile <- function(x, ...) {
  print(format(x, ...), quote = FALSE)
  invisible(x)
}

#' @export
format.quantile <- function(x, ...){
  args <- environment(x)$args %>%
    imap(~ paste0(ifelse(nchar(.y)>0, paste0(.y, " = "), ""), format(.x, trim = TRUE, digits = 2))) %>%
    invoke(paste, ., sep = ", ")
  out <- paste0(
    x%@%"qname",
    "(", args, ")"
  )
  if(x%@%"trans"){
    paste0("t(", out, ")")
  }
  else{
    out
  }
}

length.quantile <- function(x){
  environment(x)$args %>%
    map(length) %>%
    invoke(max, .)
}

#' @importFrom ggplot2 aes_
autoplot.quantile <- function(q_fn, q_range = c(0.0001, 0.9999), precision = 0.01){
  tibble(x = seq(q_range[1], q_range[2], by = precision)) %>%
    mutate(!!"density":=q_fn(x)) %>%
    ggplot(aes_(x=~x, y=~density)) + 
    geom_line()
}
