new_quantile <- function(f, transformation = ~ .x, ...){
  q_expr <- enexpr(f)
  t_fn <- as_mapper(transformation)
  new_function(alist(level = c(80, 95)), expr({
    if(any(level > 1)){
      level <- level/100
    }
    (!!t_fn)((!!q_expr)(level, !!!dots_list(...)))
  })
  ) %>% enclass("quantile",
                qname = expr_text(q_expr),
                trans = is.name(body(t_fn)))
}

#' @export
type_sum.quantile <- function(x){
  paste0(
    ifelse(x%@%"trans", "", "Transformed "),
    switch(x%@%"qname",
      qnorm = "Normal",
      x%@%"qname"
    )
  )
}

#' @importFrom ggplot2 aes_
autoplot.quantile <- function(q_fn, q_range = c(0.0001, 0.9999), precision = 0.01){
  tibble(x = seq(q_range[1], q_range[2], by = precision)) %>%
    mutate(!!"density":=q_fn(x)) %>%
    ggplot(aes_(x=~x, y=~density)) + 
    geom_line()
}