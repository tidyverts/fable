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
                dist_name = expr_text(q_expr),
                trans = is.name(body(t_fn)))
}

#' @export
type_sum.quantile <- function(x){
  paste0(
    ifelse(x%@%"trans", "", "Transformed "),
    switch(x%@%"dist_name",
      qnorm = "Normal",
      x%@%"dist_name"
    )
  )
}