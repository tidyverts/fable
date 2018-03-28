traverse <-  function(x, f = ~ .x, g = ~.x, h = ~.x, base = ~ is_syntactic_literal(.x) || is_symbol(.x)){
  .f <- as_mapper(f)
  .g <- as_mapper(g)
  .h <- as_mapper(h)
  .base <- as_mapper(base)
  
  # base case
  if(.base(x))
    return(.h(x))
  
  # recursive case
  .f(map(.g(x), traverse, f=f, g=g, h=h, base=base), .h(x))
}

traverse_list <- function(x, f = ~ as.list(.x), g = ~.x, h = ~.x, base = ~ !is.list(.x)){
  traverse(x, f=f, g=g, h=h, base=base)
}

traverse_call <- function(x,
                          f = as_mapper(~ map(.x, quo_get_expr) %>% as.call %>% new_quosure(env = get_env(.x[[1]]))), # (quo_listcall_to_call)
                          g = as_mapper(~ .x %>% get_expr %>% as.list %>% map(new_quosure, env = get_env(.x))), # (quo_call_to_listcall)
                          h = ~ .x,
                          base = ~ !quo_is_call(.x)){
  x <- enquo(x)
  traverse(x, f=f, g=g, h=h, base=base)
}
