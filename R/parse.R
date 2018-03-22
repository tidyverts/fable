traverse_list <-  function(x, f = ~ .x, base = ~ !is.list(.x)){
  f <- as_mapper(f)
  base <- as_mapper(base)
  
  # base case
  if(base(x)){
    return(x)
  }
  
  # recursive case
  f(lapply(x, traverse_list))
}

quo_listcall_to_call <- as_mapper(~ map(.x, quo_get_expr) %>% as.call %>% new_quosure(env = get_env(.x[[1]])))
quo_call_to_listcall <- as_mapper(~ .x %>% rlang::get_expr(.) %>% as.list %>% map(new_quosure, env = get_env(.x)))
# quo_is_known <- as_mapper(.x %>% rlang::get_expr(.) %>% {rlang::is_syntactic_literal(.) || is_symbol(.) || (!is.pairlist(.) && !is.call(.))})

traverse_call <- function(x,
                          f = quo_listcall_to_call,
                          g = quo_call_to_listcall,
                          h = enquo,
                          base = ~ !quo_is_call(.x)){
  .f <- as_mapper(f)
  .g <- as_mapper(g)
  .h <- as_mapper(h)
  .base <- as_mapper(base)
  
  # base cases
  if (.base(.h(x)))
    return(.h(x))
  
  # recursive case
  .f(map(.g(.h(x)), traverse_call, f=f, g=g, base=base), .h(x))
}

# # Default call traversal
# traverse_call(a+b)
# 
# # Max length object finding
# traverse_call(log(USAccDeaths + 1),
#               f = ~ .x[[1]],
#               g = ~ .x %>% rlang::get_expr(.) %>%
#                 as.list %>% map(new_quosure, env = get_env(.x)) %>%
#                 .[which.max(map(., ~ length(eval_tidy(.x))))])
# 
# # Basic operation ordering
# traverse_call(log(USAccDeaths + 1),
#               f = ~ .x[[1]],
#               g = ~ .x %>% rlang::get_expr(.) %>%
#                 as.list %>% map(new_quosure, env = get_env(.x)) %>%
#                 .[which.max(map(., ~ length(eval_tidy(.x))))],
#               base = )
