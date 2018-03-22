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

traverse_call <- function(x,
                          f = ~ map(.x, quo_get_expr) %>% as.call %>% new_quosure(env = get_env(x)),
                          g = ~ .x %>% rlang::get_expr(.) %>% as.list %>% map(new_quosure, env = get_env(x)),
                          base = ~ .x %>% rlang::get_expr(.) %>%
                            {rlang::is_syntactic_literal(.) || is_symbol(.) || (!is.pairlist(.) && !is.call(.))}){
  x <- enquo(x)
  
  .f <- as_mapper(f)
  .g <- as_mapper(g)
  .base <- as_mapper(base)
  
  # base cases
  if (.base(x))
    return(x)
  
  # recursive case
  .f(map(.g(x), traverse_call, f=f, g=g, base=base))
}

# Basic call traversal
traverse_call(a+b)

# Basic max length object finding
traverse_call(log(USAccDeaths + 1),
              f = ~ .x[[1]],
              g = ~ .x %>% rlang::get_expr(.) %>%
                as.list %>% map(new_quosure, env = get_env(.x)) %>%
                .[which.max(map(., ~ length(eval_tidy(.x))))])
