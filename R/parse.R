# Small function to combine two named lists
merge_named_list <- function(x,y){
  all_names <- union(names(x), names(y))
  all_names %>%
    map(~ c(x[[.x]], y[[.x]])) %>%
    set_names(all_names)
}

parse_specials <- function(call, specials = NULL, xreg = TRUE){
  call <- enexpr(call)
  parsed <- traverse_call(!!call,
                g = ~ .x %>%
                      get_expr %>% # Extract call
                      as.list %>% # Split call components
                      .[-1] %>% # Drop function operator (as it is known to be "+")
                      map(expr), # Rebuild quosure for recursive map
                h = function(x){ # Base types
                  if(!is_call(x) || !(call_name(x) %in% specials)){
                    if(!xreg) stop("Exogenous regressors are not supported for this model type")
                    list(xreg = get_expr(x))
                  }
                  else{# Current call is a special function
                    list(get_expr(x)) %>% set_names(call_name(x))
                  }
                },
                f = function(.x, ...) {
                  merge_named_list(.x[[1]], .x[[2]])},
                base = ~ !is_call(.x) || call_name(.x) != "+"
  )
  # Recursively combine list using "+" in-order
  parsed$xreg <- list(traverse_list(parsed$xreg, 
                f = ~ call2("+", .x[[1]], .y),
                g = ~ list(.x[-length(.x)]),
                h = ~ .x[[length(.x)]],
                base = ~ length(.x) <= 1))
  parsed
}