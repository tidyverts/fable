# Create evaluation environment for specials
new_specials_env <- function(..., parent_env = caller_env()){
  child_env(parent_env, !!!dots_list(...))
}

specials_xreg <- function(x){
  list(xreg = expr(tibble(!!!x)))
}