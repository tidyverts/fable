model_lhs <- function(model){
  model <- get_expr(model)
  if(is_formula(model)){
    transform_spec <- f_lhs(model)
  }
  else{
    transform_spec <- model
  }
}

model_rhs <- function(model){
  model <- get_expr(model)
  if(is_formula(model)){
    f_rhs(model)
  }
  else{
    expr(NULL)
  }
}