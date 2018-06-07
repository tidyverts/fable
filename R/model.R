#' Extract the left hand side of a model
#' 
#' @param model A formula
#' 
#' @export
model_lhs <- function(model){
  if(is_formula(model)){
    transform_spec <- f_lhs(model)
  }
  else{
    transform_spec <- model
  }
}

#' Extract the right hand side of a model
#' 
#' @param model A formula
#' 
#' @export
model_rhs <- function(model){
  if(is_formula(model)){
    f_rhs(model)
  }
  else{
    expr(NULL)
  }
}