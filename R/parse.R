parse_specials <- function(call = NULL, specials = NULL, xreg = TRUE){
  if(!is.null(call)){ # No model specified
    call <- enexpr(call)
    parsed <- traverse_call(!!call,
                            g = ~ .x %>%
                              get_expr %>% # Extract call
                              as.list %>% # Split call components
                              .[-1] %>% # Drop function operator (as it is known to be "+")
                              map(expr), # Rebuild quosure for recursive map
                            h = function(x){ # Base types
                              x <- get_expr(x)
                              if(!is_call(x) || !(call_name(x) %in% names(specials))){
                                if(!xreg) stop("Exogenous regressors are not supported for this model type")
                                list(xreg = list(x))
                              }
                              else{# Current call is a special function
                                list(list(x)) %>% set_names(call_name(x))
                              }
                            },
                            f = function(.x, ...) {
                              merge_named_list(.x[[1]], .x[[2]])},
                            base = ~ .x %>% get_expr %>% {!is_call(.) || call_name(.) != "+"}
    )
  } else {
    parsed <- list()
  }
  

  # Add required_specials
  missing_specials <- attr(specials, "required_specials") %>% 
    .[!(.%in%names(parsed))]
  parsed <- parsed %>%
    append(
      missing_specials %>%
        map(~ list(call2(.x))) %>%
        set_names(missing_specials)
    )
  
  if(!is.null(parsed$xreg)){
    parsed$xreg <- list(call2("xreg", !!!parsed$xreg))
  }
  parsed
}

parse_response <- function(model_lhs){
  model_lhs <- enquo(model_lhs)
  
  # Traverse call along longest argument (hopefully, the response)
  traverse_call(!!model_lhs,
                f = ~ .x[[1]],
                g = ~ .x %>%
                  get_expr %>%
                  as.list %>% 
                  map(new_quosure, env = get_env(.x)) %>%
                  .[which.max(map(., ~ length(eval_tidy(.x))))],
                h = ~ .x) %>%
    get_expr
}

#' Validate the user provided model
#' 
#' Appropriately format the user's model for evaluation. Typically ran as one of the first steps
#' in a model function.
#' @param model The user's model specification (unevaluated)
#' @param data A dataset used for automatic response selection
#' 
#' @importFrom purrr compose
#' 
#' @export
validate_model <- function(model, data = NULL){
  # Capture parent expression
  model <- enquo(model)
  model_expr <- eval_tidy(set_expr(model, expr(rlang::enexpr(!!get_expr(model)))))
  # Clean inputs
  if(is_missing(model_expr)){
    model <- guess_response(data)
  }
  else{
    if(possibly(compose(is_formula, eval_tidy), FALSE)(model)){
      model <- eval_tidy(model)
      
      # Add response if missing
      if(is.null(model_lhs(model))){
        model <- new_formula(lhs = guess_response(data), rhs = model_rhs(model), env = get_env(model))
      }
    }
    else{
      model <- model_expr
    }
  }
  
  model
}

#' Parse the model specification for specials
#' 
#' Using a list of defined special functions, the user's formula specification and data
#' is parsed to extract important modelling components.
#' 
#' @param data A dataset
#' @param model A validated model (from `validate_model`)
#' @param specials A list of functions used to be evaluated from the formula.
#' If specials is NULL, no specials are computed
#' 
#' @importFrom tibble tibble
#' @export
parse_model <- function(data, model, specials = NULL){
  # Parse model
  list2(
    model = model,
    !!!parse_model_lhs(model_lhs(model), data),
    !!!parse_model_rhs(model_rhs(model), data, specials)
  )
}

#' Parse the RHS of the model formula for specials
#' 
#' @param model_rhs The expression for the rhs of the model (from `model_rhs()`)
#' @param data Data provided by the user to evaluate the special functions within
#' @param specials The environment containing specials (from `new_specials_env()`). 
#' If specials is NULL, no specials are computed
#' 
#' @export
parse_model_rhs <- function(model_rhs, data, specials = NULL){
  if(is.null(specials)){
    return(list())
  }
  model_rhs %>%
    parse_specials(specials = specials) %>%
    map(~ .x %>%
          map(
            function(special){
              eval_tidy(special, data = data, env = specials)
            }
          )
    ) %>%
    list(args = .)
}

#' Parse the RHS of the model formula for transformations
#' 
#' @param model_lhs The expression for the lhs of the model (from `model_lhs()`)
#' @param data Data to be used to find the response
#' 
#' @export
parse_model_lhs <- function(model_lhs, data){
  transformation <- as_transformation(model_lhs, data=data)
  response <- fn_fmls(transformation)$x
  list(
    response = response,
    transformation = transformation
  )
}