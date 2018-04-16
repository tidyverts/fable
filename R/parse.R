parse_specials <- function(call = NULL, specials = NULL, xreg = TRUE){
  if(is.null(call)){ # No model specified
    return(list())
  }
  call <- enexpr(call)
  parsed <- traverse_call(!!call,
                g = ~ .x %>%
                      get_expr %>% # Extract call
                      as.list %>% # Split call components
                      .[-1] %>% # Drop function operator (as it is known to be "+")
                      map(expr), # Rebuild quosure for recursive map
                h = function(x){ # Base types
                  x <- get_expr(x)
                  if(!is_call(x) || !(call_name(x) %in% specials)){
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

  # Recursively combine list of exogenous regressors using "+" in-order
  # if(!is.null(parsed$xreg)){
  #   parsed$xreg <- list(
  #     call2("xreg",
  #           traverse_list(parsed$xreg, 
  #                         f = ~ call2("+", .x[[1]], .y),
  #                         g = ~ list(.x[-length(.x)]),
  #                         h = ~ .x[[length(.x)]],
  #                         base = ~ length(.x) <= 1)
  #           )
  #     )
  # }
  if(!is.null(parsed$xreg)){
    parsed$xreg <- list(call2("xreg", parsed$xreg))
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

#' @importFrom tibble tibble
#' @importFrom tsibble measured_vars
parse_model <- function(data, model, specials){
  # Clean inputs
  if(quo_is_missing(model)){
    model <- set_expr(model, sym(measured_vars(data)[1]))
    inform(sprintf(
      "Model not specified, using defaults set by `formula = %s`. Override this using `formula`.",
      quo_text(model)
      ))
  }
  if(is.null(model_lhs(model))){
    model <- set_expr(model, new_formula(lhs = sym(measured_vars(data)[1]), rhs = model_rhs(model)))
    inform(sprintf(
      "Response not specified, automatically selected `%s`. Override this in the `formula`.",
      expr_text(measured_vars(data)[1])
    ))
  }
  
  list(model = model) %>%
    append(parse_model_lhs(model_lhs(model), data)) %>%
    append(parse_model_rhs(model_rhs(model), specials))
}

parse_model_rhs <- function(model_rhs, specials){
  model_rhs %>%
    parse_specials(specials = names(specials)) %>%
    map(~ .x %>%
          map(
            function(special){
              eval_tidy(special, env = specials)
            }
          )
    ) %>%
    list(args = .)
}

parse_model_lhs <- function(expr, data){
  list(
    response = eval_tidy(expr(parse_response(!!expr)), data = data),
    transformation = as_transformation(expr, data=data)
  )
}