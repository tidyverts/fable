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
  if(is.null(expr)){
    expr <- sym(measured_vars(data)[1])
    inform(sprintf(
      "Response not specified, automatically selected `%s`. Override this in the `formula`.",
      expr_text(expr)
    ))
  }
  transformation_stack <- eval_tidy(expr(traverse_transformation(!!expr)), data = data)
  list(
    response = get_expr(last(transformation_stack)),
    transform = new_function(eval_tidy(expr(alist(x = !!get_expr(last(transformation_stack))))), eval_tidy(expr(substitute(!!expr, set_names(list(sym("x")), expr_text(expr)))))),
    backtransform = new_function(eval_tidy(expr(alist(x = !!get_expr(last(transformation_stack))))), invert_transformation(transformation_stack))
  )
}