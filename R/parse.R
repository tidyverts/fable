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
parse_model <- function(data, model, specials, univariate = TRUE){
  # Capture user call
  cl <- call_standardise(sys.call(-1))
  
  # Coerce data
  data <- as_tsibble(data)

  # Handle multivariate inputs
  if(univariate){
    if(n_keys(data) > 1){
      return(multi_univariate(data, cl))
    }
  }

  # Clean inputs
  if(missing(model)){
    model <- sym(measured_vars(data)[1])
    inform(sprintf(
      "Model not specified, using defaults set by `formula = %s`. Override this using `formula`.",
      expr_text(model)
      ))
  }
  else{
    # Capture model expression
    model <- enquo(model)
    model <- eval_tidy(expr(enexpr(!!quo_squash(model))), env = get_env(model))
  }
  
  if(is.null(model_lhs(model))){
    model <- new_formula(lhs = sym(measured_vars(data)[1]), rhs = model_rhs(model))
    inform(sprintf(
      "Response not specified, automatically selected `%s`. Override this in the `formula`.",
      expr_text(measured_vars(data)[1])
    ))
  }
  
  # Parse model
  quos(data = data,
       model = model,
       !!!parse_model_lhs(model_lhs(model), data),
       !!!parse_model_rhs(model_rhs(model), specials),
       cl = cl)
}

parse_model_rhs <- function(model_rhs, specials){
  model_rhs %>%
    parse_specials(specials = specials) %>%
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
  transformation <- as_transformation(expr, data=data)
  response <- fn_fmls(transformation)$x
  quos(
    response = response,
    transformation = transformation
  )
}