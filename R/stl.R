#' @inherit forecast::mstl
#' 
#' @param data A tsibble.
#' @param formula Decomposition specification.
#' @param ... Other arguments passed to \code{\link[forecast]{mstl}}. 
#' 
#' @examples 
#' USAccDeaths %>% STL(value ~ season(window = 10))
#' @export
STL <- function(data, formula, ...){
  # Capture user call
  cl <- call_standardise(match.call())
  
  # Coerce data
  data <- as_tsibble(data)
  
  formula <- validate_model(formula, data)
  
  # Handle multivariate inputs
  if(n_keys(data) > 1){
    return(multi_univariate(data, cl))
  }
  
  # Define specials
  specials <- new_specials_env(
    trend = function(window = NULL, degree = 1, jump = ceiling(window/10)){
      list(t.window=window, t.degree=degree, t.jump=jump)
    },
    season = function(window = 13, degree = 0, jump = ceiling(window/10), period = "all"){
      list(s.window=window, s.degree=degree, s.jump=jump, period = period)
    },
    
    .env = caller_env(),
    .required_specials = c("trend", "season")
  )
  
  # Parse model
  model_inputs <- parse_model(data, formula, specials = specials) %>% 
    flatten_first_args
  
  eval_tidy(quo(model_STL(!!!model_inputs, ...)))
}

#' @importFrom dplyr select bind_cols
#' @importFrom tibble as_tibble
#' @importFrom forecast msts mstl
model_STL <- function(data, model, response, transformation, args, period = "all", ...){
  period <- get_frequencies(args$period, data)
  args$period <- NULL
  
  # Drop unnecessary data
  data <- data %>%
    select(!!index(.), !!response)
  
  # Decompose data
  decomp <- eval_tidy(quo(mstl(msts(!!model_lhs(model), !!period), !!!args)), data = data)
  # Output tsibble decomposition
  dable(
    key_vals = as.list(data)[key_vars(data)],
    data = (data %>%
              grouped_df(key_vars(.)) %>%
              nest)$data,
    decomposition = list(
      enclass(
        data %>% select(!!index(.)) %>% bind_cols(as_tibble(decomp[,-1])),
        model = model,
        response = response,
        transformation = transformation,
        subclass = "STL"
      )
    )
  )
}

modelsplit.STL <- function(object, formula, ...){
  # y = seasadj + seas
  # y = trend + seas + rem
}

# model_STL <- function(data, model){
#   # Format args for call
#   args <- model$args %>% flatten_first_args
#   if(args$period == "all"){
#     seasonal.periods <- data %>% pull(!!index(.)) %>% guess_frequency
#     seasonal.periods <- seasonal.periods[seasonal.periods < NROW(data)/2]
#     message(paste0("Automatically selected `seasonal.periods = ", paste0(seasonal.periods, collapse=", "), "`, specify better frequencies using `seasonal.periods`."))
#   }
#   else{
#     seasonal.periods <- args$period
#   }
#   args$period <- NULL
#   seasonal.periods <- sort(seasonal.periods)
#   if(length(seasonal.periods) == 1){
#     iterate <- 1
#   }
#   
#   # Transform
#   out <- data %>%
#     transmute(
#       !!expr_deparse(f_lhs(get_expr(formula))) := !!f_lhs(formula)
#     )
#   
#   # Fill NA
#   
#   # Decompose
#   deseas <- out %>%
#     pull(!!quo_sym(f_lhs(formula)))
#   if(seasonal.periods[1] > 1){
#     seas <- rep(0, length(seasonal.periods)) %>%
#       as.list %>%
#       set_names(paste0("Seasonal", seasonal.periods))
#     
#     # out <- out %>%
#     #   mutate(
#     #     !!!set_names(list(rep(0, length(seasonal.periods))), paste0("Seasonal", seasonal.periods))
#     #   )
#     for(i in seq_len(pmax(1L, iterate))){
#       for (i in seq_along(seasonal.periods))
#       {
#         deseas <- deseas + seas[[i]]
#         fit <- stl(ts(deseas, frequency = seasonal.periods[i]), s.window = s.window[i], ...)
#         seas[[i]] <- seasonal(fit)
#         deseas <- deseas - seas[[i]]
#       }
#       # for (freq in seasonal.periods)
#       # {
#       #   deseas <- deseas + (out %>% pull(!!sym(paste0("Seasonal", freq))))
#       #   fit <- stl(ts(deseas, frequency = seasonal.periods[i]), s.window = s.window[i], ...)
#       #   out <- out %>%
#       #     mutate(!!paste0("Seasonal", freq))
#       #   seas[[i]] <- msts(seasonal(fit), seasonal.periods = msts)
#       #   attributes(seas[[i]]) <- attributes(x)
#       #   deseas <- deseas - seas[[i]]
#       # }
#       trend <- trendcycle(fit)
#     }
#     
#     out <- out %>%
#       mutate(
#         !!"Trend" := !!trend,
#         !!!seas
#       )
#   }
#   else{
#     out <- out %>%
#       mutate(
#         !!"Trend" := stats::supsmu(seq_len(NROW(data)), !!quo_sym(x))$y
#       )
#   }
#   
#   # Estimate remainder
#   out <- out %>%
#     mutate(Remainder = deseas - !!sym("Trend"))
#   
#   return(out)
# }

#' @export
components.STL <- function(object, ...){
  object
}