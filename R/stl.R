#' @inherit forecast::mstl
#' 
#' @param data A tsibble.
#' @param formula Decomposition specification.
#' @param ... Other arguments passed to \code{\link[forecast]{mstl}}. 
#' 
#' @examples 
#' USAccDeaths %>% as_tsibble %>% STL(value ~ season(window = 10))
#' elecdemand %>% STL(Demand ~ season(period = "day"))
#' @export
STL <- function(data, formula, ...){
  # Capture call
  cl <- new_quosure(match.call())
  formula <- enquo(formula)
  
  # Coerce data
  data <- as_tsibble(data)
  
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
    
    parent_env = get_env(cl)
  )
  
  # Parse Model
  model <- data %>% 
    parse_model(formula, specials = specials)
  
  eval_tidy(expr(model_STL(data, model, !!!flatten_first_args(model$args), ...)))
}

#' @importFrom dplyr select bind_cols
#' @importFrom tibble as_tibble
#' @importFrom forecast msts mstl
model_STL <- function(data, model, period = "all", ...){
  period <- get_frequencies(period, data)
  
  # Drop unnecessary data
  data <- data %>%
    select(!!index(.), !!model$response)
  
  # Decompose data
  decomp <- eval_tidy(call2("mstl", expr(msts(!!model_lhs(model$model), !!period)), !!!dots_list(...)), data = data)
  # Output tsibble decomposition
  new_tibble(list(x = list(data),
                  decomposition = list(
                    enclass(
                      data %>% select(!!index(.)) %>% bind_cols(as_tibble(decomp[,-1])),
                      !!!model,
                      subclass = "STL"
                    )
                  ),
             subclass = "dable"))
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