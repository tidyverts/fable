#' @inherit forecast::mstl
#' @importFrom forecast seasonal trendcycle
#' @importFrom stats ts stl
#' @importFrom tsibble guess_frequency index
#' @importFrom dplyr pull select mutate transmute
#' 
#' @param data A tsibble.
#' @param x The column from \code{data} to decompose.
#' @param seasonal.periods A vector of seasonal frequencies used for the seasonal decomposition.
#' 
#' @examples 
#' \dontrun{
#' USAccDeaths %>% as_tsibble %>% STL
#' }
#' @export
STL <- function(data, x = NULL, seasonal.periods = NULL, iterate = 2, s.window = 13, ...){
  x <- enquo(x)
  if(is.null(eval_tidy(x, data=data))){
    potential_x <- colnames(data %>% select(-!!index(.), drop=TRUE))
    message("Decomposing data from column `x = ", potential_x[1], "`. Override this using `x`.")
    x <- set_expr(x, sym(potential_x[1]))
  }
  if(is.null(seasonal.periods)){
    seasonal.periods <- data %>% pull(!!index(.)) %>% guess_frequency
    seasonal.periods <- seasonal.periods[seasonal.periods < NROW(data)/2]
    message(paste0("Automatically selected `seasonal.periods = ", paste0(seasonal.periods, collapse=", "), "`, specify better frequencies using `seasonal.periods`."))
  }
  seasonal.periods <- sort(seasonal.periods)
  if(length(seasonal.periods) == 1){
    iterate <- 1
  }

  # Transform
  out <- data %>%
    transmute(
      !!expr_deparse(get_expr(index(data))) := !!index(data),
      !!expr_deparse(get_expr(x)) := !!x
    )
  
  # Fill NA
  
  # Decompose
  deseas <- out %>%
    pull(!!quo_sym(x))
  if(seasonal.periods[1] > 1){
    seas <- rep(0, length(seasonal.periods)) %>%
      as.list %>%
      set_names(paste0("Seasonal", seasonal.periods))
      
    # out <- out %>%
    #   mutate(
    #     !!!set_names(list(rep(0, length(seasonal.periods))), paste0("Seasonal", seasonal.periods))
    #   )
    for(i in seq_len(pmax(1L, iterate))){
      for (i in seq_along(seasonal.periods))
      {
        deseas <- deseas + seas[[i]]
        fit <- stl(ts(deseas, frequency = seasonal.periods[i]), s.window = s.window[i], ...)
        seas[[i]] <- seasonal(fit)
        deseas <- deseas - seas[[i]]
      }
      # for (freq in seasonal.periods)
      # {
      #   deseas <- deseas + (out %>% pull(!!sym(paste0("Seasonal", freq))))
      #   fit <- stl(ts(deseas, frequency = seasonal.periods[i]), s.window = s.window[i], ...)
      #   out <- out %>%
      #     mutate(!!paste0("Seasonal", freq))
      #   seas[[i]] <- msts(seasonal(fit), seasonal.periods = msts)
      #   attributes(seas[[i]]) <- attributes(x)
      #   deseas <- deseas - seas[[i]]
      # }
      trend <- trendcycle(fit)
    }
    
    out <- out %>%
      mutate(
        !!"Trend" := !!trend,
        !!!seas
      )
  }
  else{
    out <- out %>%
      mutate(
        !!"Trend" := stats::supsmu(seq_len(NROW(data)), !!quo_sym(x))$y
      )
  }
  
  # Estimate remainder
  out <- out %>%
    mutate(Remainder = deseas - !!sym("Trend"))
  
  return(out)
}