tbl_xreg <- function(...){
  tibble(...)
}

exprs_xreg <- function(...){
  exprs(...)
}

#' @importFrom stats model.frame
model_xreg <- function(...){
  model_formula <- new_formula(
    lhs = NULL,
    rhs = reduce(enexprs(...), ~ call2("+", .x, .y))
  )
  model.frame(model_formula, data = .data)
}

no_xreg <- function(...){
  abort("Exogenous regressors are not supported for this model type.")
}

origin <- NULL

trend <- function(data, knots = NULL, origin = NULL){
  idx_num <- data[[expr_text(tsibble::index(data))]] %>% units_since
  knots_num <- if(is.null(knots)){NULL} else {knots %>% units_since}
  if(!is.null(origin)){
    origin_num <- units_since(origin)
    idx_num <- idx_num - origin_num
    knots_num <- knots_num - origin_num
  }
  index_interval <- idx_num %>% time_unit()
  idx_num <- idx_num/index_interval
  knots_num <- knots_num/index_interval
  knots_exprs <- knots_num %>%
    map(~ pmax(0, idx_num-.x)) %>%
    set_names(map_chr(knots, ~ paste0("trend_",format(.x))))
  tibble(trend = idx_num,
         !!!knots_exprs)
}

season <- function(data, period){
  idx_num <- data[[expr_text(tsibble::index(data))]] %>% units_since
  index_interval <- idx_num %>% time_unit()
  idx_num <- idx_num/index_interval
  period <- get_frequencies(period, data)
  season_exprs <- period %>% 
    map(~ expr(as.factor((idx_num%%(!!period))+1))) %>%
    set_names(names(period))
  tibble(!!!season_exprs)
}

#' @importFrom purrr imap
fourier <- function(data, period, K, origin = NULL){ 
  idx_num <- data[[expr_text(tsibble::index(data))]] %>% units_since
  if(!is.null(origin)){
    idx_num <- idx_num - units_since(origin)
  }
  index_interval <- idx_num %>% time_unit()
  idx_num <- idx_num/index_interval
  period <- get_frequencies(period, data)
  if (length(period) != length(K)) {
    abort("Number of periods does not match number of orders")
  }
  if (any(2 * K > period)) {
    abort("K must be not be greater than period/2")
  }
  
  fourier_exprs <- map2(as.numeric(period), K,
   function(period, K){
     set_names(seq_len(K)/period, paste0(seq_len(K), "_", round(period)))
   }) %>%
    invoke(c, .) %>%
    .[!duplicated(.)] %>%
    imap(function(p, name){
      out <- exprs(C = cospi(2 * !!p * idx_num))
      if(abs(2 * p - round(2 * p)) > .Machine$double.eps){
        out <- c(out, exprs(S = sinpi(2 * !!p * idx_num)))
      }
      names(out) <- paste0(names(out), name)
      out
    }) %>%
    set_names(NULL) %>%
    unlist(recursive = FALSE)
  
  tibble(!!!fourier_exprs)
}