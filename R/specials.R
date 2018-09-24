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
    rhs = reduce(enexprs(...), function(.x, .y) call2("+", .x, .y))
  )
  model.frame(model_formula, data = .data)
}

no_xreg <- function(...){
  abort("Exogenous regressors are not supported for this model type.")
}

origin <- NULL

trend <- function(x, knots = NULL, origin = NULL){
  UseMethod("trend")
}

trend.tbl_ts <- function(x, knots = NULL, origin = NULL){
  idx_num <- x[[expr_text(tsibble::index(x))]] %>% units_since
  knots_num <- if(is.null(knots)){NULL} else {knots %>% units_since}
  index_interval <- idx_num %>% time_unit()
  idx_num <- idx_num/index_interval
  knots_num <- knots_num/index_interval
  if(!is.null(origin)){
    origin <- units_since(origin)/index_interval
  }
  
  trend(idx_num, knots_num, origin)
}

trend.numeric <- function(x, knots = NULL, origin = NULL){
  if(!is.null(origin)){
    x <- x - origin
    knots <- knots - origin
  }
  knots_exprs <- knots %>%
    map(function(.x) pmax(0, x-.x)) %>%
    set_names(map_chr(knots, function(.x) paste0("trend_",format(.x))))
  tibble(trend = x,
         !!!knots_exprs)
}

season <- function(x, period){
  UseMethod("season")
}

season.tbl_ts <- function(x, period){
  idx_num <- x[[expr_text(tsibble::index(x))]] %>% units_since
  index_interval <- idx_num %>% time_unit()
  idx_num <- idx_num/index_interval
  period <- get_frequencies(period, x)
  
  season(idx_num, period)
}

season.numeric <- function(x, period){
  season_exprs <- period %>% 
    map(function(.x) expr(as.factor((x%%(!!.x))+1))) %>%
    set_names(names(period)%||%paste0("season_", period))
  tibble(!!!season_exprs)
}

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
    map2(., names(.), function(p, name){
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