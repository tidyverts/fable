globalVariables(c("self", "origin"))

#' @importFrom stats model.frame
model_xreg <- function(...) {
  model_formula <- new_formula(
    lhs = NULL,
    rhs = reduce(enexprs(...), function(.x, .y) call2("+", .x, .y))
  )
  env <- map(enquos(...), get_env)
  env[map_lgl(env, compose(is_empty, env_parents))] <- NULL
  env <- if (!is_empty(env)) get_env(env[[1]]) else base_env()
  out <- model.frame(model_formula, data = env, na.action = stats::na.pass)
}

no_xreg <- function(...) {
  abort("Exogenous regressors are not supported for this model type.")
}

trend <- function(x, knots = NULL, origin = NULL) {
  UseMethod("trend")
}

trend.tbl_ts <- function(x, knots = NULL, origin = NULL) {
  idx_num <- x[[index_var(x)]] %>% as.double()
  knots_num <- if (is.null(knots)) {
    NULL
  } else {
    knots %>% as.double()
  }
  index_interval <- interval(x) %>% default_time_units()
  idx_num <- idx_num / index_interval
  knots_num <- knots_num / index_interval
  if (!is.null(origin)) {
    origin <- as.double(origin) / index_interval
  }

  trend(idx_num, knots_num, origin)
}

trend.numeric <- function(x, knots = NULL, origin = NULL) {
  if (!is.null(origin)) {
    origin <- origin - 1 # trend should count from 1
    x <- x - origin
    knots <- knots - origin
  }
  knots_exprs <- knots %>%
    map(function(.x) pmax(0, x - .x)) %>%
    set_names(map_chr(knots, function(.x) paste0("trend_", format(.x))))
  tibble(
    trend = x,
    !!!knots_exprs
  )
}

season <- function(x, period) {
  UseMethod("season")
}

season.tbl_ts <- function(x, period) {
  idx_num <- x[[index_var(x)]] %>% as.double()
  index_interval <- interval(x) %>% default_time_units()
  idx_num <- idx_num / index_interval
  period <- get_frequencies(period, x, .auto = "smallest")

  season(idx_num, period)
}

season.numeric <- function(x, period) {
  season_exprs <- period %>%
    map(function(.x) expr(factor(floor((x %% (!!.x)) + 1), levels = seq_len(!!.x)))) %>%
    set_names(names(period) %||% paste0("season_", period))
  tibble(!!!season_exprs)
}

fourier <- function(x, period, K, origin = NULL) {
  UseMethod("fourier")
}

fourier.tbl_ts <- function(x, period, K, origin = NULL) {
  idx_num <- x[[index_var(x)]] %>% as.double()
  index_interval <- interval(x) %>% default_time_units()
  idx_num <- idx_num / index_interval
  if (!is.null(origin)) {
    origin <- as.double(origin) / index_interval
  }
  period <- get_frequencies(period, x, .auto = "smallest")

  fourier(idx_num, period, K, origin)
}

fourier.numeric <- function(x, period, K, origin = NULL) {
  if (length(period) != length(K)) {
    abort("Number of periods does not match number of orders")
  }
  if (any(2 * K > period)) {
    abort("K must be not be greater than period/2")
  }

  fourier_exprs <- map2(
    as.numeric(period), K,
    function(period, K) {
      set_names(seq_len(K) / period, paste0(seq_len(K), "_", round(period)))
    }
  ) %>%
    invoke(c, .) %>%
    .[!duplicated(.)] %>%
    map2(., names(.), function(p, name) {
      out <- exprs(C = cospi(2 * !!p * x))
      if (abs(2 * p - round(2 * p)) > .Machine$double.eps) {
        out <- c(out, exprs(S = sinpi(2 * !!p * x)))
      }
      names(out) <- paste0(names(out), name)
      out
    }) %>%
    set_names(NULL) %>%
    unlist(recursive = FALSE)

  tibble(!!!fourier_exprs)
}

#' Common exogenous regressors
#'
#' These special functions provide interfaces to more complicated functions within
#' the model formulae interface.
#'
#' @section Specials:
#'
#' \subsection{trend}{
#' The `trend` special includes common linear trend regressors in the model. It also supports piecewise linear trend via the `knots` argument.
#' \preformatted{
#' trend(knots = NULL, origin = NULL)
#' }
#'
#' \tabular{ll}{
#'   `knots`    \tab A vector of times (same class as the data's time index) identifying the position of knots for a piecewise linear trend.\cr
#'   `origin`   \tab An optional time value to act as the starting time for the trend.
#' }
#' }
#'
#' \subsection{season}{
#' The `season` special includes seasonal dummy variables in the model.
#' \preformatted{
#' season(period = NULL)
#' }
#'
#' \tabular{ll}{
#'   `period`   \tab The periodic nature of the seasonality. This can be either a number indicating the number of observations in each seasonal period, or text to indicate the duration of the seasonal window (for example, annual seasonality would be "1 year").
#' }
#' }
#'
#' \subsection{fourier}{
#' The `fourier` special includes seasonal fourier terms in the model. The maximum order of the fourier terms must be specified using `K`.
#' \preformatted{
#' fourier(period = NULL, K, origin = NULL)
#' }
#'
#' \tabular{ll}{
#'   `period`   \tab The periodic nature of the seasonality. This can be either a number indicating the number of observations in each seasonal period, or text to indicate the duration of the seasonal window (for example, annual seasonality would be "1 year"). \cr
#'   `K`        \tab The maximum order of the fourier terms.\cr
#'   `origin`   \tab An optional time value to act as the starting time for the fourier series.
#' }
#' }
#'
#' @format NULL
#'
#' @rdname common_xregs
common_xregs <- list(
  trend = function(knots = NULL, origin = NULL) {
    if (is.null(origin)) {
      if (is.null(self$origin)) {
        self$origin <- self$data[[index_var(self$data)]][[1]]
      }
      origin <- self$origin
    }
    fable:::trend(self$data, knots, origin) %>% as.matrix()
  },
  season = function(period = NULL) {
    fable:::season(self$data, period) %>% as_model_matrix()
  },
  fourier = function(period = NULL, K, origin = NULL) {
    if (is.null(origin)) {
      if (is.null(self$origin)) {
        self$origin <- self$data[[index_var(self$data)]][[1]]
      }
      origin <- self$origin
    }
    fable:::fourier(self$data, period, K, origin) %>% as.matrix()
  }
)

as_model_matrix <- function(tbl) {
  stats::model.matrix(~., data = tbl)[, -1, drop = FALSE]
}
