
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fable <a href='https://fable.tidyverts.org'><img src='man/figures/logo.png' align="right" height="138.5" /></a>

![R build
status](https://github.com/tidyverts/fable/workflows/R-CMD-check/badge.svg)
[![Coverage
status](https://codecov.io/gh/tidyverts/fable/branch/master/graph/badge.svg)](https://codecov.io/github/tidyverts/fable?branch=master)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/fable)](https://cran.r-project.org/package=fable)
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)

The R package *fable* provides a collection of commonly used univariate
and multivariate time series forecasting models including exponential
smoothing via state space models and automatic ARIMA modelling. These
models work within the fable framework, which provides the tools to
evaluate, visualise, and combine models in a workflow consistent with
the tidyverse.

## Installation

The can install the **stable** version from
[CRAN](https://cran.r-project.org/package=fable):

``` r
install.packages("fable")
```

You can install the **development** version from
[GitHub](https://github.com/tidyverts/fable)

``` r
# install.packages("remotes")
remotes::install_github("tidyverts/fable")
```

Installing this software requires a compiler

## Example

``` r
library(fable)
library(tsibble)
library(tsibbledata)
library(lubridate)
#> Warning: package 'lubridate' was built under R version 3.6.3
library(dplyr)
aus_retail %>%
  filter(
    State %in% c("New South Wales", "Victoria"),
    Industry == "Department stores"
  ) %>% 
  model(
    ets = ETS(box_cox(Turnover, 0.3)),
    arima = ARIMA(log(Turnover)),
    snaive = SNAIVE(Turnover)
  ) %>%
  forecast(h = "2 years") %>% 
  autoplot(filter(aus_retail, year(Month) > 2010), level = NULL)
```

<img src="man/figures/README-example-1.png" width="100%" />

## Learning to forecast with fable

  - The pkgdown site describes all models provided by fable, and how
    they are used: <http://fable.tidyverts.org/>
  - The forecasting principles and practices online textbook provides an
    introduction to time series forecasting using fable:
    <https://otexts.com/fpp3/> (WIP)
    <!-- - A quick start functionality guide can be found here: https://tidyverts.github.io/tidy-forecasting-principles/ (WIP) -->

## Getting help

  - Questions about forecasting can be asked on [Cross
    Validated](http://stats.stackexchange.com/tags/forecasting).

  - Common questions about the fable package are often found on [Stack
    Overflow](http://stackoverflow.com/tags/fable+r). You can use this
    to ask for help if the question isn’t already answered. A [minimally
    reproducible example](https://www.tidyverse.org/help/) that
    describes your issue is the best way to ask for help\!
