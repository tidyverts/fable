
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fable <a href='https://fable.tidyverts.org'><img src='man/figures/logo.png' align="right" height="138.5" /></a>

[![Travis build
status](https://travis-ci.org/tidyverts/fable.svg?branch=master)](https://travis-ci.org/tidyverts/fable)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/tidyverts/fable?branch=master&svg=true)](https://ci.appveyor.com/project/tidyverts/fable)
[![Coverage
status](https://codecov.io/gh/tidyverts/fable/branch/master/graph/badge.svg)](https://codecov.io/github/tidyverts/fable?branch=master)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/fable)](https://cran.r-project.org/package=fable)
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)

The R package *fable* provides methods and tools for displaying and
analysing univariate time series forecasts including exponential
smoothing via state space models and automatic ARIMA modelling. Data,
model and forecast objects are all stored in a tidy format.

## Installation

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
  forecast %>% 
  autoplot(filter(aus_retail, year(Month) > 2010), level = NULL)
```

<img src="man/figures/README-example-1.png" width="100%" />

You can read more about the functionality of this package and the ideas
behind it here:
<https://tidyverts.github.io/tidy-forecasting-principles/>
