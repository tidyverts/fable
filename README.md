
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fable

[![Travis build
status](https://travis-ci.org/tidyverts/fable.svg?branch=master)](https://travis-ci.org/tidyverts/fable)
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Licence](https://img.shields.io/badge/licence-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
[![Coverage
status](https://codecov.io/gh/tidyverts/fable/branch/master/graph/badge.svg)](https://codecov.io/github/tidyverts/fable?branch=master)

The R package *fable* provides methods and tools for displaying and
analysing univariate time series forecasts including exponential
smoothing via state space models and automatic ARIMA modelling. Data,
model and forecast objects are all stored in a tidy format.

## Installation

You can install the **development** version from
[Github](https://github.com/tidyverts/fable)

``` r
# install.packages("devtools")
devtools::install_github("tidyverts/fable")
```

## Example

``` r
library(fable)
library(tsibbledata)
UKLungDeaths %>%
  model(ets = ETS(log(mdeaths))) %>%
  forecast
#> # A fable: 24 x 4 [1M]
#> # Key:     .model [1]
#>    .model    index mdeaths .distribution    
#>    <chr>     <mth>   <dbl> <dist>           
#>  1 ets    1980 Jan   1832. t(N(7.5, 0.0095))
#>  2 ets    1980 Feb   1854. t(N(7.5, 0.0095))
#>  3 ets    1980 Mar   1732. t(N(7.5, 0.0094))
#>  4 ets    1980 Apr   1444. t(N(7.3, 0.0089))
#>  5 ets    1980 May   1155. t(N(7.0, 0.0084))
#>  6 ets    1980 Jun   1050. t(N(7.0, 0.0082))
#>  7 ets    1980 Jul   1000. t(N(6.9, 0.0080))
#>  8 ets    1980 Aug    915. t(N(6.8, 0.0078))
#>  9 ets    1980 Sep    915. t(N(6.8, 0.0078))
#> 10 ets    1980 Oct   1081. t(N(7.0, 0.0082))
#> # ... with 14 more rows
```

You can read more about the functionality of this package and the ideas
behind it here:
<https://tidyverts.github.io/tidy-forecasting-principles/>
