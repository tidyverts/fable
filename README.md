
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
  ETS(log(mdeaths)) %>%
  forecast
#> # A tsibble: 24 x 3 [1M]
#>       index  mean distribution         
#>       <mth> <dbl> <dist>               
#>  1 1980 Jan 1832. t(N(7.5, sd = 0.098))
#>  2 1980 Feb 1854. t(N(7.5, sd = 0.098))
#>  3 1980 Mar 1732. t(N(7.5, sd = 0.097))
#>  4 1980 Apr 1444. t(N(7.3, sd = 0.094))
#>  5 1980 May 1155. t(N(7, sd = 0.092))  
#>  6 1980 Jun 1050. t(N(7, sd = 0.09))   
#>  7 1980 Jul 1000. t(N(6.9, sd = 0.09)) 
#>  8 1980 Aug  915. t(N(6.8, sd = 0.089))
#>  9 1980 Sep  915. t(N(6.8, sd = 0.089))
#> 10 1980 Oct 1081. t(N(7, sd = 0.091))  
#> # ... with 14 more rows
```

You can read more about the functionality of this package and the ideas
behind it here:
<https://tidyverts.github.io/tidy-forecasting-principles/>
