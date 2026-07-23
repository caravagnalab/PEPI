# plotting functions Plot multivariate VAF distributions with cluster associated to tree nodes.

A multivariate plot is generated from a labelled dataset

## Usage

``` r
plot_multivariate(spectrum)
```

## Arguments

- spectrum:

  VAF spectrum with cluster labels

## Value

A multivariate plot

## Examples

``` r
library(dplyr)
library(ggplot2)
set.seed(1)
spectrum = data.frame(
  Nx = rbinom(50, 100, 0.3), DPx = rep(100, 50),
  Ny = rbinom(50, 100, 0.1), DPy = rep(100, 50),
  node = sample(c("-","+"), 50, replace = TRUE)
)
plot_multivariate(spectrum)
```
