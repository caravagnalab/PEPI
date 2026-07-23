# Autocorrelation of the sampled chains.

Autocorrelation of the sampled chains.

## Usage

``` r
plot_acf(x, params = NULL, lags = 20)
```

## Arguments

- x:

  PEPI_Multirates object fit with `method = "sample"`.

- params:

  Vector of Stan parameter base names. Defaults to a handful of global
  scalar parameters.

- lags:

  Number of lags to show.

## Value

A ggplot object (see
[`bayesplot::mcmc_acf()`](https://mc-stan.org/bayesplot/reference/MCMC-diagnostics.html)).

## Examples

``` r
if (FALSE) { # \dontrun{
plot_acf(x)
} # }
```
