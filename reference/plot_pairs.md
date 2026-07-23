# Pairs plot with divergences marked.

Useful to spot the funnel-shaped regions of parameter space that cause
divergences.

## Usage

``` r
plot_pairs(x, params = NULL)
```

## Arguments

- x:

  PEPI_Multirates object fit with `method = "sample"`.

- params:

  Vector of Stan parameter base names. Defaults to a handful of global
  scalar parameters (pairs plots grow quadratically with the number of
  parameters).

## Value

A ggplot/grid object (see
[`bayesplot::mcmc_pairs()`](https://mc-stan.org/bayesplot/reference/MCMC-scatterplots.html)).

## Examples

``` r
if (FALSE) { # \dontrun{
plot_pairs(x, params = c("lambda_n","s_epi"))
} # }
```
