# Trace plots of the sampled chains.

Trace plots of the sampled chains.

## Usage

``` r
plot_trace(x, params = NULL, divergences = TRUE)
```

## Arguments

- x:

  PEPI_Multirates object fit with `method = "sample"`.

- params:

  Vector of Stan parameter base names (no brackets - all indexed
  elements of a vector/array parameter are included automatically).
  Defaults to a handful of global scalar parameters.

- divergences:

  If TRUE (default), divergent transitions are marked on the trace.

## Value

A ggplot object (see
[`bayesplot::mcmc_trace()`](https://mc-stan.org/bayesplot/reference/MCMC-traces.html)).

## Examples

``` r
if (FALSE) { # \dontrun{
plot_trace(x, params = c("lambda_n","s_epi"))
} # }
```
