# Effective sample size (ratio to total draws) overview.

One point per parameter; low ratios indicate highly autocorrelated
chains that need more draws for the same effective precision.

## Usage

``` r
plot_ess(x, params = NULL)
```

## Arguments

- x:

  PEPI_Multirates object fit with `method = "sample"`.

- params:

  Vector of Stan parameter base names to include. Defaults to every
  sampled parameter (excludes internal latent-state arrays and generated
  quantities).

## Value

A ggplot object (see
[`bayesplot::mcmc_neff()`](https://mc-stan.org/bayesplot/reference/MCMC-diagnostics.html)).

## Examples

``` r
if (FALSE) { # \dontrun{
plot_ess(x)
} # }
```
