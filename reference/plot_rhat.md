# Rhat (potential scale reduction factor) overview.

One point per parameter; values close to 1 indicate the chains have
converged to a common distribution. Values above 1.05 are a red flag.

## Usage

``` r
plot_rhat(x, params = NULL)
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
[`bayesplot::mcmc_rhat()`](https://mc-stan.org/bayesplot/reference/MCMC-diagnostics.html)).

## Examples

``` r
if (FALSE) { # \dontrun{
plot_rhat(x)
} # }
```
