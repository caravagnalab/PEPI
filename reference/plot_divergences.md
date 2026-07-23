# Divergent transitions.

Divergences concentrating in a region of parameter space indicate the
sampler struggled there and posterior estimates may be biased.

## Usage

``` r
plot_divergences(x)
```

## Arguments

- x:

  PEPI_Multirates object fit with `method = "sample"`.

## Value

A ggplot object (see
[`bayesplot::mcmc_nuts_divergence()`](https://mc-stan.org/bayesplot/reference/MCMC-nuts.html)).

## Examples

``` r
if (FALSE) { # \dontrun{
plot_divergences(x)
} # }
```
