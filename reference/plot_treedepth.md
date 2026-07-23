# NUTS treedepth diagnostic.

Chains repeatedly hitting `max_treedepth` indicate the sampler is being
cut off before making a U-turn, which hurts efficiency (though not
necessarily validity).

## Usage

``` r
plot_treedepth(x)
```

## Arguments

- x:

  PEPI_Multirates object fit with `method = "sample"`.

## Value

A ggplot object (see
[`bayesplot::mcmc_nuts_treedepth()`](https://mc-stan.org/bayesplot/reference/MCMC-nuts.html)).

## Examples

``` r
if (FALSE) { # \dontrun{
plot_treedepth(x)
} # }
```
