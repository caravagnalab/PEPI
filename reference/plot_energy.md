# Hamiltonian Monte Carlo energy diagnostic.

Overlapping marginal and transition energy distributions indicate
efficient exploration; a large mismatch signals poor mixing.

## Usage

``` r
plot_energy(x)
```

## Arguments

- x:

  PEPI_Multirates object fit with `method = "sample"`.

## Value

A ggplot object (see
[`bayesplot::mcmc_nuts_energy()`](https://mc-stan.org/bayesplot/reference/MCMC-nuts.html)).

## Examples

``` r
if (FALSE) { # \dontrun{
plot_energy(x)
} # }
```
