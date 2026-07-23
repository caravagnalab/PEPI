# One-call convergence dashboard.

Arranges Rhat, effective sample size, divergences and the energy
diagnostic into a single 2x2 figure - a quick first look at whether an
MCMC fit is trustworthy.

## Usage

``` r
plot_diagnostics(x, params = NULL)
```

## Arguments

- x:

  PEPI_Multirates object fit with `method = "sample"`.

- params:

  Vector of Stan parameter base names for the Rhat/ESS panels. Defaults
  to every sampled parameter.

## Value

A combined plot (see
[`ggpubr::ggarrange()`](https://rpkgs.datanovia.com/ggpubr/reference/ggarrange.html)).

## Examples

``` r
if (FALSE) { # \dontrun{
plot_diagnostics(x)
} # }
```
