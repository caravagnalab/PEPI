# Plot predicted vs observed mutation counts for the unified multirates model.

Posterior-predictive distributions (violin) vs observed counts (point),
using
[`get_posterior_multirates()`](https://caravagnalab.github.io/PEPI/reference/get_posterior_multirates.md)'s
"m" family.

## Usage

``` r
plot_mutations(x)
```

## Arguments

- x:

  PEPI_Multirates object, after
  [`get_posterior_multirates()`](https://caravagnalab.github.io/PEPI/reference/get_posterior_multirates.md).

## Value

Plot of predicted vs observed mutation counts.

## Examples

``` r
if (FALSE) { # \dontrun{
plot_mutations(x)
} # }
```
