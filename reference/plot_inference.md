# Plot posterior and prior distributions of the unified multirates model's parameters.

Posterior and prior draws histograms are plotted for any required
parameter, using the raw ("param" family) rows of
[`get_posterior_multirates()`](https://caravagnalab.github.io/PEPI/reference/get_posterior_multirates.md)'s
output.

## Usage

``` r
plot_inference(x, params = NULL, groups = NULL)
```

## Arguments

- x:

  PEPI_Multirates object, after
  [`get_posterior_multirates()`](https://caravagnalab.github.io/PEPI/reference/get_posterior_multirates.md).

- params:

  A vector of canonical parameter names (e.g. "lambda_n", "s_driver").

- groups:

  A vector of groups to include
  (wt/driver/driver_n/driver_p/dc/clade_wt/global).

## Value

A plot with posterior and prior distributions.

## Examples

``` r
if (FALSE) { # \dontrun{
plot_inference(x,params = c("lambda_n","s_epi"))
} # }
```
