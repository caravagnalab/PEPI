# Fit a Gamma distribution to two quantiles.

Finds the shape and rate of the Gamma distribution whose \`p1\` and
\`p2\` quantiles are \`q1\` and \`q2\`, by least squares on the CDF.
Useful to turn a plausible range of a rate (e.g. a switch rate) into
Gamma prior hyperparameters.

## Usage

``` r
fit_gamma_quantiles(q1, q2, p1 = 0.05, p2 = 0.95)
```

## Arguments

- q1, q2:

  Lower and upper quantile values.

- p1, p2:

  Probabilities of the two quantiles.

## Value

A list with elements \`shape\` and \`rate\`.

## Examples

``` r
g = PEPI:::fit_gamma_quantiles(0.01, 0.5)
stats::qgamma(c(0.05, 0.95), shape = g$shape, rate = g$rate)
#> [1] 0.009999387 0.499998982
```
