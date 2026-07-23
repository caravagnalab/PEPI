# Clusters mutation with VIBER, a package that implements a variational non parametric Bayesian model to fit multi-variate Binomial mixtures

Clusters with centroids VAF coordiates and relative proportion are
given.

## Usage

``` r
get_viber_clusters(data, K = 10, alpha = 10, samples = 1, pi_cutoff = 0.01)
```

## Arguments

- data:

  Mutation data with number of variants and depth for any mutation and
  sample

- K:

  Maximum number of clusters

- alpha:

  Dirichlet concentration parameter

- samples:

  Number of fits computed by the algorithm

- pi_cutoff:

  Cutoff on mixing proportions to filter clusters and reassign mutations

## Value

a tibble with clusters names, mixing proportion and vaf coordinates

## Examples

``` r
if (FALSE) { # \dontrun{
# requires the VIBER package (github.com/caravagn/VIBER), not on CRAN
set.seed(1)
data = data.frame(
  Nx = rbinom(50, 100, 0.3), DPx = rep(100, 50),
  Ny = rbinom(50, 100, 0.1), DPy = rep(100, 50)
)
get_viber_clusters(data, K = 10, alpha = 1, samples = 1, pi_cutoff = 0.01)
} # }
```
