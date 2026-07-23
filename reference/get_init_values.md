# Provide a initialization parameters for a PEPI VAF fit

Return a list of initialization values for some parameters through some
euristics.

## Usage

``` r
get_init_values(spectrum, K = 10, alpha = 10, samples = 1, pi_cutoff = 0.01)
```

## Arguments

- spectrum:

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

a list with initialization values

## Examples

``` r
if (FALSE) { # \dontrun{
get_init_values(spectrum,K = 10,alpha = 10,samples = 1, pi_cutoff = 0.01)
} # }
```
