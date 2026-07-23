# Build the Stan data list for the unified multirates model.

Turns the tables stored on a \`PEPI_Multirates\` object (see
[`init_multirates()`](https://caravagnalab.github.io/PEPI/reference/init_multirates.md))
into the exact data list required by `inst/multirates_positive_s.stan`,
plus an \`index_maps\` list that lets getters translate Stan array
indices back into ids/times.

## Usage

``` r
build_stan_data_multirates(
  x,
  mu,
  l,
  t_min,
  ms_epi,
  sigma_epi,
  alpha_lambda,
  beta_lambda,
  alpha_n_wt,
  beta_n_wt,
  alpha_p_wt,
  beta_p_wt,
  min_kappa,
  max_kappa,
  min_sigma_count,
  max_sigma_count,
  include_poisson = 1L,
  include_trunk = 1L
)
```

## Arguments

- x:

  PEPI_Multirates object.

- mu:

  Mutation rate per division per bp per allele.

- l:

  Length of the genome.

- t_min:

  Lower bound for the tmrca prior.

- ms_epi, sigma_epi:

  Lognormal prior hyperparameters for s_epi.

- alpha_lambda, beta_lambda:

  Gamma prior hyperparameters for lambda_n.

- alpha_n_wt, beta_n_wt, alpha_p_wt, beta_p_wt:

  Gamma prior hyperparameters for the wt switch rates
  omega_n_wt/omega_p_wt.

- min_kappa, max_kappa:

  Uniform prior bounds for kappa.

- min_sigma_count, max_sigma_count:

  Uniform prior bounds for sigma_count.

- include_poisson, include_trunk:

  Stan model switches (0/1).

## Value

A list with elements \`data\` (the Stan data list) and \`index_maps\`.

## Examples

``` r
if (FALSE) { # \dontrun{
build_stan_data_multirates(x, mu = 1e-7, l = 2.7e9, t_min = 0,
  ms_epi = 0, sigma_epi = 0.5, alpha_lambda = 1, beta_lambda = 1,
  alpha_n_wt = 1, beta_n_wt = 10, alpha_p_wt = 1, beta_p_wt = 10,
  min_kappa = 10, max_kappa = 1000, min_sigma_count = 0.01, max_sigma_count = 1)
} # }
```
