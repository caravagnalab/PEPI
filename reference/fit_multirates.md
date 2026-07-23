# Fit the unified multirates model.

Compiles and runs `inst/multirates_positive_s.stan` on a
`PEPI_Multirates` object built with
[`init_multirates()`](https://caravagnalab.github.io/PEPI/reference/init_multirates.md).

## Usage

``` r
fit_multirates(
  x,
  cmdstan_path = cmdstanr::cmdstan_path(),
  method = c("variational", "sample"),
  ndraws = 1000,
  chains = 4,
  seed = 45,
  init = NULL,
  mu = 1e-07,
  l = 2.7e+09,
  t_min = 0,
  ms_epi = 0,
  sigma_epi = 0.5,
  alpha_lambda = 1,
  beta_lambda = 1,
  alpha_n_wt = 1,
  beta_n_wt = 10,
  alpha_p_wt = 1,
  beta_p_wt = 10,
  min_kappa = 10,
  max_kappa = 1000,
  min_sigma_count = 0.01,
  max_sigma_count = 1,
  include_poisson = 1L,
  include_trunk = 1L
)
```

## Arguments

- x:

  PEPI_Multirates object.

- cmdstan_path:

  String specifying the path to the cmdstan installation.

- method:

  "variational" (default, fast approximate posterior) or "sample" (full
  MCMC).

- ndraws:

  Number of posterior draws (variational: output_samples; sample:
  iter_sampling).

- chains:

  Number of MCMC chains (only used when method = "sample").

- seed:

  Seed of the computation.

- init:

  List of initialization values, or NULL to use a built-in default (see
  [`.multirates_default_init()`](https://caravagnalab.github.io/PEPI/reference/dot-multirates_default_init.md)).

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

PEPI_Multirates object with \`inference\$multirates\`,
\`stan_data\$multirates\` and \`index_maps\` populated.

## Examples

``` r
if (FALSE) { # \dontrun{
fit_multirates(x, cmdstan_path = cmdstanr::cmdstan_path(),
  method = "variational", ndraws = 1000, seed = 45,
  mu = 1e-7, l = 2.7e9, t_min = 0, ms_epi = 0, sigma_epi = 0.5,
  alpha_lambda = 1, beta_lambda = 1, alpha_n_wt = 1, beta_n_wt = 10,
  alpha_p_wt = 1, beta_p_wt = 10, min_kappa = 10, max_kappa = 1000,
  min_sigma_count = 0.01, max_sigma_count = 1)
} # }
```
