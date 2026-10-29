# Fit the unified multirates model.

Compiles and runs `inst/multirates_positive_s.stan` on a
`PEPI_Multirates` object built with
[`init_multirates()`](https://caravagnalab.github.io/PEPI/reference/init_multirates.md).

## Usage

``` r
fit_multirates(
  x,
  cmdstan_path = cmdstanr::cmdstan_path(),
  seed = 45,
  ndraws = 1000,
  chains = 4,
  parallel_chains = 4,
  adapt_delta = 0.99,
  max_treedepth = 15,
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
  include_poisson_driver = include_poisson,
  include_poisson_driver_n = include_poisson,
  include_poisson_driver_p = include_poisson,
  include_poisson_dc = 1L,
  include_trunk = 1L
)
```

## Arguments

- x:

  PEPI_Multirates object.

- cmdstan_path:

  String specifying the path to the cmdstan installation.

- seed:

  Seed of the computation.

- ndraws:

  Number of posterior draws.

- chains:

  Number of MCMC chains.

- parallel_chains:

  Number of parallel MCMC chains.

- adapt_delta:

  Target average acceptance probability.

- max_treedepth:

  Maximum number of trajectory steps the NUTS sampler can take per
  iteration.

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

- include_poisson:

  Default (0/1) of the Poisson terms on the driver mutations, used by
  the per-type switches below when they are not given.

- include_poisson_driver, include_poisson_driver_n,
  include_poisson_driver_p:

  Switches (0/1) of the Poisson term on the mutations of the driver,
  driver_n and driver_p populations (for driver_p, an upper bound).

- include_poisson_dc:

  Switch (0/1) of the Poisson term on the driver mutations of the dc
  populations; the clade terms are always included.

- include_trunk:

  Switch (0/1) of the trunk regularization of tmrca.

## Value

PEPI_Multirates object with \`inference\$multirates\`,
\`stan_data\$multirates\` and \`index_maps\` populated.

## Examples

``` r
# \donttest{
source(system.file("examples", "tables_example.R", package = "PEPI"))
x = fit_multirates(x, seed = 31, chains = 2, parallel_chains = 2, ndraws = 200,
  init = rep(list(list(omega_p_wt = 3e-3, omega_p_driver = 3e-3)), 2),
  mu = mu$mean, l = 51e06, t_min = -20, ms_epi = -0.8, sigma_epi = 0.9,
  alpha_lambda = 10, beta_lambda = 10, alpha_n_wt = 1, beta_n_wt = 1e4,
  alpha_p_wt = 0.2, beta_p_wt = 10, min_kappa = 3000, max_kappa = 5000,
  min_sigma_count = 0.05, max_sigma_count = 0.15, include_poisson = 1, include_trunk = 0)
#> CmdStan path set to: /home/runner/.cmdstan/cmdstan-2.35.0
#> Init values were only set for a subset of parameters. 
#> Missing init values for the following parameters:
#>  - chain 1: lambda_n, s_epi, omega_n_wt, bc_wt, U_wt, W_wt, tmrca, t_clade_wt, bc_clade_wt, s_driver, t_driver, bc_driver, U_driver, W_driver, s_driver_n, t_driver_n, bc_driver_n, s_driver_p, t_driver_p, bc_driver_p, s_dc, delta_t_dc, bc_driver_dc, bc_clade_dc, U_dc, W_dc, omega_p_dc, omega_n_dc, omega_n_driver, kappa, sigma_count
#>  - chain 2: lambda_n, s_epi, omega_n_wt, bc_wt, U_wt, W_wt, tmrca, t_clade_wt, bc_clade_wt, s_driver, t_driver, bc_driver, U_driver, W_driver, s_driver_n, t_driver_n, bc_driver_n, s_driver_p, t_driver_p, bc_driver_p, s_dc, delta_t_dc, bc_driver_dc, bc_clade_dc, U_dc, W_dc, omega_p_dc, omega_n_dc, omega_n_driver, kappa, sigma_count
#> 
#> To disable this message use options(cmdstanr_warn_inits = FALSE).
#> Running MCMC with 2 parallel chains...
#> 
#> Chain 1 Iteration:    1 / 1200 [  0%]  (Warmup) 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: gamma_lpdf: Random variable is inf, but must be positive finite! (in '/tmp/Rtmp9AUyja/model-1be44a13bd0d.stan', line 647, column 0 to column 44)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: gamma_lpdf: Random variable is inf, but must be positive finite! (in '/tmp/Rtmp9AUyja/model-1be44a13bd0d.stan', line 647, column 0 to column 44)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: gamma_lpdf: Random variable is inf, but must be positive finite! (in '/tmp/Rtmp9AUyja/model-1be44a13bd0d.stan', line 647, column 0 to column 44)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/Rtmp9AUyja/model-1be44a13bd0d.stan', line 652, column 0 to column 42)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: beta_proportion_lpdf: Location parameter is -nan, but must be positive! (in '/tmp/Rtmp9AUyja/model-1be44a13bd0d.stan', line 675, column 2 to line 677, column 4)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: beta_proportion_lpdf: Location parameter is -nan, but must be positive! (in '/tmp/Rtmp9AUyja/model-1be44a13bd0d.stan', line 675, column 2 to line 677, column 4)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: beta_proportion_lpdf: Location parameter is -nan, but must be positive! (in '/tmp/Rtmp9AUyja/model-1be44a13bd0d.stan', line 675, column 2 to line 677, column 4)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: gamma_lpdf: Random variable is inf, but must be positive finite! (in '/tmp/Rtmp9AUyja/model-1be44a13bd0d.stan', line 647, column 0 to column 44)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: gamma_lpdf: Random variable is inf, but must be positive finite! (in '/tmp/Rtmp9AUyja/model-1be44a13bd0d.stan', line 647, column 0 to column 44)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 2 Iteration:    1 / 1200 [  0%]  (Warmup) 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: gamma_lpdf: Random variable is inf, but must be positive finite! (in '/tmp/Rtmp9AUyja/model-1be44a13bd0d.stan', line 647, column 0 to column 44)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: gamma_lpdf: Random variable is inf, but must be positive finite! (in '/tmp/Rtmp9AUyja/model-1be44a13bd0d.stan', line 647, column 0 to column 44)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: gamma_lpdf: Random variable is inf, but must be positive finite! (in '/tmp/Rtmp9AUyja/model-1be44a13bd0d.stan', line 647, column 0 to column 44)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: beta_proportion_lpdf: Location parameter is -nan, but must be positive! (in '/tmp/Rtmp9AUyja/model-1be44a13bd0d.stan', line 675, column 2 to line 677, column 4)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: beta_proportion_lpdf: Location parameter is -nan, but must be positive! (in '/tmp/Rtmp9AUyja/model-1be44a13bd0d.stan', line 675, column 2 to line 677, column 4)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: beta_proportion_lpdf: Location parameter is -nan, but must be positive! (in '/tmp/Rtmp9AUyja/model-1be44a13bd0d.stan', line 675, column 2 to line 677, column 4)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: beta_proportion_lpdf: Location parameter is -nan, but must be positive! (in '/tmp/Rtmp9AUyja/model-1be44a13bd0d.stan', line 675, column 2 to line 677, column 4)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: gamma_lpdf: Random variable is inf, but must be positive finite! (in '/tmp/Rtmp9AUyja/model-1be44a13bd0d.stan', line 647, column 0 to column 44)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: gamma_lpdf: Random variable is inf, but must be positive finite! (in '/tmp/Rtmp9AUyja/model-1be44a13bd0d.stan', line 647, column 0 to column 44)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: beta_proportion_lpdf: Location parameter is -nan, but must be positive! (in '/tmp/Rtmp9AUyja/model-1be44a13bd0d.stan', line 675, column 2 to line 677, column 4)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 1 Iteration:  100 / 1200 [  8%]  (Warmup) 
#> Chain 2 Iteration:  100 / 1200 [  8%]  (Warmup) 
#> Chain 1 Iteration:  200 / 1200 [ 16%]  (Warmup) 
#> Chain 2 Iteration:  200 / 1200 [ 16%]  (Warmup) 
#> Chain 1 Iteration:  300 / 1200 [ 25%]  (Warmup) 
#> Chain 2 Iteration:  300 / 1200 [ 25%]  (Warmup) 
#> Chain 1 Iteration:  400 / 1200 [ 33%]  (Warmup) 
#> Chain 2 Iteration:  400 / 1200 [ 33%]  (Warmup) 
#> Chain 1 Iteration:  500 / 1200 [ 41%]  (Warmup) 
#> Chain 1 Iteration:  600 / 1200 [ 50%]  (Warmup) 
#> Chain 2 Iteration:  500 / 1200 [ 41%]  (Warmup) 
#> Chain 2 Iteration:  600 / 1200 [ 50%]  (Warmup) 
#> Chain 1 Iteration:  700 / 1200 [ 58%]  (Warmup) 
#> Chain 2 Iteration:  700 / 1200 [ 58%]  (Warmup) 
#> Chain 1 Iteration:  800 / 1200 [ 66%]  (Warmup) 
#> Chain 2 Iteration:  800 / 1200 [ 66%]  (Warmup) 
#> Chain 1 Iteration:  900 / 1200 [ 75%]  (Warmup) 
#> Chain 2 Iteration:  900 / 1200 [ 75%]  (Warmup) 
#> Chain 1 Iteration: 1000 / 1200 [ 83%]  (Warmup) 
#> Chain 1 Iteration: 1001 / 1200 [ 83%]  (Sampling) 
#> Chain 2 Iteration: 1000 / 1200 [ 83%]  (Warmup) 
#> Chain 2 Iteration: 1001 / 1200 [ 83%]  (Sampling) 
#> Chain 2 Iteration: 1100 / 1200 [ 91%]  (Sampling) 
#> Chain 1 Iteration: 1100 / 1200 [ 91%]  (Sampling) 
#> Chain 2 Iteration: 1200 / 1200 [100%]  (Sampling) 
#> Chain 2 finished in 14.7 seconds.
#> Chain 1 Iteration: 1200 / 1200 [100%]  (Sampling) 
#> Chain 1 finished in 15.8 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 15.3 seconds.
#> Total execution time: 16.0 seconds.
#> 
# }
```
