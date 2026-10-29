# Get posterior predictive trajectories of the sorted populations.

For every posterior draw and every time of a grid, the latent "-" and
"+" sizes of the wild type and of every driver (and dc/driver-only)
clade are computed with the model's growth/switch approximation.
Observations are then generated as in the likelihood of the Stan model:
for each population, its total count is drawn from \`lognormal(log(Z- +
Z+), sigma_count)\` and its split in "+" cells from \`beta_proportion(Z+
/ (Z- + Z+), kappa)\`. The private populations driver_n and driver_p
live in one compartment and only get the lognormal count, as do wt/dc
clades whose observed "-" count is 0. CCFs are evaluated on the sampled
counts, within the "-" and "+" samples; wt clades and dc clades are
subpopulations of wt and dc and do not enter the sample sizes. The total
population size (\`z_tot\` of \`pop = "total"\`) follows the likelihood
of the intermediate total counts, \`lognormal(log(ztot_lat),
sigma_count)\`, with \`ztot_lat\` the latent size of all populations
except the clades.

## Usage

``` r
get_predicted_trajectories(
  x,
  times = NULL,
  n_grid = 30,
  t_start = 0,
  noise = TRUE
)
```

## Arguments

- x:

  PEPI_Multirates object fit with
  [`fit_multirates()`](https://caravagnalab.github.io/PEPI/reference/fit_multirates.md).

- times:

  Time grid. If NULL, \`n_grid\` equally spaced times from \`t_start\`
  to the last observed time + 0.25, plus the observed times.

- n_grid:

  Number of points of the default time grid.

- t_start:

  First time of the default time grid.

- noise:

  If TRUE (default), sample the observations; if FALSE, return the
  latent trajectories.

## Value

A tibble with one row per draw (\`samp\`), \`time\` and population
(\`pop\`, \`id\`; \`pop = "total"\` is the whole population) and columns
\`z_neg\`, \`z_pos\` (counts), \`z_tot\` (total size, only for
\`"total"\`), \`frac_pos\` (fraction of "+" cells within the
population), \`ccf_neg\`, \`ccf_pos\` (CCF in the "-" and "+" samples).

## Examples

``` r
# \donttest{
source(system.file("examples", "fit_example.R", package = "PEPI"))
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
#> Chain 1 Iteration:  300 / 1200 [ 25%]  (Warmup) 
#> Chain 2 Iteration:  200 / 1200 [ 16%]  (Warmup) 
#> Chain 1 Iteration:  400 / 1200 [ 33%]  (Warmup) 
#> Chain 2 Iteration:  300 / 1200 [ 25%]  (Warmup) 
#> Chain 2 Iteration:  400 / 1200 [ 33%]  (Warmup) 
#> Chain 1 Iteration:  500 / 1200 [ 41%]  (Warmup) 
#> Chain 1 Iteration:  600 / 1200 [ 50%]  (Warmup) 
#> Chain 2 Iteration:  500 / 1200 [ 41%]  (Warmup) 
#> Chain 1 Iteration:  700 / 1200 [ 58%]  (Warmup) 
#> Chain 2 Iteration:  600 / 1200 [ 50%]  (Warmup) 
#> Chain 1 Iteration:  800 / 1200 [ 66%]  (Warmup) 
#> Chain 2 Iteration:  700 / 1200 [ 58%]  (Warmup) 
#> Chain 1 Iteration:  900 / 1200 [ 75%]  (Warmup) 
#> Chain 2 Iteration:  800 / 1200 [ 66%]  (Warmup) 
#> Chain 2 Iteration:  900 / 1200 [ 75%]  (Warmup) 
#> Chain 1 Iteration: 1000 / 1200 [ 83%]  (Warmup) 
#> Chain 1 Iteration: 1001 / 1200 [ 83%]  (Sampling) 
#> Chain 2 Iteration: 1000 / 1200 [ 83%]  (Warmup) 
#> Chain 2 Iteration: 1001 / 1200 [ 83%]  (Sampling) 
#> Chain 1 Iteration: 1100 / 1200 [ 91%]  (Sampling) 
#> Chain 2 Iteration: 1100 / 1200 [ 91%]  (Sampling) 
#> Chain 2 Iteration: 1200 / 1200 [100%]  (Sampling) 
#> Chain 2 finished in 15.3 seconds.
#> Chain 1 Iteration: 1200 / 1200 [100%]  (Sampling) 
#> Chain 1 finished in 15.6 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 15.4 seconds.
#> Total execution time: 15.6 seconds.
#> 
#> Error in read_cmdstan_csv(files = self$output_files(include_failed = FALSE),     variables = variables, sampler_diagnostics = sampler_diagnostics,     format = format): Assertion on 'files' failed: File does not exist: '/tmp/Rtmp9AUyja/multirates_positive_s-202610291009-1-8ea04a.csv'.
get_predicted_trajectories(x)
#> Error in get_predicted_trajectories(x): no multirates inference
# }
```
