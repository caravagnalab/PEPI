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
#> Chain 1 Exception: gamma_lpdf: Random variable is inf, but must be positive finite! (in '/tmp/RtmpKs1G8e/model-1d8638cbc177.stan', line 641, column 0 to column 44)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: gamma_lpdf: Random variable is inf, but must be positive finite! (in '/tmp/RtmpKs1G8e/model-1d8638cbc177.stan', line 641, column 0 to column 44)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: gamma_lpdf: Random variable is inf, but must be positive finite! (in '/tmp/RtmpKs1G8e/model-1d8638cbc177.stan', line 641, column 0 to column 44)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpKs1G8e/model-1d8638cbc177.stan', line 646, column 0 to column 42)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: beta_proportion_lpdf: Location parameter is -nan, but must be positive! (in '/tmp/RtmpKs1G8e/model-1d8638cbc177.stan', line 669, column 2 to line 671, column 4)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: beta_proportion_lpdf: Location parameter is -nan, but must be positive! (in '/tmp/RtmpKs1G8e/model-1d8638cbc177.stan', line 669, column 2 to line 671, column 4)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: beta_proportion_lpdf: Location parameter is -nan, but must be positive! (in '/tmp/RtmpKs1G8e/model-1d8638cbc177.stan', line 669, column 2 to line 671, column 4)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: gamma_lpdf: Random variable is inf, but must be positive finite! (in '/tmp/RtmpKs1G8e/model-1d8638cbc177.stan', line 641, column 0 to column 44)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: gamma_lpdf: Random variable is inf, but must be positive finite! (in '/tmp/RtmpKs1G8e/model-1d8638cbc177.stan', line 641, column 0 to column 44)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 2 Iteration:    1 / 1200 [  0%]  (Warmup) 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: gamma_lpdf: Random variable is inf, but must be positive finite! (in '/tmp/RtmpKs1G8e/model-1d8638cbc177.stan', line 641, column 0 to column 44)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: gamma_lpdf: Random variable is inf, but must be positive finite! (in '/tmp/RtmpKs1G8e/model-1d8638cbc177.stan', line 641, column 0 to column 44)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: gamma_lpdf: Random variable is inf, but must be positive finite! (in '/tmp/RtmpKs1G8e/model-1d8638cbc177.stan', line 641, column 0 to column 44)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: beta_proportion_lpdf: Location parameter is -nan, but must be positive! (in '/tmp/RtmpKs1G8e/model-1d8638cbc177.stan', line 669, column 2 to line 671, column 4)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: beta_proportion_lpdf: Location parameter is -nan, but must be positive! (in '/tmp/RtmpKs1G8e/model-1d8638cbc177.stan', line 669, column 2 to line 671, column 4)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: beta_proportion_lpdf: Location parameter is -nan, but must be positive! (in '/tmp/RtmpKs1G8e/model-1d8638cbc177.stan', line 669, column 2 to line 671, column 4)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: beta_proportion_lpdf: Location parameter is -nan, but must be positive! (in '/tmp/RtmpKs1G8e/model-1d8638cbc177.stan', line 669, column 2 to line 671, column 4)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: gamma_lpdf: Random variable is inf, but must be positive finite! (in '/tmp/RtmpKs1G8e/model-1d8638cbc177.stan', line 641, column 0 to column 44)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: gamma_lpdf: Random variable is inf, but must be positive finite! (in '/tmp/RtmpKs1G8e/model-1d8638cbc177.stan', line 641, column 0 to column 44)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: beta_proportion_lpdf: Location parameter is -nan, but must be positive! (in '/tmp/RtmpKs1G8e/model-1d8638cbc177.stan', line 669, column 2 to line 671, column 4)
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
#> Chain 2 finished in 19.6 seconds.
#> Chain 1 Iteration: 1200 / 1200 [100%]  (Sampling) 
#> Chain 1 finished in 19.8 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 19.7 seconds.
#> Total execution time: 20.0 seconds.
#> 
get_predicted_trajectories(x)
#> # A tibble: 38,400 × 10
#>     samp  time pop   id    z_neg   z_pos z_tot frac_pos ccf_neg ccf_pos
#>    <int> <dbl> <chr> <chr> <dbl>   <dbl> <dbl>    <dbl>   <dbl>   <dbl>
#>  1     1     0 total NA    7.03  0.0202  7.42   0.00286       1       1
#>  2     2     0 total NA    0.158 0.00188 0.165  0.0118        1       1
#>  3     3     0 total NA    2.52  0.00493 2.03   0.00195       1       1
#>  4     4     0 total NA    5.97  0.0375  6.17   0.00624       1       1
#>  5     5     0 total NA    1.01  0.00278 0.990  0.00276       1       1
#>  6     6     0 total NA    0.408 0.00408 0.394  0.00992       1       1
#>  7     7     0 total NA    2.09  0.0157  2.12   0.00742       1       1
#>  8     8     0 total NA    2.69  0.00745 2.66   0.00277       1       1
#>  9     9     0 total NA    0.920 0.00258 0.909  0.00280       1       1
#> 10    10     0 total NA    1.26  0.0110  1.42   0.00872       1       1
#> # ℹ 38,390 more rows
# }
```
