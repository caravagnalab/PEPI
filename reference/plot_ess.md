# Effective sample size (ratio to total draws) overview.

One point per parameter; low ratios indicate highly autocorrelated
chains that need more draws for the same effective precision.

## Usage

``` r
plot_ess(x, params = NULL)
```

## Arguments

- x:

  PEPI_Multirates object fit with
  [`fit_multirates()`](https://caravagnalab.github.io/PEPI/reference/fit_multirates.md).

- params:

  Vector of Stan parameter base names to include. Defaults to every
  sampled parameter (excludes internal latent-state arrays and generated
  quantities).

## Value

A ggplot object (see
[`bayesplot::mcmc_neff()`](https://mc-stan.org/bayesplot/reference/MCMC-diagnostics.html)).

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
#> Chain 2 finished in 18.6 seconds.
#> Chain 1 Iteration: 1200 / 1200 [100%]  (Sampling) 
#> Chain 1 finished in 19.7 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 19.1 seconds.
#> Total execution time: 19.8 seconds.
#> 
#> Error in read_cmdstan_csv(files = self$output_files(include_failed = FALSE),     variables = variables, sampler_diagnostics = sampler_diagnostics,     format = format): Assertion on 'files' failed: File does not exist: '/tmp/RtmpKs1G8e/multirates_positive_s-202610201935-1-8ea1ec.csv'.
plot_ess(x)
#> Error in .multirates_mcmc_fit(x): no multirates inference
# }
```
