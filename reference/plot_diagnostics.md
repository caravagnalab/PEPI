# One-call convergence dashboard.

Arranges Rhat, effective sample size, divergences and the energy
diagnostic into a single 2x2 figure - a quick first look at whether an
MCMC fit is trustworthy.

## Usage

``` r
plot_diagnostics(x, params = NULL)
```

## Arguments

- x:

  PEPI_Multirates object fit with `method = "sample"`.

- params:

  Vector of Stan parameter base names for the Rhat/ESS panels. Defaults
  to every sampled parameter.

## Value

A combined plot (see
[`ggpubr::ggarrange()`](https://rpkgs.datanovia.com/ggpubr/reference/ggarrange.html)).

## Examples

``` r
# \donttest{
sim = simulate_multirates_tree(
  sampling_times = c(3, 6, 9), tmrca = 1, t_min = 0,
  lambda_n = 1, s_epi = 0.15, omega_n_wt = 5e-3, omega_p_wt = 2e-3,
  drivers = list(list(id = "KRAS", t_driver = 2, s_driver = 0.3,
      omega_n_driver = 5e-3, omega_p_driver = 2e-3,
      ms_driver = -0.5, sigma_driver = 0.5,
      alpha_n_driver = 1, beta_n_driver = 10, alpha_p_driver = 1, beta_p_driver = 10)),
  clades_wt = list(list(id = "c1", t_clade_wt = 1.5)),
  mu = 1e-7, l = 2.7e9, kappa = 20, sigma_count = 0.1, seed = 42
)
x = init_multirates(sim$tables, m_trunk = sim$truth$m_trunk)
x = fit_multirates(x, cmdstan_path = cmdstanr::cmdstan_path(),
  method = "sample", chains = 2, ndraws = 200, seed = 45,
  mu = 1e-7, l = 2.7e9, t_min = 0, ms_epi = 0, sigma_epi = 0.5,
  alpha_lambda = 1, beta_lambda = 1, alpha_n_wt = 1, beta_n_wt = 10,
  alpha_p_wt = 1, beta_p_wt = 10, min_kappa = 5, max_kappa = 200,
  min_sigma_count = 0.01, max_sigma_count = 1)
#> CmdStan path set to: /home/runner/.cmdstan/cmdstan-2.39.0
#> Init values were only set for a subset of parameters. 
#> Missing init values for the following parameters:
#>  - chain 1: lambda_n, s_epi, omega_p_wt, omega_n_wt, bc_wt, U_wt, W_wt, bc_clade_wt, s_driver, bc_driver, U_driver, W_driver, s_driver_n, t_driver_n, bc_driver_n, s_driver_p, t_driver_p, bc_driver_p, s_dc, delta_t_dc, bc_driver_dc, bc_clade_dc, U_dc, W_dc, omega_p_dc, omega_n_dc, omega_p_driver, omega_n_driver, kappa, sigma_count
#>  - chain 2: lambda_n, s_epi, omega_p_wt, omega_n_wt, bc_wt, U_wt, W_wt, bc_clade_wt, s_driver, bc_driver, U_driver, W_driver, s_driver_n, t_driver_n, bc_driver_n, s_driver_p, t_driver_p, bc_driver_p, s_dc, delta_t_dc, bc_driver_dc, bc_clade_dc, U_dc, W_dc, omega_p_dc, omega_n_dc, omega_p_driver, omega_n_driver, kappa, sigma_count
#> 
#> To disable this message use options(cmdstanr_warn_inits = FALSE).
#> Running MCMC with 2 sequential chains...
#> 
#> Chain 1 Rejecting initial value:
#> Chain 1   Error evaluating the log probability at the initial value.
#> Chain 1 Exception: lub_constrain: lb is 4.87392, but must be less than 3.000000 (in '/tmp/RtmpAjDVJG/model-1be660c054e.stan', line 353, column 2 to column 81)
#> Chain 1 Exception: lub_constrain: lb is 4.87392, but must be less than 3.000000 (in '/tmp/RtmpAjDVJG/model-1be660c054e.stan', line 353, column 2 to column 81)
#> Chain 1 Rejecting initial value:
#> Chain 1   Error evaluating the log probability at the initial value.
#> Chain 1 Exception: lub_constrain: lb is 5.59741, but must be less than 3.000000 (in '/tmp/RtmpAjDVJG/model-1be660c054e.stan', line 353, column 2 to column 81)
#> Chain 1 Exception: lub_constrain: lb is 5.59741, but must be less than 3.000000 (in '/tmp/RtmpAjDVJG/model-1be660c054e.stan', line 353, column 2 to column 81)
#> Chain 1 Rejecting initial value:
#> Chain 1   Error evaluating the log probability at the initial value.
#> Chain 1 Exception: lub_constrain: lb is 6.88754, but must be less than 3.000000 (in '/tmp/RtmpAjDVJG/model-1be660c054e.stan', line 353, column 2 to column 81)
#> Chain 1 Exception: lub_constrain: lb is 6.88754, but must be less than 3.000000 (in '/tmp/RtmpAjDVJG/model-1be660c054e.stan', line 353, column 2 to column 81)
#> Chain 1 Iteration:    1 / 1200 [  0%]  (Warmup) 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: lub_constrain: lb is 9, but must be less than 3.000000 (in '/tmp/RtmpAjDVJG/model-1be660c054e.stan', line 353, column 2 to column 81)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: lub_constrain: lb is 9, but must be less than 3.000000 (in '/tmp/RtmpAjDVJG/model-1be660c054e.stan', line 353, column 2 to column 81)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: lub_constrain: lb is 9, but must be less than 3.000000 (in '/tmp/RtmpAjDVJG/model-1be660c054e.stan', line 353, column 2 to column 81)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: lub_constrain: lb is 9, but must be less than 3.000000 (in '/tmp/RtmpAjDVJG/model-1be660c054e.stan', line 353, column 2 to column 81)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: lub_constrain: lb is 9, but must be less than 3.000000 (in '/tmp/RtmpAjDVJG/model-1be660c054e.stan', line 353, column 2 to column 81)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: lub_constrain: lb is 8.99695, but must be less than 3.000000 (in '/tmp/RtmpAjDVJG/model-1be660c054e.stan', line 353, column 2 to column 81)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: lub_constrain: lb is 5.17598, but must be less than 3.000000 (in '/tmp/RtmpAjDVJG/model-1be660c054e.stan', line 353, column 2 to column 81)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: lub_constrain: lb is 9, but must be less than 3.000000 (in '/tmp/RtmpAjDVJG/model-1be660c054e.stan', line 353, column 2 to column 81)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: lub_constrain: lb is 9, but must be less than 3.000000 (in '/tmp/RtmpAjDVJG/model-1be660c054e.stan', line 353, column 2 to column 81)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: lub_constrain: lb is 9, but must be less than 3.000000 (in '/tmp/RtmpAjDVJG/model-1be660c054e.stan', line 353, column 2 to column 81)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: lub_constrain: lb is 3.00815, but must be less than 3.000000 (in '/tmp/RtmpAjDVJG/model-1be660c054e.stan', line 353, column 2 to column 81)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Iteration:  100 / 1200 [  8%]  (Warmup) 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: lub_constrain: lb is 3.29367, but must be less than 3.000000 (in '/tmp/RtmpAjDVJG/model-1be660c054e.stan', line 353, column 2 to column 81)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Iteration:  200 / 1200 [ 16%]  (Warmup) 
#> Chain 1 Iteration:  300 / 1200 [ 25%]  (Warmup) 
#> Chain 1 Iteration:  400 / 1200 [ 33%]  (Warmup) 
#> Chain 1 Iteration:  500 / 1200 [ 41%]  (Warmup) 
#> Chain 1 Iteration:  600 / 1200 [ 50%]  (Warmup) 
#> Chain 1 Iteration:  700 / 1200 [ 58%]  (Warmup) 
#> Chain 1 Iteration:  800 / 1200 [ 66%]  (Warmup) 
#> Chain 1 Iteration:  900 / 1200 [ 75%]  (Warmup) 
#> Chain 1 Iteration: 1000 / 1200 [ 83%]  (Warmup) 
#> Chain 1 Iteration: 1001 / 1200 [ 83%]  (Sampling) 
#> Chain 1 Iteration: 1100 / 1200 [ 91%]  (Sampling) 
#> Chain 1 Iteration: 1200 / 1200 [100%]  (Sampling) 
#> Chain 1 finished in 1.2 seconds.
#> Chain 2 Iteration:    1 / 1200 [  0%]  (Warmup) 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: lub_constrain: lb is 9, but must be less than 3.000000 (in '/tmp/RtmpAjDVJG/model-1be660c054e.stan', line 353, column 2 to column 81)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: lub_constrain: lb is 9, but must be less than 3.000000 (in '/tmp/RtmpAjDVJG/model-1be660c054e.stan', line 353, column 2 to column 81)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: lub_constrain: lb is 9, but must be less than 3.000000 (in '/tmp/RtmpAjDVJG/model-1be660c054e.stan', line 353, column 2 to column 81)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: lub_constrain: lb is 9, but must be less than 3.000000 (in '/tmp/RtmpAjDVJG/model-1be660c054e.stan', line 353, column 2 to column 81)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: lub_constrain: lb is 8.99953, but must be less than 3.000000 (in '/tmp/RtmpAjDVJG/model-1be660c054e.stan', line 353, column 2 to column 81)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: lub_constrain: lb is 6.51592, but must be less than 3.000000 (in '/tmp/RtmpAjDVJG/model-1be660c054e.stan', line 353, column 2 to column 81)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: lub_constrain: lb is 9, but must be less than 3.000000 (in '/tmp/RtmpAjDVJG/model-1be660c054e.stan', line 353, column 2 to column 81)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: lub_constrain: lb is 9, but must be less than 3.000000 (in '/tmp/RtmpAjDVJG/model-1be660c054e.stan', line 353, column 2 to column 81)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: lub_constrain: lb is 5.75093, but must be less than 3.000000 (in '/tmp/RtmpAjDVJG/model-1be660c054e.stan', line 353, column 2 to column 81)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Iteration:  100 / 1200 [  8%]  (Warmup) 
#> Chain 2 Iteration:  200 / 1200 [ 16%]  (Warmup) 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: lub_constrain: lb is 3.23029, but must be less than 3.000000 (in '/tmp/RtmpAjDVJG/model-1be660c054e.stan', line 353, column 2 to column 81)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Iteration:  300 / 1200 [ 25%]  (Warmup) 
#> Chain 2 Iteration:  400 / 1200 [ 33%]  (Warmup) 
#> Chain 2 Iteration:  500 / 1200 [ 41%]  (Warmup) 
#> Chain 2 Iteration:  600 / 1200 [ 50%]  (Warmup) 
#> Chain 2 Iteration:  700 / 1200 [ 58%]  (Warmup) 
#> Chain 2 Iteration:  800 / 1200 [ 66%]  (Warmup) 
#> Chain 2 Iteration:  900 / 1200 [ 75%]  (Warmup) 
#> Chain 2 Iteration: 1000 / 1200 [ 83%]  (Warmup) 
#> Chain 2 Iteration: 1001 / 1200 [ 83%]  (Sampling) 
#> Chain 2 Iteration: 1100 / 1200 [ 91%]  (Sampling) 
#> Chain 2 Iteration: 1200 / 1200 [100%]  (Sampling) 
#> Chain 2 finished in 1.3 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 1.2 seconds.
#> Total execution time: 2.6 seconds.
#> 
#> Warning: 1 of 400 (0.0%) transitions ended with a divergence.
#> See https://mc-stan.org/misc/warnings for details.
plot_diagnostics(x)
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> `stat_bin()` using `bins = 30`. Pick better value `binwidth`.
#> `stat_bin()` using `bins = 30`. Pick better value `binwidth`.

# }
```
