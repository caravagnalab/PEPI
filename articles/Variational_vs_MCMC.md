# Variational vs MCMC inference

``` r

library(PEPI)
library(cmdstanr)
#> This is cmdstanr version 0.8.0
#> - CmdStanR documentation and vignettes: mc-stan.org/cmdstanr
#> - CmdStan path: /home/runner/.cmdstan/cmdstan-2.39.0
#> - CmdStan version: 2.39.0
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(ggplot2)
```

[`fit_multirates()`](https://caravagnalab.github.io/PEPI/reference/fit_multirates.md)
supports two inference methods:

- `method = "variational"` (the default) runs ADVI, cmdstanr’s automatic
  differentiation variational inference. It is fast and a good choice
  while exploring data or iterating on priors, but it approximates the
  posterior with a (full-rank Gaussian) parametric family, has no
  convergence diagnostics, and can under- or over-estimate uncertainty
  in ways that are hard to detect from the fit alone.
- `method = "sample"` runs full Hamiltonian Monte Carlo (NUTS). It is
  slower but asymptotically exact, and comes with the convergence
  diagnostics covered in `vignette("Convergence_diagnostics")` (R-hat,
  effective sample size, divergences, …).

This vignette fits the same data both ways and compares the two
posteriors directly, to build intuition for when the difference matters.

## Simulate data and fit both ways

Stan’s own sampler console output is suppressed below for readability -
it normally includes informational messages about rejected initial
values during warmup, which are expected and not a sign of a problem
here.

``` r


set.seed(42)

sim = simulate_multirates_tree(
  sampling_times = c(3, 6, 9),
  tmrca = 1, t_min = 0,
  lambda_n = 1, s_epi = 0.15,
  omega_n_wt = 5e-3, omega_p_wt = 2e-3,
  drivers = list(
    list(id = "KRAS", t_driver = 2, s_driver = 0.3,
        omega_n_driver = 5e-3, omega_p_driver = 2e-3,
        ms_driver = -0.5, sigma_driver = 0.5,
        alpha_n_driver = 1, beta_n_driver = 10,
        alpha_p_driver = 1, beta_p_driver = 10)
  ),
  clades_wt = list(
    list(id = "c1", t_clade_wt = 1.5)
  ),
  mu = 1e-7, l = 2.7e9, kappa = 20, sigma_count = 0.1,
  seed = 42
)

fit_args = list(
  cmdstan_path = cmdstanr::cmdstan_path(), seed = 45,
  mu = 1e-7, l = 2.7e9, t_min = 0, ms_epi = 0, sigma_epi = 0.5,
  alpha_lambda = 1, beta_lambda = 1, alpha_n_wt = 1, beta_n_wt = 10,
  alpha_p_wt = 1, beta_p_wt = 10, min_kappa = 5, max_kappa = 200,
  min_sigma_count = 0.01, max_sigma_count = 1
)

x_vb = init_multirates(sim$tables, m_trunk = sim$truth$m_trunk)
x_vb = do.call(fit_multirates, c(list(x = x_vb, method = "variational", ndraws = 1000), fit_args))
x_vb = get_posterior_multirates(x_vb)

x_mc = init_multirates(sim$tables, m_trunk = sim$truth$m_trunk)
x_mc = do.call(fit_multirates, c(list(x = x_mc, method = "sample", chains = 2, ndraws = 200), fit_args))
x_mc = get_posterior_multirates(x_mc)
```

## Compare the two posteriors

``` r


params = c("lambda_n", "s_epi", "s_driver", "omega_n_wt", "omega_p_wt")

comparison = bind_rows(
  x_vb$posterior$multirates %>% filter(family == "param", type == "posterior") %>% mutate(method = "variational"),
  x_mc$posterior$multirates %>% filter(family == "param", type == "posterior") %>% mutate(method = "sample")
) %>%
  filter(base %in% params) %>%
  mutate(facet = ifelse(is.na(id), base, paste0(base, "[", id, "]")))

ggplot(comparison) +
  geom_density(aes(x = value, fill = method), alpha = 0.5) +
  facet_wrap(~facet, scales = "free") +
  theme_light(base_size = 10) +
  theme(legend.position = "bottom")
```

![](Variational_vs_MCMC_files/figure-html/unnamed-chunk-3-1.png)

For `lambda_n`, `omega_n_wt`, `s_driver` and `s_epi` the two methods
agree closely. `omega_p_wt` is a different story: ADVI’s Gaussian
approximation gives it a longer right tail than MCMC does, because this
switch rate’s posterior is concentrated near a boundary of its support -
exactly the kind of skewed, non-Gaussian shape ADVI struggles to
approximate well. Had we only fit with `method = "variational"`, we
would have overstated the uncertainty on `omega_p_wt` without any way to
notice from the ADVI fit alone.

## Confirm the MCMC fit is trustworthy before trusting the comparison

The comparison above is only informative if the MCMC fit itself has
converged - otherwise “MCMC disagrees with ADVI” could just as easily
mean “the sampler hasn’t mixed yet”. Always check:

``` r

plot_diagnostics(x_mc)
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> `stat_bin()` using `bins = 30`. Pick better value `binwidth`.
#> `stat_bin()` using `bins = 30`. Pick better value `binwidth`.
```

![](Variational_vs_MCMC_files/figure-html/unnamed-chunk-4-1.png)

## Practical guidance

- Use `method = "variational"` while exploring data, iterating on
  priors, or checking that the model runs at all - it is typically an
  order of magnitude faster.
- Before reporting final results, refit the same data with
  `method = "sample"` and inspect
  [`plot_diagnostics()`](https://caravagnalab.github.io/PEPI/reference/plot_diagnostics.md).
- If the two posteriors disagree noticeably for a parameter, trust the
  MCMC fit (once it has converged) - ADVI’s approximation error is
  exactly what full sampling is designed to avoid.
