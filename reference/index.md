# Package index

## Input & data construction

Functions to validate a multirates input and build the corresponding
Stan data list.

- [`init_multirates()`](https://caravagnalab.github.io/PEPI/reference/init_multirates.md)
  : Create a PEPI object for the unified multirates model.
- [`check_input_multirates()`](https://caravagnalab.github.io/PEPI/reference/check_input_multirates.md)
  : Check if multirates input data has the correct form.
- [`build_stan_data_multirates()`](https://caravagnalab.github.io/PEPI/reference/build_stan_data_multirates.md)
  : Build the Stan data list for the unified multirates model.

## Fitting

Fit the unified multirates model with Hamiltonian Monte Carlo.

- [`fit_multirates()`](https://caravagnalab.github.io/PEPI/reference/fit_multirates.md)
  : Fit the unified multirates model.

## Getters

Functions to extract posterior draws and posterior predictive
trajectories from a fit.

- [`get_posterior_multirates()`](https://caravagnalab.github.io/PEPI/reference/get_posterior_multirates.md)
  : Extract and tidy the full posterior (and built-in prior draws) from
  a multirates fit.
- [`get_predicted_trajectories()`](https://caravagnalab.github.io/PEPI/reference/get_predicted_trajectories.md)
  : Get posterior predictive trajectories of the sorted populations.
- [`get_init_values()`](https://caravagnalab.github.io/PEPI/reference/get_init_values.md)
  : Provide a initialization parameters for a PEPI VAF fit

## Inference plots

Prior vs posterior distributions and posterior predictive checks.

- [`plot_prior_posterior()`](https://caravagnalab.github.io/PEPI/reference/plot_prior_posterior.md)
  : plotting functions Plot prior and posterior distributions of the
  model parameters.
- [`plot_posterior_predictive()`](https://caravagnalab.github.io/PEPI/reference/plot_posterior_predictive.md)
  : Plot all posterior predictive checks.
- [`plot_fractions()`](https://caravagnalab.github.io/PEPI/reference/plot_fractions.md)
  : Plot the posterior predictive fraction of "+" cells in the whole
  population.
- [`plot_counts()`](https://caravagnalab.github.io/PEPI/reference/plot_counts.md)
  : Plot the posterior predictive total population size.
- [`plot_split()`](https://caravagnalab.github.io/PEPI/reference/plot_split.md)
  : Plot the posterior predictive "+"/"-" split of each population.
- [`plot_ccf()`](https://caravagnalab.github.io/PEPI/reference/plot_ccf.md)
  : Plot the posterior predictive CCF of each population in the sorted
  samples.

## Convergence diagnostics

Thin bayesplot wrappers for the MCMC fit returned by
[`fit_multirates()`](https://caravagnalab.github.io/PEPI/reference/fit_multirates.md).

- [`plot_diagnostics()`](https://caravagnalab.github.io/PEPI/reference/plot_diagnostics.md)
  : One-call convergence dashboard.
- [`plot_trace()`](https://caravagnalab.github.io/PEPI/reference/plot_trace.md)
  : Trace plots of the sampled chains.
- [`plot_rhat()`](https://caravagnalab.github.io/PEPI/reference/plot_rhat.md)
  : Rhat (potential scale reduction factor) overview.
- [`plot_ess()`](https://caravagnalab.github.io/PEPI/reference/plot_ess.md)
  : Effective sample size (ratio to total draws) overview.
- [`plot_acf()`](https://caravagnalab.github.io/PEPI/reference/plot_acf.md)
  : Autocorrelation of the sampled chains.
- [`plot_divergences()`](https://caravagnalab.github.io/PEPI/reference/plot_divergences.md)
  : Divergent transitions.
- [`plot_energy()`](https://caravagnalab.github.io/PEPI/reference/plot_energy.md)
  : Hamiltonian Monte Carlo energy diagnostic.
- [`plot_treedepth()`](https://caravagnalab.github.io/PEPI/reference/plot_treedepth.md)
  : NUTS treedepth diagnostic.
- [`plot_pairs()`](https://caravagnalab.github.io/PEPI/reference/plot_pairs.md)
  : Pairs plot with divergences marked.
