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

Functions to fit the unified multirates model.

- [`fit_multirates()`](https://caravagnalab.github.io/PEPI/reference/fit_multirates.md)
  : Fit the unified multirates model.

## Synthetic data

Functions to simulate synthetic driver/clade CCF trajectories for
testing and demonstration.

- [`simulate_multirates_tree()`](https://caravagnalab.github.io/PEPI/reference/simulate_multirates_tree.md)
  : Simulate synthetic input tables for \`init_multirates()\`.

## Getters

Functions to extract posterior draws and predicted quantities from a
fit.

- [`get_posterior_multirates()`](https://caravagnalab.github.io/PEPI/reference/get_posterior_multirates.md)
  : Extract and tidy the full posterior (and built-in prior draws) from
  a multirates fit.
- [`get_predicted_counts()`](https://caravagnalab.github.io/PEPI/reference/get_predicted_counts.md)
  : Get posterior-predictive counts (z\_\*\_pred) with observed counts
  joined in.
- [`get_predicted_fractions()`](https://caravagnalab.github.io/PEPI/reference/get_predicted_fractions.md)
  : Get posterior-predictive fractions (frac\_\*\_pred) with observed
  fractions joined in.
- [`get_predicted_ccf()`](https://caravagnalab.github.io/PEPI/reference/get_predicted_ccf.md)
  : Get posterior-predictive CCF (ccf\_\*\_pred) with observed CCF
  joined in.
- [`get_predicted_mutations()`](https://caravagnalab.github.io/PEPI/reference/get_predicted_mutations.md)
  : Get posterior-predictive mutation counts (m\_\*\_pred) with observed
  counts joined in.
- [`get_viber_clusters()`](https://caravagnalab.github.io/PEPI/reference/get_viber_clusters.md)
  : Clusters mutation with VIBER, a package that implements a
  variational non parametric Bayesian model to fit multi-variate
  Binomial mixtures
- [`get_init_values()`](https://caravagnalab.github.io/PEPI/reference/get_init_values.md)
  : Provide a initialization parameters for a PEPI VAF fit

## Plotting functions

Functions to plot data, predictions and posterior/prior distributions.

- [`plot_counts_multirates()`](https://caravagnalab.github.io/PEPI/reference/plot_counts_multirates.md)
  : Plot predicted vs observed counts for the unified multirates model.
- [`plot_fractions()`](https://caravagnalab.github.io/PEPI/reference/plot_fractions.md)
  : Plot predicted vs observed fractions for the unified multirates
  model.
- [`plot_ccf()`](https://caravagnalab.github.io/PEPI/reference/plot_ccf.md)
  : Plot predicted vs observed CCF for the unified multirates model.
- [`plot_mutations()`](https://caravagnalab.github.io/PEPI/reference/plot_mutations.md)
  : Plot predicted vs observed mutation counts for the unified
  multirates model.
- [`plot_inference()`](https://caravagnalab.github.io/PEPI/reference/plot_inference.md)
  : Plot posterior and prior distributions of the unified multirates
  model's parameters.
- [`plot_multivariate()`](https://caravagnalab.github.io/PEPI/reference/plot_multivariate.md)
  : plotting functions Plot multivariate VAF distributions with cluster
  associated to tree nodes.
- [`plot_marginal()`](https://caravagnalab.github.io/PEPI/reference/plot_marginal.md)
  : Plot marginal VAF distributions with cluster associated to tree
  nodes.

## Convergence diagnostics

Thin bayesplot wrappers for fits obtained with
`fit_multirates(..., method = "sample")`.

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
