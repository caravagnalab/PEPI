# A short HMC fit of the polyclonal example shipped with PEPI (see the
# "Inference" article for the full fit).

source(system.file("examples", "tables_example.R", package = "PEPI"))

x = fit_multirates(x, seed = 31, chains = 2, parallel_chains = 2, ndraws = 200,
  init = rep(list(list(omega_p_wt = 3e-3, omega_p_driver = 3e-3)), 2),
  mu = mu$mean, l = 51e06, t_min = -20, ms_epi = -0.8, sigma_epi = 0.9,
  alpha_lambda = 10, beta_lambda = 10, alpha_n_wt = 1, beta_n_wt = 1e4,
  alpha_p_wt = 0.2, beta_p_wt = 10, min_kappa = 3000, max_kappa = 5000,
  min_sigma_count = 0.05, max_sigma_count = 0.15, include_poisson = 1, include_trunk = 0)
