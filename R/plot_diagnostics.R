# MCMC convergence diagnostics for the unified multirates model
#
# Thin wrappers around bayesplot, which has built-in support for cmdstanr's
# CmdStanMCMC fit objects (bayesplot::rhat(), neff_ratio(), nuts_params(),
# log_posterior() all have CmdStanMCMC methods) - the diagnostics below just
# extract the right pieces and hand them to the corresponding mcmc_*() plot.
# These are only meaningful for a fit obtained with
# fit_multirates(..., method = "sample"); the default method = "variational"
# runs ADVI, which has no chains, Rhat, divergences, or treedepth to speak of.

.multirates_default_diag_params = c("lambda_n","s_epi","omega_n_wt","omega_p_wt","tmrca","kappa","sigma_count")

#' @keywords internal
.multirates_sampled_param_names = function(){

  bases = names(.multirates_dims)
  bases[vapply(bases, .multirates_family, character(1)) == "param" & !grepl("_prior$", bases)]

}

#' @keywords internal
.multirates_mcmc_fit = function(x){

  fit = x$inference$multirates

  if(is.null(fit)){
    stop("no multirates inference")
  }

  if(!inherits(fit, "CmdStanMCMC")){
    stop("convergence diagnostics require a fit obtained with fit_multirates(..., method = 'sample'); ",
         "the default method = 'variational' (ADVI) has no chains, Rhat, divergences, or treedepth.")
  }

  fit

}

#' @keywords internal
.multirates_available_params = function(fit, params){

  intersect(params, fit$metadata()$stan_variables)

}

#' Trace plots of the sampled chains.
#'
#' @param x PEPI_Multirates object fit with \code{method = "sample"}.
#' @param params Vector of Stan parameter base names (no brackets - all
#'   indexed elements of a vector/array parameter are included automatically).
#'   Defaults to a handful of global scalar parameters.
#' @param divergences If TRUE (default), divergent transitions are marked on the trace.
#' @return A ggplot object (see \code{bayesplot::mcmc_trace()}).
#' @examples
#' \donttest{
#' sim = simulate_multirates_tree(
#'   sampling_times = c(3, 6, 9), tmrca = 1, t_min = 0,
#'   lambda_n = 1, s_epi = 0.15, omega_n_wt = 5e-3, omega_p_wt = 2e-3,
#'   drivers = list(list(id = "KRAS", t_driver = 2, s_driver = 0.3,
#'       omega_n_driver = 5e-3, omega_p_driver = 2e-3,
#'       ms_driver = -0.5, sigma_driver = 0.5,
#'       alpha_n_driver = 1, beta_n_driver = 10, alpha_p_driver = 1, beta_p_driver = 10)),
#'   clades_wt = list(list(id = "c1", t_clade_wt = 1.5)),
#'   mu = 1e-7, l = 2.7e9, kappa = 20, sigma_count = 0.1, seed = 42
#' )
#' x = init_multirates(sim$tables, m_trunk = sim$truth$m_trunk)
#' x = fit_multirates(x, cmdstan_path = cmdstanr::cmdstan_path(),
#'   method = "sample", chains = 2, ndraws = 200, seed = 45,
#'   mu = 1e-7, l = 2.7e9, t_min = 0, ms_epi = 0, sigma_epi = 0.5,
#'   alpha_lambda = 1, beta_lambda = 1, alpha_n_wt = 1, beta_n_wt = 10,
#'   alpha_p_wt = 1, beta_p_wt = 10, min_kappa = 5, max_kappa = 200,
#'   min_sigma_count = 0.01, max_sigma_count = 1)
#' plot_trace(x, params = c("lambda_n","s_epi"))
#' }
#' @export

plot_trace = function(x, params = NULL, divergences = TRUE){

  fit = .multirates_mcmc_fit(x)

  if(is.null(params)) params = .multirates_default_diag_params
  params = .multirates_available_params(fit, params)

  np = if(divergences) bayesplot::nuts_params(fit) else NULL

  bayesplot::mcmc_trace(fit$draws(variables = params), np = np)

}

#' Rhat (potential scale reduction factor) overview.
#'
#' One point per parameter; values close to 1 indicate the chains have
#' converged to a common distribution. Values above 1.05 are a red flag.
#'
#' @param x PEPI_Multirates object fit with \code{method = "sample"}.
#' @param params Vector of Stan parameter base names to include. Defaults
#'   to every sampled parameter (excludes internal latent-state arrays and
#'   generated quantities).
#' @return A ggplot object (see \code{bayesplot::mcmc_rhat()}).
#' @examples
#' \donttest{
#' sim = simulate_multirates_tree(
#'   sampling_times = c(3, 6, 9), tmrca = 1, t_min = 0,
#'   lambda_n = 1, s_epi = 0.15, omega_n_wt = 5e-3, omega_p_wt = 2e-3,
#'   drivers = list(list(id = "KRAS", t_driver = 2, s_driver = 0.3,
#'       omega_n_driver = 5e-3, omega_p_driver = 2e-3,
#'       ms_driver = -0.5, sigma_driver = 0.5,
#'       alpha_n_driver = 1, beta_n_driver = 10, alpha_p_driver = 1, beta_p_driver = 10)),
#'   clades_wt = list(list(id = "c1", t_clade_wt = 1.5)),
#'   mu = 1e-7, l = 2.7e9, kappa = 20, sigma_count = 0.1, seed = 42
#' )
#' x = init_multirates(sim$tables, m_trunk = sim$truth$m_trunk)
#' x = fit_multirates(x, cmdstan_path = cmdstanr::cmdstan_path(),
#'   method = "sample", chains = 2, ndraws = 200, seed = 45,
#'   mu = 1e-7, l = 2.7e9, t_min = 0, ms_epi = 0, sigma_epi = 0.5,
#'   alpha_lambda = 1, beta_lambda = 1, alpha_n_wt = 1, beta_n_wt = 10,
#'   alpha_p_wt = 1, beta_p_wt = 10, min_kappa = 5, max_kappa = 200,
#'   min_sigma_count = 0.01, max_sigma_count = 1)
#' plot_rhat(x)
#' }
#' @export

plot_rhat = function(x, params = NULL){

  fit = .multirates_mcmc_fit(x)

  if(is.null(params)) params = .multirates_sampled_param_names()
  params = .multirates_available_params(fit, params)

  bayesplot::mcmc_rhat(bayesplot::rhat(fit, pars = params))

}

#' Effective sample size (ratio to total draws) overview.
#'
#' One point per parameter; low ratios indicate highly autocorrelated
#' chains that need more draws for the same effective precision.
#'
#' @param x PEPI_Multirates object fit with \code{method = "sample"}.
#' @param params Vector of Stan parameter base names to include. Defaults
#'   to every sampled parameter (excludes internal latent-state arrays and
#'   generated quantities).
#' @return A ggplot object (see \code{bayesplot::mcmc_neff()}).
#' @examples
#' \donttest{
#' sim = simulate_multirates_tree(
#'   sampling_times = c(3, 6, 9), tmrca = 1, t_min = 0,
#'   lambda_n = 1, s_epi = 0.15, omega_n_wt = 5e-3, omega_p_wt = 2e-3,
#'   drivers = list(list(id = "KRAS", t_driver = 2, s_driver = 0.3,
#'       omega_n_driver = 5e-3, omega_p_driver = 2e-3,
#'       ms_driver = -0.5, sigma_driver = 0.5,
#'       alpha_n_driver = 1, beta_n_driver = 10, alpha_p_driver = 1, beta_p_driver = 10)),
#'   clades_wt = list(list(id = "c1", t_clade_wt = 1.5)),
#'   mu = 1e-7, l = 2.7e9, kappa = 20, sigma_count = 0.1, seed = 42
#' )
#' x = init_multirates(sim$tables, m_trunk = sim$truth$m_trunk)
#' x = fit_multirates(x, cmdstan_path = cmdstanr::cmdstan_path(),
#'   method = "sample", chains = 2, ndraws = 200, seed = 45,
#'   mu = 1e-7, l = 2.7e9, t_min = 0, ms_epi = 0, sigma_epi = 0.5,
#'   alpha_lambda = 1, beta_lambda = 1, alpha_n_wt = 1, beta_n_wt = 10,
#'   alpha_p_wt = 1, beta_p_wt = 10, min_kappa = 5, max_kappa = 200,
#'   min_sigma_count = 0.01, max_sigma_count = 1)
#' plot_ess(x)
#' }
#' @export

plot_ess = function(x, params = NULL){

  fit = .multirates_mcmc_fit(x)

  if(is.null(params)) params = .multirates_sampled_param_names()
  params = .multirates_available_params(fit, params)

  bayesplot::mcmc_neff(bayesplot::neff_ratio(fit, pars = params))

}

#' Autocorrelation of the sampled chains.
#'
#' @param x PEPI_Multirates object fit with \code{method = "sample"}.
#' @param params Vector of Stan parameter base names. Defaults to a
#'   handful of global scalar parameters.
#' @param lags Number of lags to show.
#' @return A ggplot object (see \code{bayesplot::mcmc_acf()}).
#' @examples
#' \donttest{
#' sim = simulate_multirates_tree(
#'   sampling_times = c(3, 6, 9), tmrca = 1, t_min = 0,
#'   lambda_n = 1, s_epi = 0.15, omega_n_wt = 5e-3, omega_p_wt = 2e-3,
#'   drivers = list(list(id = "KRAS", t_driver = 2, s_driver = 0.3,
#'       omega_n_driver = 5e-3, omega_p_driver = 2e-3,
#'       ms_driver = -0.5, sigma_driver = 0.5,
#'       alpha_n_driver = 1, beta_n_driver = 10, alpha_p_driver = 1, beta_p_driver = 10)),
#'   clades_wt = list(list(id = "c1", t_clade_wt = 1.5)),
#'   mu = 1e-7, l = 2.7e9, kappa = 20, sigma_count = 0.1, seed = 42
#' )
#' x = init_multirates(sim$tables, m_trunk = sim$truth$m_trunk)
#' x = fit_multirates(x, cmdstan_path = cmdstanr::cmdstan_path(),
#'   method = "sample", chains = 2, ndraws = 200, seed = 45,
#'   mu = 1e-7, l = 2.7e9, t_min = 0, ms_epi = 0, sigma_epi = 0.5,
#'   alpha_lambda = 1, beta_lambda = 1, alpha_n_wt = 1, beta_n_wt = 10,
#'   alpha_p_wt = 1, beta_p_wt = 10, min_kappa = 5, max_kappa = 200,
#'   min_sigma_count = 0.01, max_sigma_count = 1)
#' plot_acf(x)
#' }
#' @export

plot_acf = function(x, params = NULL, lags = 20){

  fit = .multirates_mcmc_fit(x)

  if(is.null(params)) params = .multirates_default_diag_params
  params = .multirates_available_params(fit, params)

  bayesplot::mcmc_acf(fit$draws(variables = params), lags = lags)

}

#' Divergent transitions.
#'
#' Divergences concentrating in a region of parameter space indicate the
#' sampler struggled there and posterior estimates may be biased.
#'
#' @param x PEPI_Multirates object fit with \code{method = "sample"}.
#' @return A ggplot object (see \code{bayesplot::mcmc_nuts_divergence()}).
#' @examples
#' \donttest{
#' sim = simulate_multirates_tree(
#'   sampling_times = c(3, 6, 9), tmrca = 1, t_min = 0,
#'   lambda_n = 1, s_epi = 0.15, omega_n_wt = 5e-3, omega_p_wt = 2e-3,
#'   drivers = list(list(id = "KRAS", t_driver = 2, s_driver = 0.3,
#'       omega_n_driver = 5e-3, omega_p_driver = 2e-3,
#'       ms_driver = -0.5, sigma_driver = 0.5,
#'       alpha_n_driver = 1, beta_n_driver = 10, alpha_p_driver = 1, beta_p_driver = 10)),
#'   clades_wt = list(list(id = "c1", t_clade_wt = 1.5)),
#'   mu = 1e-7, l = 2.7e9, kappa = 20, sigma_count = 0.1, seed = 42
#' )
#' x = init_multirates(sim$tables, m_trunk = sim$truth$m_trunk)
#' x = fit_multirates(x, cmdstan_path = cmdstanr::cmdstan_path(),
#'   method = "sample", chains = 2, ndraws = 200, seed = 45,
#'   mu = 1e-7, l = 2.7e9, t_min = 0, ms_epi = 0, sigma_epi = 0.5,
#'   alpha_lambda = 1, beta_lambda = 1, alpha_n_wt = 1, beta_n_wt = 10,
#'   alpha_p_wt = 1, beta_p_wt = 10, min_kappa = 5, max_kappa = 200,
#'   min_sigma_count = 0.01, max_sigma_count = 1)
#' plot_divergences(x)
#' }
#' @export

plot_divergences = function(x){

  fit = .multirates_mcmc_fit(x)

  bayesplot::mcmc_nuts_divergence(bayesplot::nuts_params(fit), lp = bayesplot::log_posterior(fit))

}

#' Hamiltonian Monte Carlo energy diagnostic.
#'
#' Overlapping marginal and transition energy distributions indicate
#' efficient exploration; a large mismatch signals poor mixing.
#'
#' @param x PEPI_Multirates object fit with \code{method = "sample"}.
#' @return A ggplot object (see \code{bayesplot::mcmc_nuts_energy()}).
#' @examples
#' \donttest{
#' sim = simulate_multirates_tree(
#'   sampling_times = c(3, 6, 9), tmrca = 1, t_min = 0,
#'   lambda_n = 1, s_epi = 0.15, omega_n_wt = 5e-3, omega_p_wt = 2e-3,
#'   drivers = list(list(id = "KRAS", t_driver = 2, s_driver = 0.3,
#'       omega_n_driver = 5e-3, omega_p_driver = 2e-3,
#'       ms_driver = -0.5, sigma_driver = 0.5,
#'       alpha_n_driver = 1, beta_n_driver = 10, alpha_p_driver = 1, beta_p_driver = 10)),
#'   clades_wt = list(list(id = "c1", t_clade_wt = 1.5)),
#'   mu = 1e-7, l = 2.7e9, kappa = 20, sigma_count = 0.1, seed = 42
#' )
#' x = init_multirates(sim$tables, m_trunk = sim$truth$m_trunk)
#' x = fit_multirates(x, cmdstan_path = cmdstanr::cmdstan_path(),
#'   method = "sample", chains = 2, ndraws = 200, seed = 45,
#'   mu = 1e-7, l = 2.7e9, t_min = 0, ms_epi = 0, sigma_epi = 0.5,
#'   alpha_lambda = 1, beta_lambda = 1, alpha_n_wt = 1, beta_n_wt = 10,
#'   alpha_p_wt = 1, beta_p_wt = 10, min_kappa = 5, max_kappa = 200,
#'   min_sigma_count = 0.01, max_sigma_count = 1)
#' plot_energy(x)
#' }
#' @export

plot_energy = function(x){

  fit = .multirates_mcmc_fit(x)

  bayesplot::mcmc_nuts_energy(bayesplot::nuts_params(fit))

}

#' NUTS treedepth diagnostic.
#'
#' Chains repeatedly hitting \code{max_treedepth} indicate the sampler is
#' being cut off before making a U-turn, which hurts efficiency (though
#' not necessarily validity).
#'
#' @param x PEPI_Multirates object fit with \code{method = "sample"}.
#' @return A ggplot object (see \code{bayesplot::mcmc_nuts_treedepth()}).
#' @examples
#' \donttest{
#' sim = simulate_multirates_tree(
#'   sampling_times = c(3, 6, 9), tmrca = 1, t_min = 0,
#'   lambda_n = 1, s_epi = 0.15, omega_n_wt = 5e-3, omega_p_wt = 2e-3,
#'   drivers = list(list(id = "KRAS", t_driver = 2, s_driver = 0.3,
#'       omega_n_driver = 5e-3, omega_p_driver = 2e-3,
#'       ms_driver = -0.5, sigma_driver = 0.5,
#'       alpha_n_driver = 1, beta_n_driver = 10, alpha_p_driver = 1, beta_p_driver = 10)),
#'   clades_wt = list(list(id = "c1", t_clade_wt = 1.5)),
#'   mu = 1e-7, l = 2.7e9, kappa = 20, sigma_count = 0.1, seed = 42
#' )
#' x = init_multirates(sim$tables, m_trunk = sim$truth$m_trunk)
#' x = fit_multirates(x, cmdstan_path = cmdstanr::cmdstan_path(),
#'   method = "sample", chains = 2, ndraws = 200, seed = 45,
#'   mu = 1e-7, l = 2.7e9, t_min = 0, ms_epi = 0, sigma_epi = 0.5,
#'   alpha_lambda = 1, beta_lambda = 1, alpha_n_wt = 1, beta_n_wt = 10,
#'   alpha_p_wt = 1, beta_p_wt = 10, min_kappa = 5, max_kappa = 200,
#'   min_sigma_count = 0.01, max_sigma_count = 1)
#' plot_treedepth(x)
#' }
#' @export

plot_treedepth = function(x){

  fit = .multirates_mcmc_fit(x)

  bayesplot::mcmc_nuts_treedepth(bayesplot::nuts_params(fit), lp = bayesplot::log_posterior(fit))

}

#' Pairs plot with divergences marked.
#'
#' Useful to spot the funnel-shaped regions of parameter space that cause
#' divergences.
#'
#' @param x PEPI_Multirates object fit with \code{method = "sample"}.
#' @param params Vector of Stan parameter base names. Defaults to a
#'   handful of global scalar parameters (pairs plots grow quadratically
#'   with the number of parameters).
#' @return A ggplot/grid object (see \code{bayesplot::mcmc_pairs()}).
#' @examples
#' \donttest{
#' sim = simulate_multirates_tree(
#'   sampling_times = c(3, 6, 9), tmrca = 1, t_min = 0,
#'   lambda_n = 1, s_epi = 0.15, omega_n_wt = 5e-3, omega_p_wt = 2e-3,
#'   drivers = list(list(id = "KRAS", t_driver = 2, s_driver = 0.3,
#'       omega_n_driver = 5e-3, omega_p_driver = 2e-3,
#'       ms_driver = -0.5, sigma_driver = 0.5,
#'       alpha_n_driver = 1, beta_n_driver = 10, alpha_p_driver = 1, beta_p_driver = 10)),
#'   clades_wt = list(list(id = "c1", t_clade_wt = 1.5)),
#'   mu = 1e-7, l = 2.7e9, kappa = 20, sigma_count = 0.1, seed = 42
#' )
#' x = init_multirates(sim$tables, m_trunk = sim$truth$m_trunk)
#' x = fit_multirates(x, cmdstan_path = cmdstanr::cmdstan_path(),
#'   method = "sample", chains = 2, ndraws = 200, seed = 45,
#'   mu = 1e-7, l = 2.7e9, t_min = 0, ms_epi = 0, sigma_epi = 0.5,
#'   alpha_lambda = 1, beta_lambda = 1, alpha_n_wt = 1, beta_n_wt = 10,
#'   alpha_p_wt = 1, beta_p_wt = 10, min_kappa = 5, max_kappa = 200,
#'   min_sigma_count = 0.01, max_sigma_count = 1)
#' plot_pairs(x, params = c("lambda_n","s_epi"))
#' }
#' @export

plot_pairs = function(x, params = NULL){

  fit = .multirates_mcmc_fit(x)

  if(is.null(params)) params = .multirates_default_diag_params
  params = .multirates_available_params(fit, params)

  bayesplot::mcmc_pairs(fit$draws(variables = params), np = bayesplot::nuts_params(fit))

}

#' One-call convergence dashboard.
#'
#' Arranges Rhat, effective sample size, divergences and the energy
#' diagnostic into a single 2x2 figure - a quick first look at whether an
#' MCMC fit is trustworthy.
#'
#' @param x PEPI_Multirates object fit with \code{method = "sample"}.
#' @param params Vector of Stan parameter base names for the Rhat/ESS
#'   panels. Defaults to every sampled parameter.
#' @return A combined plot (see \code{ggpubr::ggarrange()}).
#' @examples
#' \donttest{
#' sim = simulate_multirates_tree(
#'   sampling_times = c(3, 6, 9), tmrca = 1, t_min = 0,
#'   lambda_n = 1, s_epi = 0.15, omega_n_wt = 5e-3, omega_p_wt = 2e-3,
#'   drivers = list(list(id = "KRAS", t_driver = 2, s_driver = 0.3,
#'       omega_n_driver = 5e-3, omega_p_driver = 2e-3,
#'       ms_driver = -0.5, sigma_driver = 0.5,
#'       alpha_n_driver = 1, beta_n_driver = 10, alpha_p_driver = 1, beta_p_driver = 10)),
#'   clades_wt = list(list(id = "c1", t_clade_wt = 1.5)),
#'   mu = 1e-7, l = 2.7e9, kappa = 20, sigma_count = 0.1, seed = 42
#' )
#' x = init_multirates(sim$tables, m_trunk = sim$truth$m_trunk)
#' x = fit_multirates(x, cmdstan_path = cmdstanr::cmdstan_path(),
#'   method = "sample", chains = 2, ndraws = 200, seed = 45,
#'   mu = 1e-7, l = 2.7e9, t_min = 0, ms_epi = 0, sigma_epi = 0.5,
#'   alpha_lambda = 1, beta_lambda = 1, alpha_n_wt = 1, beta_n_wt = 10,
#'   alpha_p_wt = 1, beta_p_wt = 10, min_kappa = 5, max_kappa = 200,
#'   min_sigma_count = 0.01, max_sigma_count = 1)
#' plot_diagnostics(x)
#' }
#' @export

plot_diagnostics = function(x, params = NULL){

  p1 = plot_rhat(x, params = params) + ggplot2::labs(title = "Rhat")
  p2 = plot_ess(x, params = params) + ggplot2::labs(title = "Effective sample size ratio")
  p3 = plot_divergences(x) + ggplot2::labs(title = "Divergences")
  p4 = plot_energy(x) + ggplot2::labs(title = "Energy (E-BFMI)")

  ggpubr::ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)

}
