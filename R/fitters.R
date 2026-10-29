#' Build a default set of initial values for the multirates model.
#'
#' Several time parameters have bounds that depend on the *sampled* value of
#' `tmrca` (e.g. `t_clade_wt` is bounded above by the earliest sampling time,
#' below by `tmrca`), so Stan's default random initialization frequently
#' draws an infeasible `tmrca` and fails before sampling starts. This builds
#' a single feasible starting point instead.
#'
#' @param data Stan data list, as built by \code{build_stan_data_multirates()}.
#' @return A list (of length 1) suitable for cmdstanr's `init` argument.
#' @keywords internal

.multirates_default_init = function(data){

  t_min = data$t_min
  tmax = data$times[1]
  earliest_sample = data$times[data$sampling_index[1]]

  tmrca0 = t_min + 0.3 * (earliest_sample - t_min)

  init = list(tmrca = tmrca0)

  if(data$N_clades_wt > 0){
    init$t_clade_wt = rep(tmrca0 + 0.5 * (earliest_sample - tmrca0), data$N_clades_wt)
  }
  if(data$N_driver > 0){
    init$t_driver = rep(tmrca0 + 0.5 * (tmax - tmrca0), data$N_driver)
  }
  if(data$N_driver_n > 0){
    init$t_driver_n = rep(tmrca0 + 0.5 * (tmax - tmrca0), data$N_driver_n)
  }
  if(data$N_driver_p > 0){
    init$t_driver_p = rep(tmrca0 + 0.5 * (tmax - tmrca0), data$N_driver_p)
  }

  list(init)

}


#' Fit the unified multirates model.
#'
#' Compiles and runs \code{inst/multirates_positive_s.stan} on a
#' \code{PEPI_Multirates} object built with \code{init_multirates()}.
#'
#' @param x PEPI_Multirates object.
#' @param cmdstan_path String specifying the path to the cmdstan installation.
#' @param seed Seed of the computation.
#' @param ndraws Number of posterior draws.
#' @param chains Number of MCMC chains.
#' @param parallel_chains Number of parallel MCMC chains.
#' @param adapt_delta Target average acceptance probability.
#' @param max_treedepth Maximum number of trajectory steps the NUTS sampler can take per iteration.
#' @param init List of initialization values, or NULL to use a built-in
#'   default (see \code{.multirates_default_init()}).
#' @param mu Mutation rate per division per bp per allele.
#' @param l Length of the genome.
#' @param t_min Lower bound for the tmrca prior.
#' @param ms_epi,sigma_epi Lognormal prior hyperparameters for s_epi.
#' @param alpha_lambda,beta_lambda Gamma prior hyperparameters for lambda_n.
#' @param alpha_n_wt,beta_n_wt,alpha_p_wt,beta_p_wt Gamma prior hyperparameters
#'   for the wt switch rates omega_n_wt/omega_p_wt.
#' @param min_kappa,max_kappa Uniform prior bounds for kappa.
#' @param min_sigma_count,max_sigma_count Uniform prior bounds for sigma_count.
#' @param include_poisson Default (0/1) of the Poisson terms on the driver
#'   mutations, used by the per-type switches below when they are not given.
#' @param include_poisson_driver,include_poisson_driver_n,include_poisson_driver_p
#'   Switches (0/1) of the Poisson term on the mutations of the driver,
#'   driver_n and driver_p populations (for driver_p, an upper bound).
#' @param include_poisson_dc Switch (0/1) of the Poisson term on the driver
#'   mutations of the dc populations; the clade terms are always included.
#' @param include_trunk Switch (0/1) of the trunk regularization of tmrca.
#' @return PEPI_Multirates object with `inference$multirates`, `stan_data$multirates`
#'   and `index_maps` populated.
#' @examples
#' \donttest{
#' source(system.file("examples", "tables_example.R", package = "PEPI"))
#' x = fit_multirates(x, seed = 31, chains = 2, parallel_chains = 2, ndraws = 200,
#'   init = rep(list(list(omega_p_wt = 3e-3, omega_p_driver = 3e-3)), 2),
#'   mu = mu$mean, l = 51e06, t_min = -20, ms_epi = -0.8, sigma_epi = 0.9,
#'   alpha_lambda = 10, beta_lambda = 10, alpha_n_wt = 1, beta_n_wt = 1e4,
#'   alpha_p_wt = 0.2, beta_p_wt = 10, min_kappa = 3000, max_kappa = 5000,
#'   min_sigma_count = 0.05, max_sigma_count = 0.15, include_poisson = 1, include_trunk = 0)
#' }
#' @export

fit_multirates = function(x, cmdstan_path = cmdstanr::cmdstan_path(),
                          seed = 45, ndraws = 1000, chains = 4, parallel_chains = 4,  
                          adapt_delta = 0.99,max_treedepth = 15,
                          init = NULL,
                          mu = 1e-7, l = 2.7e9, t_min = 0,
                          ms_epi = 0, sigma_epi = 0.5,
                          alpha_lambda = 1, beta_lambda = 1,
                          alpha_n_wt = 1, beta_n_wt = 10, alpha_p_wt = 1, beta_p_wt = 10,
                          min_kappa = 10, max_kappa = 1000,
                          min_sigma_count = 0.01, max_sigma_count = 1,
                          include_poisson = 1L,
    include_poisson_driver = include_poisson, include_poisson_driver_n = include_poisson,
    include_poisson_driver_p = include_poisson, include_poisson_dc = 1L,
    include_trunk = 1L){

  cmdstanr::set_cmdstan_path(cmdstan_path)

  if(is.null(x$multirates)){
    stop("no multirates input tables")
  }

  built = build_stan_data_multirates(x, mu = mu, l = l, t_min = t_min,
            ms_epi = ms_epi, sigma_epi = sigma_epi,
            alpha_lambda = alpha_lambda, beta_lambda = beta_lambda,
            alpha_n_wt = alpha_n_wt, beta_n_wt = beta_n_wt,
            alpha_p_wt = alpha_p_wt, beta_p_wt = beta_p_wt,
            min_kappa = min_kappa, max_kappa = max_kappa,
            min_sigma_count = min_sigma_count, max_sigma_count = max_sigma_count,
            include_poisson = include_poisson,
            include_poisson_driver = include_poisson_driver,
            include_poisson_driver_n = include_poisson_driver_n,
            include_poisson_driver_p = include_poisson_driver_p,
            include_poisson_dc = include_poisson_dc,
            include_trunk = include_trunk)

  # compile into a user cache directory, never next to the .stan file: an
  # executable shipped in the package (or built for another platform) would
  # otherwise be considered up to date and used as is
  exe_dir = tools::R_user_dir("PEPI", which = "cache")
  dir.create(exe_dir, recursive = TRUE, showWarnings = FALSE)

  mod = cmdstanr::cmdstan_model(system.file("multirates_positive_s.stan", package = "PEPI"), dir = exe_dir)

  if(is.null(init)) init = rep(.multirates_default_init(built$data), chains)

  # fit = mod$sample(data = built$data, seed = seed, init = init,
  #             chains = chains,  parallel_chains = parallel_chains,
  #             iter_sampling = ndraws,
  #             adapt_delta = adapt_delta, max_treedepth = max_treedepth)
  # cmdstanr >= 0.9 reads the fit metadata from the CmdStan config files,
  # which are written only with save_cmdstan_config = TRUE
  fit = mod$sample(data = built$data, seed = seed, init = init,
              chains = chains,  parallel_chains = parallel_chains,
              iter_sampling = ndraws,
              adapt_delta = adapt_delta, max_treedepth = max_treedepth,
              save_cmdstan_config = TRUE)


  x$inference$multirates = fit
  x$stan_data$multirates = built$data
  x$index_maps = built$index_maps

  return(x)

}
