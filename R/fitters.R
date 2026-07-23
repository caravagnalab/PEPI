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
#' @param method "variational" (default, fast approximate posterior) or "sample" (full MCMC).
#' @param ndraws Number of posterior draws (variational: output_samples; sample: iter_sampling).
#' @param chains Number of MCMC chains (only used when method = "sample").
#' @param seed Seed of the computation.
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
#' @param include_poisson,include_trunk Stan model switches (0/1).
#' @return PEPI_Multirates object with `inference$multirates`, `stan_data$multirates`
#'   and `index_maps` populated.
#' @examples
#' \dontrun{
#' fit_multirates(x, cmdstan_path = cmdstanr::cmdstan_path(),
#'   method = "variational", ndraws = 1000, seed = 45,
#'   mu = 1e-7, l = 2.7e9, t_min = 0, ms_epi = 0, sigma_epi = 0.5,
#'   alpha_lambda = 1, beta_lambda = 1, alpha_n_wt = 1, beta_n_wt = 10,
#'   alpha_p_wt = 1, beta_p_wt = 10, min_kappa = 10, max_kappa = 1000,
#'   min_sigma_count = 0.01, max_sigma_count = 1)
#' }
#' @export

fit_multirates = function(x, cmdstan_path = cmdstanr::cmdstan_path(),
                          method = c("variational","sample"),
                          ndraws = 1000, chains = 4, seed = 45, init = NULL,
                          mu = 1e-7, l = 2.7e9, t_min = 0,
                          ms_epi = 0, sigma_epi = 0.5,
                          alpha_lambda = 1, beta_lambda = 1,
                          alpha_n_wt = 1, beta_n_wt = 10, alpha_p_wt = 1, beta_p_wt = 10,
                          min_kappa = 10, max_kappa = 1000,
                          min_sigma_count = 0.01, max_sigma_count = 1,
                          include_poisson = 1L, include_trunk = 1L){

  method = match.arg(method)

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
            include_poisson = include_poisson, include_trunk = include_trunk)

  if(is.null(init)){
    init = .multirates_default_init(built$data)
  }

  mod = cmdstanr::cmdstan_model(system.file("multirates_positive_s.stan", package = "PEPI"))

  fit = if(method == "variational"){

    mod$variational(data = built$data, seed = seed, init = init,
                    output_samples = ndraws, algorithm = "fullrank")

  }else{

    mod$sample(data = built$data, seed = seed, init = init,
              chains = chains, iter_sampling = ndraws)

  }

  x$inference$multirates = fit
  x$stan_data$multirates = built$data
  x$index_maps = built$index_maps

  return(x)

}
