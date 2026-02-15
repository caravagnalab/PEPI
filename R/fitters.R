#' Fit a PEPI model using Stan and update the PEPI object
#'
#' This function takes a PEPI object, fits a tumor evolutionary model
#' using Stan, and updates the PEPI object with the fitted model and Stan data.
#'
#' Fit the PEPI Bayesian model
#'
#' This function fits the PEPI Bayesian model to bulk sequencing and
#' longitudinal population counts, jointly inferring clonal structure,
#' driver expansions, epigenetic switching, and fitness advantages of
#' the positive epistate (⊕).
#'
#' @param pepi A PEPI object containing clade statistics, counts, and genomic constants.
#'
#' @param model_type Character. Functional form of the ω_⊕ prior.
#'   Options: "logistic" (default), "gumbel", or "log".
#'
#' @param ms_driver_n Numeric vector. Mean (μ) of the lognormal prior for
#'   drivers of type \code{driver_n}. Default automatically set if \code{NULL}.
#'
#' @param sigma_driver_n Numeric vector. Standard deviation (σ) of the
#'   lognormal prior for drivers of type \code{driver_n}. Default automatically set if \code{NULL}.
#'
#' @param ms_dc Numeric vector. Mean of the lognormal prior for dc drivers.
#'   Default automatically set if \code{NULL}.
#'
#' @param sigma_dc Numeric vector. Standard deviation of the lognormal prior for dc drivers.
#'   Default automatically set if \code{NULL}.
#'
#' @param ms_cd Numeric vector. Mean of the lognormal prior for cd drivers.
#'   Default automatically set if \code{NULL}.
#'
#' @param sigma_cd Numeric vector. Standard deviation of the lognormal prior for cd drivers.
#'   Default automatically set if \code{NULL}.
#'
#' @param alpha_lambda Numeric. Shape parameter of the Gamma prior on the baseline
#'   population growth rate λ.
#'
#' @param beta_lambda Numeric. Rate parameter of the Gamma prior on the baseline
#'   population growth rate λ.
#'
#' @param alpha_plus Numeric. Shape parameter of the Gamma prior on the variance
#'   of Gaussian observation noise for counts in the ⊕ epistate.
#'
#' @param beta_plus Numeric. Rate parameter of the Gamma prior on the variance
#'   of Gaussian observation noise for counts in the ⊕ epistate.
#'
#' @param alpha_minus Numeric. Shape parameter of the Gamma prior on the variance
#'   of Gaussian observation noise for counts in the ⊖ epistate.
#'
#' @param beta_minus Numeric. Rate parameter of the Gamma prior on the variance
#'   of Gaussian observation noise for counts in the ⊖ epistate.
#'
#' @param alpha_n Numeric. Shape parameter of the Gamma prior on the epigenetic
#'   switching rate ω_n (⊖ → ⊕).
#'
#' @param beta_n Numeric. Rate parameter of the Gamma prior on ω_n.
#'
#' @param alpha_p Numeric. Shape parameter of the Gamma prior on the epigenetic
#'   switching rate ω_p (⊕ → ⊖).
#'
#' @param beta_p Numeric. Rate parameter of the Gamma prior on ω_p.
#'
#' @param t_min Numeric. Lower bound for the prior on the time to the most recent
#'   common ancestor (MRCA) of clades and driver lineages, relative to the most recent
#'   sampling time.
#'
#' @param ms_epi Numeric. Mean (μ) of the lognormal prior on the fitness
#'   advantage of the positive epistate (⊕).
#'
#' @param sigma_epi Numeric. Standard deviation (σ) of the lognormal prior
#'   on the fitness advantage of the positive epistate (⊕).
#'
#' @param ccf_thr_clade Numeric. Minimum cancer cell fraction (CCF) threshold for
#'   clades to be included in the inference.
#'
#' @param ccf_thr_count Numeric. Minimum CCF threshold for population counts to be
#'   included in the inference.
#'
#' @param include_bp Logical. If TRUE, includes branching-process–based terms for the
#'   total positive (⊕) counts at the first time point in the likelihood.
#'   Default FALSE.
#'
#' @param n_chains Integer. Number of Markov chains to run in Stan. Default 4.
#'
#' @param adapt_delta Numeric. Target acceptance probability for the Stan NUTS sampler.
#'   Default 0.8.
#'
#' @param iter_warmup Integer. Number of warmup (burn-in) iterations per chain.
#'   Default 1000.
#'
#' @param iter_sampling Integer. Number of sampling iterations per chain.
#'   Default 1000.
#'
#' @param seed Integer. Random seed for reproducibility. Default 1.
#'
#' @param parallel_chains Integer. Number of chains to run in parallel. Default 1.
#'
#' @return A PEPI object updated with:
#' \item{stan_data}{The Stan data list used for fitting.}
#' \item{fit}{The fitted \code{cmdstanr::CmdStanMCMC} object.}
#' 
#' @export



fit_pepi <- function(
    pepi,
    model_type = "logistic",
    ms_driver_n = NULL,
    sigma_driver_n = NULL,
    ms_dc = NULL,
    sigma_dc = NULL,
    ms_cd = NULL,
    sigma_cd = NULL,
    alpha_lambda = 1, beta_lambda = 1,
    alpha_plus = 10, beta_plus = 250,
    alpha_minus = 10, beta_minus = 250,
    alpha_n = 0.5, beta_n = 20,
    alpha_p = 0.2, beta_p = 10,
    t_min = 0,
    ms_epi = 0,
    sigma_epi = 0.5,
    ccf_thr_clade = 0.05,
    ccf_thr_count = 0.02,
    include_bp = FALSE,
    n_chains = 4,
    parallel_chains = 4,
    adapt_delta = 0.8,
    iter_warmup = 1000,
    iter_sampling = 1000,
    seed = 1
) {
  
  # --- Detect cmdstan path and Stan model file ---
  if (!requireNamespace("cmdstanr", quietly = TRUE)) stop("Package 'cmdstanr' required.")
  cmdstanr::set_cmdstan_path()
  
  models_path <- file.path(getwd(), "models")
  if (!model_type %in% c("logistic","gumbel","log")) stop("model_type must be 'logistic', 'gumbel', or 'log'.")
  stan_file <- file.path(models_path, model_type, "model.stan")
  if (!file.exists(stan_file)) stop("Stan model file not found: ", stan_file)
  
  # --- Extract Stan data from pepi ---
  stan_data <- get_stan_data_pepi(pepi)
  
  # --- Fill prior parameters ---
  N_driver_n <- stan_data$N_driver_n
  N_dc <- stan_data$N_dc
  N_cd <- stan_data$N_cd
  
  if (is.null(ms_driver_n)) ms_driver_n <- if (N_driver_n>0) rep(0,N_driver_n) else numeric(0)
  if (is.null(sigma_driver_n)) sigma_driver_n <- if (N_driver_n>0) rep(1.5,N_driver_n) else numeric(0)
  
  if (is.null(ms_dc)) ms_dc <- if (N_dc>0) rep(0,N_dc) else numeric(0)
  if (is.null(sigma_dc)) sigma_dc <- if (N_dc>0) rep(1.5,N_dc) else numeric(0)
  
  if (is.null(ms_cd)) ms_cd <- if (N_cd>0) rep(0,N_cd) else numeric(0)
  if (is.null(sigma_cd)) sigma_cd <- if (N_cd>0) rep(1.5,N_cd) else numeric(0)
  
  # --- Combine data and priors for Stan ---
  stan_data <- c(
    stan_data,
    list(
      ms_driver_n = ms_driver_n,
      sigma_driver_n = sigma_driver_n,
      ms_dc = ms_dc,
      sigma_dc = sigma_dc,
      ms_cd = ms_cd,
      sigma_cd = sigma_cd,
      alpha_lambda = alpha_lambda,
      beta_lambda = beta_lambda,
      alpha_plus = alpha_plus,
      beta_plus = beta_plus,
      alpha_minus = alpha_minus,
      beta_minus = beta_minus,
      t_min = t_min,
      ms_epi = ms_epi,
      sigma_epi = sigma_epi,
      ccf_thr_count = ccf_thr_count,
      ccf_thr_clade = ccf_thr_clade,
      alpha_n = alpha_n,
      beta_n = beta_n,
      alpha_p = alpha_p,
      beta_p = beta_p,
      include_bp = include_bp
    )
  )
  
  # --- Compile and fit Stan model ---
  mod <- cmdstanr::cmdstan_model(stan_file)
  fit <- mod$sample(
    data = stan_data,
    chains = n_chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta,
    seed = seed
  )
  
  # --- Update PEPI object ---
  pepi$stan_data <- stan_data
  pepi$fit <- fit
  
  return(pepi)
}
