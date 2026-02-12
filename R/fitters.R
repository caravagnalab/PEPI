#' Fit a PEPI model using Stan and update the PEPI object
#'
#' This function takes a PEPI object, fits a tumor evolutionary model
#' using Stan, and updates the PEPI object with the fitted model and Stan data.
#'
#' @param pepi A PEPI object containing clade_statistics, Counts, and genomic_constants.
#' @param model_type Character. The model prior to use: "logistic" (default), "gumbel", or "log".
#' @param ms_driver_n Numeric vector. Prior mean for driver_n clusters. Default automatically set if NULL.
#' @param sigma_driver_n Numeric vector. Prior SD for driver_n clusters. Default automatically set if NULL.
#' @param ms_dc Numeric vector. Prior mean for DC clusters. Default automatically set if NULL.
#' @param sigma_dc Numeric vector. Prior SD for DC clusters. Default automatically set if NULL.
#' @param ms_cd Numeric vector. Prior mean for CD clusters. Default automatically set if NULL.
#' @param sigma_cd Numeric vector. Prior SD for CD clusters. Default automatically set if NULL.
#' @param alpha_lambda Numeric. Hyperparameter for lambda prior. Default 1.
#' @param beta_lambda Numeric. Hyperparameter for lambda prior. Default 1.
#' @param alpha_plus Numeric. Hyperparameter for positive epistate prior. Default 10.
#' @param beta_plus Numeric. Hyperparameter for positive epistate prior. Default 250.
#' @param alpha_minus Numeric. Hyperparameter for negative epistate prior. Default 1.
#' @param beta_minus Numeric. Hyperparameter for negative epistate prior. Default 1.
#' @param alpha_n Numeric. Hyperparameter for negative counts prior. Default 0.5.
#' @param beta_n Numeric. Hyperparameter for negative counts prior. Default 20.
#' @param alpha_p Numeric. Hyperparameter for positive counts prior. Default 0.2.
#' @param beta_p Numeric. Hyperparameter for positive counts prior. Default 10.
#' @param t_min Numeric. Minimum time prior. Default 0.
#' @param ms_epi Numeric. Prior mean for epistate. Default 0.
#' @param sigma_epi Numeric. Prior SD for epistate. Default 0.5.
#' @param ccf_thr_clade Numeric. Threshold for clade CCF. Default 0.05.
#' @param ccf_thr_count Numeric. Threshold for count CCF. Default 0.02.
#' @param include_bp Logical. Whether to include base pair information. Default FALSE.
#' @param n_chains Integer. Number of Stan chains. Default 4.
#' @param adapt_delta Numeric. Stan adapt_delta parameter. Default 0.8.
#' @param iter_warmup Integer. Number of warmup iterations. Default 1000.
#' @param iter_sampling Integer. Number of sampling iterations. Default 1000.
#'
#' @return A PEPI object updated with:
#' \item{stan_data}{The Stan data list used for fitting.}
#' \item{fit}{The fitted \code{cmdstanr::CmdStanMCMC} object.}
#'
#' @examples
#' \dontrun{
#' pepi_object <- fit_pepi(pepi_object, model_type = "logistic")
#' }
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
    alpha_minus = 1, beta_minus = 1,
    alpha_n = 0.5, beta_n = 20,
    alpha_p = 0.2, beta_p = 10,
    t_min = 0,
    ms_epi = 0,
    sigma_epi = 0.5,
    ccf_thr_clade = 0.05,
    ccf_thr_count = 0.02,
    include_bp = FALSE,
    n_chains = 4,
    adapt_delta = 0.8,
    iter_warmup = 1000,
    iter_sampling = 1000
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
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta
  )
  
  # --- Update PEPI object ---
  pepi$stan_data <- stan_data
  pepi$fit <- fit
  
  return(pepi)
}
