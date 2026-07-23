#' Stan data construction for the unified multirates model

#' Pivot a long-format CCF table into an array shaped for Stan.
#'
#' @param long_df Long-format tibble with an id column, `time`, `ccf`
#'   (and `epistate` if `epistate = TRUE`).
#' @param id_col Name of the id column in `long_df`.
#' @param id_levels Character vector of ids, in the order they should appear
#'   along the first array dimension.
#' @param time_levels Numeric vector of sampling times, in ascending
#'   chronological order (first = earliest sample).
#' @param epistate If TRUE, returns a `[length(id_levels), length(time_levels), 2]`
#'   array (epistate `-`/`+` on the third dimension); otherwise a
#'   `[length(id_levels), length(time_levels)]` matrix.
#' @return An array/matrix of CCF values.
#' @keywords internal

pivot_ccf = function(long_df, id_col, id_levels, time_levels, epistate = TRUE){

  n_id = length(id_levels)
  n_t = length(time_levels)

  if(epistate){

    arr = array(NA_real_, dim = c(n_id, n_t, 2))
    if(n_id == 0) return(arr)

    for(i in seq_along(id_levels)){
      for(t in seq_along(time_levels)){
        for(e in 1:2){

          es = c("-","+")[e]
          v = long_df$ccf[long_df[[id_col]] == id_levels[i] & long_df$time == time_levels[t] & long_df$epistate == es]

          if(length(v) != 1){
            stop(sprintf("expected exactly 1 ccf for %s=%s, time=%s, epistate=%s (found %d)",
                          id_col, id_levels[i], time_levels[t], es, length(v)))
          }

          arr[i,t,e] = v
        }
      }
    }

  }else{

    arr = matrix(NA_real_, nrow = n_id, ncol = n_t)
    if(n_id == 0) return(arr)

    for(i in seq_along(id_levels)){
      for(t in seq_along(time_levels)){

        v = long_df$ccf[long_df[[id_col]] == id_levels[i] & long_df$time == time_levels[t]]

        if(length(v) != 1){
          stop(sprintf("expected exactly 1 ccf for %s=%s, time=%s (found %d)",
                        id_col, id_levels[i], time_levels[t], length(v)))
        }

        arr[i,t] = v
      }
    }
  }

  arr
}


#' Build the Stan data list for the unified multirates model.
#'
#' Turns the tables stored on a `PEPI_Multirates` object (see
#' \code{init_multirates()}) into the exact data list required by
#' \code{inst/multirates_positive_s.stan}, plus an `index_maps` list that
#' lets getters translate Stan array indices back into ids/times.
#'
#' @param x PEPI_Multirates object.
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
#' @return A list with elements `data` (the Stan data list) and `index_maps`.
#' @examples
#' sim = simulate_multirates_tree(
#'   sampling_times = c(3, 6, 9), tmrca = 1, t_min = 0,
#'   lambda_n = 1, s_epi = 0.15, omega_n_wt = 5e-3, omega_p_wt = 2e-3,
#'   clades_wt = list(list(id = "c1", t_clade_wt = 1.5)),
#'   mu = 1e-7, l = 2.7e9, kappa = 20, sigma_count = 0.1, seed = 42
#' )
#' x = init_multirates(sim$tables, m_trunk = sim$truth$m_trunk)
#' build_stan_data_multirates(x, mu = 1e-7, l = 2.7e9, t_min = 0,
#'   ms_epi = 0, sigma_epi = 0.5, alpha_lambda = 1, beta_lambda = 1,
#'   alpha_n_wt = 1, beta_n_wt = 10, alpha_p_wt = 1, beta_p_wt = 10,
#'   min_kappa = 5, max_kappa = 200, min_sigma_count = 0.01, max_sigma_count = 1)
#' @export

build_stan_data_multirates = function(x, mu, l, t_min, ms_epi, sigma_epi,
    alpha_lambda, beta_lambda, alpha_n_wt, beta_n_wt, alpha_p_wt, beta_p_wt,
    min_kappa, max_kappa, min_sigma_count, max_sigma_count,
    include_poisson = 1L, include_trunk = 1L){

  tables = x$multirates

  if(is.null(tables)){
    stop("no multirates input tables")
  }

  check_input_multirates(tables)

  # ---- time grid ------------------------------------------------------

  sampling_rows = tables$population_sizes %>% dplyr::filter(type == "sampling") %>% dplyr::arrange(time)
  intermediate_rows = tables$population_sizes %>% dplyr::filter(type == "intermediate") %>% dplyr::arrange(time)

  n_sampling_times = nrow(sampling_rows)
  n_intermediate_times = nrow(intermediate_rows)
  n_times = n_sampling_times + n_intermediate_times

  times = numeric(n_times)
  times[1] = sampling_rows$time[n_sampling_times]
  if(n_sampling_times > 1){
    times[2:n_sampling_times] = sampling_rows$time[1:(n_sampling_times - 1)]
  }

  sampling_index = integer(n_sampling_times)
  sampling_index[n_sampling_times] = 1L
  if(n_sampling_times > 1){
    sampling_index[1:(n_sampling_times - 1)] = 2:n_sampling_times
  }

  if(n_intermediate_times > 0){
    times[(n_sampling_times + 1):n_times] = intermediate_rows$time
    intermediate_index = as.integer((n_sampling_times + 1):n_times)
  }else{
    intermediate_index = integer(0)
  }

  sampling_times = sampling_rows$time

  # ---- wt ---------------------------------------------------------------

  ccf_wt = matrix(NA_real_, nrow = n_sampling_times, ncol = 2)
  for(t in seq_len(n_sampling_times)){
    for(e in 1:2){
      es = c("-","+")[e]
      v = tables$wt_ccf %>% dplyr::filter(time == sampling_times[t], epistate == es) %>% dplyr::pull(ccf)
      if(length(v) != 1) stop(sprintf("expected exactly 1 wt_ccf for time=%s, epistate=%s", sampling_times[t], es))
      ccf_wt[t,e] = v
    }
  }

  # ---- driver -------------------------------------------------------------

  driver_muts = tables$driver_muts
  N_driver = if(is.null(driver_muts)) 0L else nrow(driver_muts)

  if(N_driver > 0){
    m_driver = as.integer(driver_muts$m_driver)
    ms_driver = as.numeric(driver_muts$ms_driver)
    sigma_driver = as.numeric(driver_muts$sigma_driver)
    alpha_n_driver = as.numeric(driver_muts$alpha_n_driver)
    beta_n_driver = as.numeric(driver_muts$beta_n_driver)
    alpha_p_driver = as.numeric(driver_muts$alpha_p_driver)
    beta_p_driver = as.numeric(driver_muts$beta_p_driver)
    ccf_driver = pivot_ccf(tables$driver_ccf, "driver_id", driver_muts$driver_id, sampling_times, epistate = TRUE)
  }else{
    m_driver = integer(0); ms_driver = numeric(0); sigma_driver = numeric(0)
    alpha_n_driver = numeric(0); beta_n_driver = numeric(0)
    alpha_p_driver = numeric(0); beta_p_driver = numeric(0)
    ccf_driver = array(numeric(0), dim = c(0, n_sampling_times, 2))
  }

  # ---- driver_n -------------------------------------------------------------

  driver_n_muts = tables$driver_n_muts
  N_driver_n = if(is.null(driver_n_muts)) 0L else nrow(driver_n_muts)

  if(N_driver_n > 0){
    m_driver_n = as.integer(driver_n_muts$m_driver_n)
    ms_driver_n = as.numeric(driver_n_muts$ms_driver_n)
    sigma_driver_n = as.numeric(driver_n_muts$sigma_driver_n)
    ccf_driver_n = pivot_ccf(tables$driver_n_ccf, "driver_n_id", driver_n_muts$driver_n_id, sampling_times, epistate = FALSE)
  }else{
    m_driver_n = integer(0); ms_driver_n = numeric(0); sigma_driver_n = numeric(0)
    ccf_driver_n = matrix(numeric(0), nrow = 0, ncol = n_sampling_times)
  }

  # ---- driver_p -------------------------------------------------------------
  # note: the Stan data{} block has no m_driver_p array (asymmetric with driver_n by design)

  driver_p_muts = tables$driver_p_muts
  N_driver_p = if(is.null(driver_p_muts)) 0L else nrow(driver_p_muts)

  if(N_driver_p > 0){
    ms_driver_p = as.numeric(driver_p_muts$ms_driver_p)
    sigma_driver_p = as.numeric(driver_p_muts$sigma_driver_p)
    ccf_driver_p = pivot_ccf(tables$driver_p_ccf, "driver_p_id", driver_p_muts$driver_p_id, sampling_times, epistate = FALSE)
  }else{
    ms_driver_p = numeric(0); sigma_driver_p = numeric(0)
    ccf_driver_p = matrix(numeric(0), nrow = 0, ncol = n_sampling_times)
  }

  # ---- dc (driver + clade combos) ------------------------------------------

  dc_muts = tables$dc_muts
  N_dc = if(is.null(dc_muts)) 0L else nrow(dc_muts)

  if(N_dc > 0){
    m_dc = cbind(as.integer(dc_muts$m_dc_driver), as.integer(dc_muts$m_dc_clade))
    ms_dc = as.numeric(dc_muts$ms_dc)
    sigma_dc = as.numeric(dc_muts$sigma_dc)
    alpha_n_dc = as.numeric(dc_muts$alpha_n_dc)
    beta_n_dc = as.numeric(dc_muts$beta_n_dc)
    alpha_p_dc = as.numeric(dc_muts$alpha_p_dc)
    beta_p_dc = as.numeric(dc_muts$beta_p_dc)

    dc_ccf_driver_long = tables$dc_ccf %>% dplyr::filter(quantity == "driver")
    dc_ccf_clade_long = tables$dc_ccf %>% dplyr::filter(quantity == "clade")

    ccf_dc_driver = pivot_ccf(dc_ccf_driver_long, "dc_id", dc_muts$dc_id, sampling_times, epistate = TRUE)
    ccf_dc_clade = pivot_ccf(dc_ccf_clade_long, "dc_id", dc_muts$dc_id, sampling_times, epistate = TRUE)
  }else{
    m_dc = matrix(integer(0), nrow = 0, ncol = 2)
    ms_dc = numeric(0); sigma_dc = numeric(0)
    alpha_n_dc = numeric(0); beta_n_dc = numeric(0)
    alpha_p_dc = numeric(0); beta_p_dc = numeric(0)
    ccf_dc_driver = array(numeric(0), dim = c(0, n_sampling_times, 2))
    ccf_dc_clade = array(numeric(0), dim = c(0, n_sampling_times, 2))
  }

  # ---- clades_wt -------------------------------------------------------------

  clades_wt = tables$clades_wt
  N_clades_wt = if(is.null(clades_wt)) 0L else nrow(clades_wt)

  if(N_clades_wt > 0){
    m_clade_wt = as.integer(clades_wt$m_clade_wt)
    ccf_clade_wt = pivot_ccf(tables$clades_wt_ccf, "clade_id", clades_wt$clade_id, sampling_times, epistate = TRUE)
  }else{
    m_clade_wt = integer(0)
    ccf_clade_wt = array(numeric(0), dim = c(0, n_sampling_times, 2))
  }

  # ---- assemble -------------------------------------------------------------

  data = list(
    N_clades_wt = N_clades_wt,
    n_times = n_times,
    n_sampling_times = n_sampling_times,
    n_intermediate_times = n_intermediate_times,
    sampling_index = sampling_index,
    intermediate_index = intermediate_index,
    N_driver = N_driver,
    N_dc = N_dc,
    N_driver_n = N_driver_n,
    N_driver_p = N_driver_p,
    include_poisson = as.integer(include_poisson),
    m_trunk = as.integer(x$m_trunk),
    include_trunk = as.integer(include_trunk),

    m_driver = m_driver,
    ccf_driver = ccf_driver,
    ms_driver = ms_driver,
    sigma_driver = sigma_driver,

    m_driver_n = m_driver_n,
    ccf_driver_n = ccf_driver_n,
    ms_driver_n = ms_driver_n,
    sigma_driver_n = sigma_driver_n,

    ccf_driver_p = ccf_driver_p,
    ms_driver_p = ms_driver_p,
    sigma_driver_p = sigma_driver_p,

    m_dc = m_dc,
    ccf_dc_driver = ccf_dc_driver,
    ccf_dc_clade = ccf_dc_clade,
    ms_dc = ms_dc,
    sigma_dc = sigma_dc,

    times = times,
    zminus = sampling_rows$zminus,
    zplus = sampling_rows$zplus,
    ztot = if(n_intermediate_times > 0) intermediate_rows$ztot else numeric(0),
    m_clade_wt = m_clade_wt,
    ccf_clade_wt = ccf_clade_wt,
    ccf_wt = ccf_wt,

    t_min = t_min,
    ms_epi = ms_epi,
    sigma_epi = sigma_epi,
    alpha_lambda = alpha_lambda,
    beta_lambda = beta_lambda,
    alpha_n_wt = alpha_n_wt,
    beta_n_wt = beta_n_wt,
    alpha_p_wt = alpha_p_wt,
    beta_p_wt = beta_p_wt,

    alpha_n_driver = alpha_n_driver,
    beta_n_driver = beta_n_driver,
    alpha_p_driver = alpha_p_driver,
    beta_p_driver = beta_p_driver,

    alpha_n_dc = alpha_n_dc,
    beta_n_dc = beta_n_dc,
    alpha_p_dc = alpha_p_dc,
    beta_p_dc = beta_p_dc,

    mu = mu,
    l = l,
    min_kappa = min_kappa,
    max_kappa = max_kappa,
    max_sigma_count = max_sigma_count,
    min_sigma_count = min_sigma_count
  )

  index_maps = list(
    sampling_times = sampling_times,
    intermediate_times = if(n_intermediate_times > 0) intermediate_rows$time else numeric(0),
    driver_ids = if(N_driver > 0) driver_muts$driver_id else character(0),
    driver_n_ids = if(N_driver_n > 0) driver_n_muts$driver_n_id else character(0),
    driver_p_ids = if(N_driver_p > 0) driver_p_muts$driver_p_id else character(0),
    dc_ids = if(N_dc > 0) dc_muts$dc_id else character(0),
    clade_wt_ids = if(N_clades_wt > 0) clades_wt$clade_id else character(0),
    epistates = c("-","+")
  )

  list(data = data, index_maps = index_maps)
}
