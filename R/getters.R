# utilities inference

# ---- variable lookup table for the unified multirates model -----------------
# Maps every Stan variable base name (as it appears before "[...]" in
# fit$draws() column names) to the sequence of index "kinds" its brackets
# encode, so get_posterior_multirates() can parse draws into a tidy table.
# Kinds: "id" (driver/dc/clade id, via index_maps), "time" (sampling time),
# "itime" (intermediate time), "epistate" (- / +), "event" (raw 1/2 index,
# e.g. the two legs of a dc combo). Variables not listed here (e.g. internal
# ODE latent states like Z_wt_lat) are intentionally dropped.

.multirates_dims = local({

  m = list()

  reg = function(bases, dims){
    for(b in bases) m[[b]] <<- dims
  }

  reg(c("lambda_n","s_epi","omega_p_wt","omega_n_wt","bc_wt","U_wt","W_wt","tmrca","kappa","sigma_count",
        "lambda_n_prior","tmrca_prior","s_epi_prior","omega_n_wt_prior","omega_p_wt_prior","m_trunk_pred"),
      character(0))

  reg(c("t_clade_wt","bc_clade_wt","t_clade_wt_prior","m_clade_wt_pred"), "id")

  reg(c("s_driver","t_driver","bc_driver","U_driver","W_driver","omega_p_driver","omega_n_driver",
        "s_driver_prior","omega_n_driver_prior","omega_p_driver_prior","t_driver_prior","m_driver_pred"), "id")

  reg(c("s_driver_n","t_driver_n","bc_driver_n","s_driver_n_prior","t_driver_n_prior",
        "omega_p_driver_n_pred","m_driver_n_pred"), "id")

  reg(c("s_driver_p","t_driver_p","bc_driver_p","s_driver_p_prior","t_driver_p_prior",
        "omega_n_driver_p_pred"), "id")

  reg(c("s_dc","bc_driver_dc","bc_clade_dc","U_dc","W_dc","omega_p_dc","omega_n_dc",
        "s_dc_prior","omega_n_dc_prior","omega_p_dc_prior"), "id")

  reg(c("delta_t_dc","t_dc_prior","m_dc_pred"), c("id","event"))

  reg("z_wt_pred", c("time","epistate"))
  reg("z_driver_pred", c("id","time","epistate"))
  reg(c("z_dc_driver_pred","z_dc_clade_pred"), c("id","time","epistate"))
  reg("z_driver_n_pred", c("id","time"))
  reg("z_driver_p_pred", c("id","time"))
  reg("z_clade_wt_pred", c("id","time","epistate"))
  reg("ztot_pred", "itime")

  reg("frac_wt_pred", c("time","epistate"))
  reg("frac_driver_pred", c("id","time","epistate"))
  reg(c("frac_dc_driver_pred","frac_dc_clade_pred"), c("id","time","epistate"))
  reg("frac_clade_wt_pred", c("id","time","epistate"))
  reg(c("frac_pos_pred","frac_neg_pred"), "time")

  reg("ccf_driver_pred", c("id","time","epistate"))
  reg("ccf_driver_n_pred", c("id","time"))
  reg("ccf_driver_p_pred", c("id","time"))
  reg("ccf_driver_dc_pred", c("id","time","epistate"))

  m
})

.multirates_family = function(base){

  if(base == "ztot_pred") return("z")
  if(grepl("^ccf_", base)) return("ccf")
  if(grepl("^frac_", base)) return("frac")
  if(grepl("^z_", base)) return("z")
  if(grepl("^m_", base)) return("m")
  "param"

}

.multirates_group = function(base){

  if(grepl("clade_wt", base)) return("clade_wt")
  if(grepl("driver_n(_|$)", base)) return("driver_n")
  if(grepl("driver_p(_|$)", base)) return("driver_p")
  if(grepl("dc", base)) return("dc")
  if(grepl("driver", base)) return("driver")
  if(grepl("wt", base)) return("wt")
  "global"

}

#' Extract and tidy the full posterior (and built-in prior draws) from a
#' multirates fit.
#'
#' Every recognized Stan variable (raw parameters and all `generated
#' quantities` families: z/frac/ccf/m posterior predictives, plus the
#' inline `*_prior` draws) is parsed into one long tidy table with columns
#' `type` (prior/posterior), `family` (z/frac/ccf/m/param), `group`
#' (wt/driver/driver_n/driver_p/dc/clade_wt/global), `id`, `sampling_time`,
#' `epistate`, `event`, `quantity` (driver/clade leg for dc entities),
#' `base`, `variable`, draw indices, and `value`.
#'
#' @param x PEPI_Multirates object with `x$inference$multirates` populated
#'   by \code{fit_multirates()}.
#' @return `x` with `x$posterior$multirates` populated.
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
#'   method = "variational", ndraws = 500, seed = 45,
#'   mu = 1e-7, l = 2.7e9, t_min = 0, ms_epi = 0, sigma_epi = 0.5,
#'   alpha_lambda = 1, beta_lambda = 1, alpha_n_wt = 1, beta_n_wt = 10,
#'   alpha_p_wt = 1, beta_p_wt = 10, min_kappa = 5, max_kappa = 200,
#'   min_sigma_count = 0.01, max_sigma_count = 1)
#' x = get_posterior_multirates(x)
#' }
#' @export

get_posterior_multirates = function(x){

  fit = x$inference$multirates
  if(is.null(fit)){
    stop("no multirates inference")
  }

  index_maps = x$index_maps

  draws_df = posterior::as_draws_df(fit$draws()) %>% tibble::as_tibble()

  long = draws_df %>%
    tidyr::pivot_longer(-dplyr::any_of(c(".chain",".iteration",".draw")), names_to = "variable", values_to = "value")

  base_vec = sub("\\[.*\\]$", "", long$variable)
  idx_vec = ifelse(grepl("\\[", long$variable), sub("^[^\\[]*\\[(.*)\\]$", "\\1", long$variable), NA_character_)

  long$base = base_vec
  long$idx = idx_vec

  long = long[long$base %in% names(.multirates_dims), ]

  if(nrow(long) == 0){
    stop("no recognized multirates variables found in fit$draws()")
  }

  long$type = ifelse(grepl("_prior$", long$base), "prior", "posterior")
  long$family = vapply(long$base, .multirates_family, character(1))
  long$group = vapply(long$base, .multirates_group, character(1))

  parsed = dplyr::bind_rows(lapply(split(long, long$base), function(df){

    base = df$base[1]
    dims = .multirates_dims[[base]]
    group = df$group[1]

    df$id = NA_character_
    df$sampling_time = NA_real_
    df$epistate = NA_character_
    df$event = NA_integer_

    if(length(dims) > 0){

      idx_mat = do.call(rbind, strsplit(df$idx, ","))
      storage.mode(idx_mat) = "integer"

      for(k in seq_along(dims)){

        kind = dims[k]
        col = idx_mat[,k]

        if(kind == "id"){
          ids = index_maps[[paste0(group, "_ids")]]
          df$id = if(!is.null(ids) && length(ids) > 0) ids[col] else as.character(col)
        }else if(kind == "time"){
          ts = index_maps$sampling_times
          df$sampling_time = if(!is.null(ts) && length(ts) > 0) ts[col] else col
        }else if(kind == "itime"){
          its = index_maps$intermediate_times
          df$sampling_time = if(!is.null(its) && length(its) > 0) its[col] else col
        }else if(kind == "epistate"){
          df$epistate = c("-","+")[col]
        }else if(kind == "event"){
          df$event = col
        }
      }
    }

    df

  }))

  parsed = parsed %>%
    dplyr::mutate(
      epistate = dplyr::case_when(
        base == "frac_pos_pred" ~ "+",
        base == "frac_neg_pred" ~ "-",
        TRUE ~ epistate
      ),
      quantity = dplyr::case_when(
        group != "dc" ~ NA_character_,
        grepl("clade", base) ~ "clade",
        family == "m" & event == 2 ~ "clade",
        TRUE ~ "driver"
      )
    ) %>%
    dplyr::select(type, family, group, id, sampling_time, epistate, event, quantity,
                  base, variable, dplyr::any_of(c(".chain",".iteration",".draw")), value)

  x$posterior$multirates = parsed

  return(x)

}


# ---- observed-data helpers (used by the get_predicted_*() wrappers) --------

.multirates_observed_ccf = function(x){

  tabs = x$multirates
  pieces = list(wt = tabs$wt_ccf %>% dplyr::mutate(group = "wt", id = NA_character_))

  if(!is.null(tabs$driver_ccf)){
    pieces$driver = tabs$driver_ccf %>% dplyr::rename(id = driver_id) %>% dplyr::mutate(group = "driver")
  }
  if(!is.null(tabs$driver_n_ccf)){
    pieces$driver_n = tabs$driver_n_ccf %>% dplyr::rename(id = driver_n_id) %>%
      dplyr::mutate(group = "driver_n", epistate = "-")
  }
  if(!is.null(tabs$driver_p_ccf)){
    pieces$driver_p = tabs$driver_p_ccf %>% dplyr::rename(id = driver_p_id) %>%
      dplyr::mutate(group = "driver_p", epistate = "+")
  }
  if(!is.null(tabs$clades_wt_ccf)){
    pieces$clade_wt = tabs$clades_wt_ccf %>% dplyr::rename(id = clade_id) %>% dplyr::mutate(group = "clade_wt")
  }
  if(!is.null(tabs$dc_ccf)){
    pieces$dc = tabs$dc_ccf %>% dplyr::rename(id = dc_id) %>% dplyr::mutate(group = "dc")
  }

  out = dplyr::bind_rows(pieces)

  if(!"quantity" %in% colnames(out)){
    out$quantity = NA_character_
  }

  out

}

.multirates_observed_counts = function(x){

  pop = x$multirates$population_sizes %>% dplyr::filter(type == "sampling") %>%
    dplyr::select(time, zminus, zplus)

  .multirates_observed_ccf(x) %>%
    dplyr::left_join(pop, by = "time") %>%
    dplyr::mutate(z_observed = ifelse(epistate == "-", zminus, zplus)) %>%
    dplyr::transmute(group, id, sampling_time = time, epistate, quantity, z_observed)

}

.multirates_observed_fractions = function(x){

  z_obs = .multirates_observed_counts(x)

  grp = z_obs %>%
    dplyr::group_by(group, id, sampling_time, quantity) %>%
    dplyr::mutate(total = sum(z_observed)) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(group, id, sampling_time, epistate, quantity, frac_observed = z_observed/total)

  pop = x$multirates$population_sizes %>% dplyr::filter(type == "sampling")

  glob = dplyr::bind_rows(
    pop %>% dplyr::transmute(group = "global", id = NA_character_, sampling_time = time,
                             epistate = "+", quantity = NA_character_, frac_observed = zplus/(zminus+zplus)),
    pop %>% dplyr::transmute(group = "global", id = NA_character_, sampling_time = time,
                             epistate = "-", quantity = NA_character_, frac_observed = zminus/(zminus+zplus))
  )

  dplyr::bind_rows(grp, glob)

}

.multirates_observed_mutations = function(x){

  tabs = x$multirates
  pieces = list(m_trunk = tibble::tibble(group = "global", id = NA_character_, event = NA_integer_, m_observed = x$m_trunk))

  if(!is.null(tabs$driver_muts)){
    pieces$driver = tabs$driver_muts %>% dplyr::transmute(group = "driver", id = driver_id, event = NA_integer_, m_observed = m_driver)
  }
  if(!is.null(tabs$driver_n_muts)){
    pieces$driver_n = tabs$driver_n_muts %>% dplyr::transmute(group = "driver_n", id = driver_n_id, event = NA_integer_, m_observed = m_driver_n)
  }
  if(!is.null(tabs$dc_muts)){
    pieces$dc1 = tabs$dc_muts %>% dplyr::transmute(group = "dc", id = dc_id, event = 1L, m_observed = m_dc_driver)
    pieces$dc2 = tabs$dc_muts %>% dplyr::transmute(group = "dc", id = dc_id, event = 2L, m_observed = m_dc_clade)
  }
  if(!is.null(tabs$clades_wt)){
    pieces$clade_wt = tabs$clades_wt %>% dplyr::transmute(group = "clade_wt", id = clade_id, event = NA_integer_, m_observed = m_clade_wt)
  }

  dplyr::bind_rows(pieces)

}

#' Get posterior-predictive counts (z_*_pred) with observed counts joined in.
#'
#' @param x PEPI_Multirates object, after \code{get_posterior_multirates()}.
#' @return A tibble with posterior mean/95%CI and the matching observed count.
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
#'   method = "variational", ndraws = 500, seed = 45,
#'   mu = 1e-7, l = 2.7e9, t_min = 0, ms_epi = 0, sigma_epi = 0.5,
#'   alpha_lambda = 1, beta_lambda = 1, alpha_n_wt = 1, beta_n_wt = 10,
#'   alpha_p_wt = 1, beta_p_wt = 10, min_kappa = 5, max_kappa = 200,
#'   min_sigma_count = 0.01, max_sigma_count = 1)
#' x = get_posterior_multirates(x)
#' get_predicted_counts(x)
#' }
#' @export

get_predicted_counts = function(x){

  if(is.null(x$posterior$multirates)){
    stop("run get_posterior_multirates() first")
  }

  pred = x$posterior$multirates %>%
    dplyr::filter(family == "z", type == "posterior") %>%
    dplyr::group_by(group, id, sampling_time, epistate, quantity) %>%
    dplyr::summarize(mean = mean(value), lower = stats::quantile(value,0.025), upper = stats::quantile(value,0.975), .groups = "drop")

  dplyr::left_join(pred, .multirates_observed_counts(x), by = c("group","id","sampling_time","epistate","quantity"))

}

#' Get posterior-predictive fractions (frac_*_pred) with observed fractions joined in.
#'
#' @param x PEPI_Multirates object, after \code{get_posterior_multirates()}.
#' @return A tibble with posterior mean/95%CI and the matching observed fraction.
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
#'   method = "variational", ndraws = 500, seed = 45,
#'   mu = 1e-7, l = 2.7e9, t_min = 0, ms_epi = 0, sigma_epi = 0.5,
#'   alpha_lambda = 1, beta_lambda = 1, alpha_n_wt = 1, beta_n_wt = 10,
#'   alpha_p_wt = 1, beta_p_wt = 10, min_kappa = 5, max_kappa = 200,
#'   min_sigma_count = 0.01, max_sigma_count = 1)
#' x = get_posterior_multirates(x)
#' get_predicted_fractions(x)
#' }
#' @export

get_predicted_fractions = function(x){

  if(is.null(x$posterior$multirates)){
    stop("run get_posterior_multirates() first")
  }

  pred = x$posterior$multirates %>%
    dplyr::filter(family == "frac", type == "posterior") %>%
    dplyr::group_by(group, id, sampling_time, epistate, quantity) %>%
    dplyr::summarize(mean = mean(value), lower = stats::quantile(value,0.025), upper = stats::quantile(value,0.975), .groups = "drop")

  dplyr::left_join(pred, .multirates_observed_fractions(x), by = c("group","id","sampling_time","epistate","quantity"))

}

#' Get posterior-predictive CCF (ccf_*_pred) with observed CCF joined in.
#'
#' @param x PEPI_Multirates object, after \code{get_posterior_multirates()}.
#' @return A tibble with posterior mean/95%CI and the matching observed CCF.
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
#'   method = "variational", ndraws = 500, seed = 45,
#'   mu = 1e-7, l = 2.7e9, t_min = 0, ms_epi = 0, sigma_epi = 0.5,
#'   alpha_lambda = 1, beta_lambda = 1, alpha_n_wt = 1, beta_n_wt = 10,
#'   alpha_p_wt = 1, beta_p_wt = 10, min_kappa = 5, max_kappa = 200,
#'   min_sigma_count = 0.01, max_sigma_count = 1)
#' x = get_posterior_multirates(x)
#' get_predicted_ccf(x)
#' }
#' @export

get_predicted_ccf = function(x){

  if(is.null(x$posterior$multirates)){
    stop("run get_posterior_multirates() first")
  }

  pred = x$posterior$multirates %>%
    dplyr::filter(family == "ccf", type == "posterior") %>%
    dplyr::group_by(group, id, sampling_time, epistate, quantity) %>%
    dplyr::summarize(mean = mean(value), lower = stats::quantile(value,0.025), upper = stats::quantile(value,0.975), .groups = "drop")

  observed = .multirates_observed_ccf(x) %>%
    dplyr::transmute(group, id, sampling_time = time, epistate, quantity, ccf_observed = ccf)

  dplyr::left_join(pred, observed, by = c("group","id","sampling_time","epistate","quantity"))

}

#' Get posterior-predictive mutation counts (m_*_pred) with observed counts joined in.
#'
#' @param x PEPI_Multirates object, after \code{get_posterior_multirates()}.
#' @return A tibble with posterior mean/95%CI and the matching observed mutation count.
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
#'   method = "variational", ndraws = 500, seed = 45,
#'   mu = 1e-7, l = 2.7e9, t_min = 0, ms_epi = 0, sigma_epi = 0.5,
#'   alpha_lambda = 1, beta_lambda = 1, alpha_n_wt = 1, beta_n_wt = 10,
#'   alpha_p_wt = 1, beta_p_wt = 10, min_kappa = 5, max_kappa = 200,
#'   min_sigma_count = 0.01, max_sigma_count = 1)
#' x = get_posterior_multirates(x)
#' get_predicted_mutations(x)
#' }
#' @export

get_predicted_mutations = function(x){

  if(is.null(x$posterior$multirates)){
    stop("run get_posterior_multirates() first")
  }

  pred = x$posterior$multirates %>%
    dplyr::filter(family == "m", type == "posterior") %>%
    dplyr::group_by(group, id, event) %>%
    dplyr::summarize(mean = mean(value), lower = stats::quantile(value,0.025), upper = stats::quantile(value,0.975), .groups = "drop")

  dplyr::left_join(pred, .multirates_observed_mutations(x), by = c("group","id","event"))

}


#' Clusters mutation with VIBER, a package that implements a variational
#' non parametric Bayesian model to fit multi-variate Binomial mixtures
#'
#' Clusters with centroids VAF coordiates and relative proportion are given.
#'
#' @param data Mutation data with number of variants and depth for any mutation and sample
#' @param K Maximum number of clusters
#' @param alpha Dirichlet concentration parameter
#' @param samples Number of fits computed by the algorithm
#' @param pi_cutoff Cutoff on mixing proportions to filter clusters and reassign mutations
#' @return a tibble with clusters names, mixing proportion and vaf coordinates
#' @examples
#' \dontrun{
#' # requires the VIBER package (github.com/caravagn/VIBER), not on CRAN
#' set.seed(1)
#' data = data.frame(
#'   Nx = rbinom(50, 100, 0.3), DPx = rep(100, 50),
#'   Ny = rbinom(50, 100, 0.1), DPy = rep(100, 50)
#' )
#' get_viber_clusters(data, K = 10, alpha = 1, samples = 1, pi_cutoff = 0.01)
#' }

get_viber_clusters = function(data,K = 10,alpha = 10,samples = 1,pi_cutoff = 0.01){

  DPs =  data %>% dplyr::select(DPx,DPy) %>% rename(X = DPx, Y = DPy)
  NVs =  data %>% dplyr::select(Nx,Ny) %>% rename(X = Nx, Y = Ny)


  fit = VIBER::variational_fit(
    NVs,
    DPs,
    K = K,
    samples = samples,
    alpha_0 = alpha
  )

  fit = VIBER::choose_clusters(
    fit,
    binomial_cutoff = 0,
    dimensions_cutoff = 0,
    pi_cutoff = pi_cutoff,
    re_assign = TRUE
  )

  pi = data.frame(cluster = names(fit$pi_k),pi = fit$pi_k %>% as.vector())
  vaf = fit$theta_k %>% t %>% as.data.frame() %>% dplyr::rename(VAFx = X, VAFy = Y)
  vaf = vaf %>% as.data.frame() %>% mutate(cluster = rownames(vaf)) %>% as_tibble()

  x = full_join(as_tibble(pi),vaf,by = "cluster")  %>% mutate(m = pi*nrow(data))

  return(x)

}

#' Provide a initialization parameters for a PEPI VAF fit
#'
#' Return a list of initialization values for some parameters through some euristics.
#'
#' @param spectrum Mutation data with number of variants and depth for any mutation and sample
#' @param K Maximum number of clusters
#' @param alpha Dirichlet concentration parameter
#' @param samples Number of fits computed by the algorithm
#' @param pi_cutoff Cutoff on mixing proportions to filter clusters and reassign mutations
#' @return a list with initialization values
#' @examples
#' \dontrun{
#' # requires the VIBER package (github.com/caravagn/VIBER), not on CRAN
#' set.seed(1)
#' spectrum = data.frame(
#'   Nx = rbinom(50, 100, 0.3), DPx = rep(100, 50),
#'   Ny = rbinom(50, 100, 0.1), DPy = rep(100, 50)
#' )
#' get_init_values(spectrum, K = 10, alpha = 10, samples = 1, pi_cutoff = 0.01)
#' }
#' @export

get_init_values = function(spectrum,K = 10,alpha = 10,samples = 1,pi_cutoff = 0.01){

  cl =  get_viber_clusters(spectrum,K = K,alpha = alpha,
                     samples = samples, pi_cutoff = pi_cutoff)

  tr = cl %>% arrange(desc(VAFx + VAFy))
  nu_n = tr[1,]$pi
  vaf_minus_n = tr[1,]$VAFx
  vaf_plus_n = tr[1,]$VAFy
  claids_x = cl %>% filter(VAFy < 0.01 & VAFx > 0.01)
  claids_y = cl %>% filter(VAFy > 0.01 & VAFx < 0.01)
  shared = cl %>% filter(VAFx > 0.01 & VAFy > 0.01 &
                           cluster != tr[1,]$cluster)

  if( abs(vaf_plus_n - max(claids_y$VAFy)) < 0.05){

    w_minus_nn = max(claids_x$VAFx)/vaf_minus_n
    w_plus_nn = 0
    nu_nn = claids_x %>% filter(VAFx == max(claids_x$VAFx)) %>% pull(pi)
    rn = 1/(nrow(spectrum))

   }else{

    m_nn = ifelse(max(claids_x$VAFx) > max(claids_y$VAFy),
                  shared$m %>% min(),shared$m %>% max())
    nu_nn = shared %>% filter(m == m_nn) %>% pull(pi)
    vaf_minus_nn = shared %>% filter(m == m_nn) %>% pull(VAFx)
    vaf_plus_nn = shared %>% filter(m == m_nn) %>% pull(VAFy)

     w_minus_nn =  vaf_minus_nn/vaf_minus_n
     w_plus_nn = vaf_plus_nn/vaf_plus_n
     rn = 1/m_nn
}

  if(abs(vaf_minus_n - max(claids_x$VAFx)) < 0.05){

    w_plus_np = max(claids_y$VAFy)/vaf_plus_n
    w_minus_np = 0
    nu_np = claids_y %>% filter(VAFy == max(claids_y$VAFy)) %>% pull(pi)
    rp = 1/(nrow(spectrum))

  }else{

    m_np = ifelse(max(claids_y$VAFy) > max(claids_x$VAFx),
                  shared$m %>% min(),shared$m %>% max())
    nu_np = shared %>% filter(m == m_np) %>% pull(pi)
    vaf_minus_np = shared %>% filter(m == m_np) %>% pull(VAFx)
    vaf_plus_np = shared %>% filter(m == m_np) %>% pull(VAFy)

    w_minus_np =  vaf_minus_np/vaf_minus_n
    w_plus_np = vaf_plus_np/vaf_plus_n
    rp = 1/m_np

 }

 return(list(list( nu_n = nu_n,
                   nu_nn = nu_nn,
                   nu_np = nu_np,
                   rn = rn,
                   rp = rp,
                   vaf_minus_n = vaf_minus_n,
                   vaf_plus_n = vaf_plus_n,
                   w_minus_nn = w_minus_nn,
                   w_plus_nn = w_plus_nn,
                   w_minus_np = w_minus_np,
                   w_plus_np =  w_plus_np)))

}



#' Associate colors to nodes.
#'
#' A list of colors labelled by nodes is generated.
#'
#' @param max_depth Maximal number of levels of the tree
#' @return Named list of colors.
#' @keywords internal
#' @examples
#' library(dplyr)
#' PEPI:::get_colors(max_depth = 2)

get_colors = function(max_depth){

  tree = data.frame(node = "-",level = 0)

  for(l in 1:max_depth){

    epsilon = tree %>% filter(level == l-1)
    node = epsilon %>% pull(node)

    new =  lapply(1:length(node), function(i){

      new_nodes = tibble(node = c(paste0(node[i],"-"),paste0(node[i],"+")), level = l)

    }) %>% bind_rows()

    tree = rbind(tree,new)

  }

  nodes = tree %>% pull(node)
  cls = ggsci::pal_igv()(nodes %>% length())
  names(cls) = nodes

  return(cls)

}


#' Generate personalized ggplot theme.
#'
#' A ggplot theme is generated.
#'
#' @return a ggplot theme.
#' @keywords internal
#' @examples
#' PEPI:::get_pepi_theme()

get_pepi_theme = function(){

  ggplot2::theme_light(base_size = 10) +
    ggplot2::theme(legend.position = "bottom",
      legend.key.size = ggplot2::unit(0.3, "cm"),
      panel.background = ggplot2::element_rect(fill = "white"))

}
