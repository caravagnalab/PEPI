# Posterior predictive trajectories for the unified multirates model.
#
# Population sizes are recomputed in R from the posterior draws over an
# arbitrary time grid, using the same closed-form approximation of the
# growth/switch process the Stan model uses (approx_Z() in
# inst/multirates_positive_s.stan). All functions are vectorised over draws.

#' R version of the Stan approx_Z(), vectorised over posterior draws.
#'
#' @return A list with the "-" (`neg`) and "+" (`pos`) population sizes.
#' @keywords internal

.approx_Z = function(lambda_n, lambda_p, omega_n, omega_p, t_1, t_2, W1, chi){

  dt = t_2 - t_1

  eps = 1e-10

  # smooth spectral gap
  delta = sqrt((lambda_p - lambda_n)^2 + eps)

  zeta = lambda_n / lambda_p

  # smooth sine factor
  sin_safe = sqrt(sin(pi * zeta)^2 + eps)

  # dominant eigenvalue
  x1 = 0.5 * (lambda_p + lambda_n) + sqrt(0.25 * (lambda_p - lambda_n)^2 + omega_p * omega_n)

  log_W2 = (log(pi) + log(omega_p + eps) + log(W1 + eps) - log(lambda_p + eps) - log(sin_safe)) / zeta +
    log(chi)

  W2 = exp(log_W2)

  exp_n = exp(lambda_n * dt)
  exp_p = exp(x1 * dt)

  list(neg = W1 * exp_n + (omega_n / delta) * W2 * exp_p,
       pos = W2 * exp_p)

}

#' @keywords internal
.draw_col = function(draws, name){

  if(!name %in% colnames(draws)){
    stop("variable ", name, " not found in the posterior draws")
  }
  as.numeric(draws[[name]])

}


#' Get posterior predictive trajectories of the sorted populations.
#'
#' For every posterior draw and every time of a grid, the latent "-" and "+"
#' sizes of the wild type and of every driver (and dc/driver-only) clade are
#' computed with the model's growth/switch approximation. Observations are
#' then generated as in the likelihood: for each population, its total count
#' is drawn from `lognormal(log(Z- + Z+), sigma_count)` and its split in "+"
#' cells from `beta_proportion(Z+ / (Z- + Z+), kappa)` (driver-only clades
#' live in one compartment and only get the lognormal count). CCFs are
#' finally evaluated on the sampled counts, within the "-" and "+" samples.
#'
#' @param x PEPI_Multirates object fit with \code{fit_multirates()}.
#' @param times Time grid. If NULL, `n_grid` equally spaced times from
#'   `t_start` to the last observed time + 0.25, plus the observed times.
#' @param n_grid Number of points of the default time grid.
#' @param t_start First time of the default time grid.
#' @param noise If TRUE (default), sample the observations; if FALSE, return
#'   the latent trajectories.
#' @return A tibble with one row per draw (`samp`), `time` and population
#'   (`pop`, `id`; `pop = "total"` is the whole population) and columns
#'   `z_neg`, `z_pos` (counts), `frac_pos` (fraction of "+" cells within the
#'   population), `ccf_neg`, `ccf_pos` (CCF in the "-" and "+" samples).
#' @examples
#' \donttest{
#' sim = simulate_multirates_tree(
#'   sampling_times = c(3, 6, 9), tmrca = 1, t_min = 0,
#'   lambda_n = 1, s_epi = 0.15, omega_n_wt = 5e-3, omega_p_wt = 2e-3,
#'   drivers = list(list(id = "KRAS", t_driver = 2, s_driver = 0.3,
#'       omega_n_driver = 5e-3, omega_p_driver = 2e-3,
#'       ms_driver = -0.5, sigma_driver = 0.5,
#'       alpha_n_driver = 1, beta_n_driver = 10, alpha_p_driver = 1, beta_p_driver = 10)),
#'   mu = 1e-7, l = 2.7e9, kappa = 20, sigma_count = 0.1, seed = 42
#' )
#' x = init_multirates(sim$tables, m_trunk = sim$truth$m_trunk)
#' x = fit_multirates(x, chains = 2, ndraws = 200, seed = 45,
#'   mu = 1e-7, l = 2.7e9, t_min = 0, ms_epi = 0, sigma_epi = 0.5,
#'   alpha_lambda = 1, beta_lambda = 1, alpha_n_wt = 1, beta_n_wt = 10,
#'   alpha_p_wt = 1, beta_p_wt = 10, min_kappa = 5, max_kappa = 200,
#'   min_sigma_count = 0.01, max_sigma_count = 1)
#' get_predicted_trajectories(x)
#' }
#' @export

get_predicted_trajectories = function(x, times = NULL, n_grid = 30, t_start = 0, noise = TRUE){

  fit = x$inference$multirates
  data = x$stan_data$multirates

  if(is.null(fit)){
    stop("no multirates inference")
  }

  if(is.null(times)){
    times = sort(c(seq(from = t_start, to = max(data$times) + 0.25, length.out = n_grid),
                   data$times[data$sampling_index]))
  }

  draws = posterior::as_draws_df(fit$draws()) %>% tibble::as_tibble()

  p = function(name) .draw_col(draws, name)

  ndraws = nrow(draws)
  lambda_n = p("lambda_n")
  s_epi = p("s_epi")
  kappa = p("kappa")
  sigma_count = p("sigma_count")

  ids = function(type, n){
    v = x$index_maps[[paste0(type, "_ids")]]
    if(length(v) == n) v else as.character(seq_len(n))
  }

  # clades as list(pop, id, compartment, size(t) -> list(neg, pos))
  clades = list()

  clades[[1]] = list(pop = "wt", id = "wt", compartment = "both", size = function(t){
    .approx_Z(lambda_n, lambda_n * (1 + s_epi), p("omega_n_wt"), p("omega_p_wt"),
              p("tmrca"), t, p("bc_wt"), p("chi_wt"))
  })

  for(j in seq_len(data$N_driver)){
    local({
      j = j
      s = p(paste0("s_driver[", j, "]"))
      clades[[length(clades) + 1]] <<- list(pop = "driver", id = ids("driver", data$N_driver)[j],
                                            compartment = "both", size = function(t){
        .approx_Z(lambda_n * (1 + s), lambda_n * (1 + s_epi) * (1 + s),
                  p(paste0("omega_n_driver[", j, "]")), p(paste0("omega_p_driver[", j, "]")),
                  p(paste0("t_driver[", j, "]")), t,
                  p(paste0("bc_driver[", j, "]")), p(paste0("chi_driver[", j, "]")))
      })
    })
  }

  for(j in seq_len(data$N_dc)){
    local({
      j = j
      s = p(paste0("s_dc[", j, "]"))
      clades[[length(clades) + 1]] <<- list(pop = "dc", id = ids("dc", data$N_dc)[j],
                                            compartment = "both", size = function(t){
        .approx_Z(lambda_n * (1 + s), lambda_n * (1 + s_epi) * (1 + s),
                  p(paste0("omega_n_dc[", j, "]")), p(paste0("omega_p_dc[", j, "]")),
                  p(paste0("t_dc[", j, ",1]")), t,
                  p(paste0("bc_driver_dc[", j, "]")), p(paste0("chi_dc[", j, "]")))
      })
    })
  }

  for(j in seq_len(data$N_driver_n)){
    local({
      j = j
      clades[[length(clades) + 1]] <<- list(pop = "driver_n", id = ids("driver_n", data$N_driver_n)[j],
                                            compartment = "neg", size = function(t){
        list(neg = exp(lambda_n * (1 + p(paste0("s_driver_n[", j, "]"))) * (t - p(paste0("t_driver_n[", j, "]"))) +
                         log(p(paste0("bc_driver_n[", j, "]")))),
             pos = rep(0, ndraws))
      })
    })
  }

  for(j in seq_len(data$N_driver_p)){
    local({
      j = j
      clades[[length(clades) + 1]] <<- list(pop = "driver_p", id = ids("driver_p", data$N_driver_p)[j],
                                            compartment = "pos", size = function(t){
        list(neg = rep(0, ndraws),
             pos = exp(lambda_n * (1 + p(paste0("s_driver_p[", j, "]"))) * (1 + s_epi) *
                         (t - p(paste0("t_driver_p[", j, "]"))) + log(p(paste0("bc_driver_p[", j, "]")))))
      })
    })
  }

  # observation model of one population: lognormal total count, beta split
  observe = function(Z, compartment){

    if(!noise){
      return(Z)
    }

    if(compartment == "neg"){
      return(list(neg = exp(stats::rnorm(ndraws, log(Z$neg), sigma_count)), pos = rep(0, ndraws)))
    }
    if(compartment == "pos"){
      return(list(neg = rep(0, ndraws), pos = exp(stats::rnorm(ndraws, log(Z$pos), sigma_count))))
    }

    tot = exp(stats::rnorm(ndraws, log(Z$neg + Z$pos), sigma_count))
    loc = pmin(1 - 1e-9, pmax(1e-9, Z$pos / (Z$neg + Z$pos)))
    split = stats::rbeta(ndraws, loc * kappa, (1 - loc) * kappa)

    list(neg = tot * (1 - split), pos = tot * split)

  }

  out = lapply(times, function(t){

    z = lapply(clades, function(cl) observe(cl$size(t), cl$compartment))

    # sizes of the "-" and "+" samples (named so that they are not masked
    # by the z_neg/z_pos columns inside tibble())
    tot_neg = Reduce(`+`, lapply(z, function(v) v$neg))
    tot_pos = Reduce(`+`, lapply(z, function(v) v$pos))

    res = list(tibble::tibble(samp = seq_len(ndraws), time = t, pop = "total", id = NA_character_,
                              z_neg = tot_neg, z_pos = tot_pos,
                              frac_pos = tot_pos / (tot_neg + tot_pos),
                              ccf_neg = 1, ccf_pos = 1))

    for(k in seq_along(clades)){
      res[[k + 1]] = tibble::tibble(samp = seq_len(ndraws), time = t, pop = clades[[k]]$pop, id = clades[[k]]$id,
                                    z_neg = z[[k]]$neg, z_pos = z[[k]]$pos,
                                    frac_pos = z[[k]]$pos / (z[[k]]$neg + z[[k]]$pos),
                                    ccf_neg = z[[k]]$neg / tot_neg,
                                    ccf_pos = z[[k]]$pos / tot_pos)
    }

    dplyr::bind_rows(res)

  }) %>% dplyr::bind_rows()

  out

}


#' Observed data in the format of \code{get_predicted_trajectories()}.
#'
#' Counts of each population in the "-"/"+" samples are its CCF times the
#' size of the sample, as in the Stan model's transformed data.
#'
#' @return A tibble with columns `time`, `pop`, `id`, `z_neg`, `z_pos`,
#'   `frac_pos`, `ccf_neg`, `ccf_pos`.
#' @keywords internal

.observed_trajectories = function(x){

  tabs = x$multirates

  pop = tabs$population_sizes %>% dplyr::filter(type == "sampling") %>%
    dplyr::select(time, zminus, zplus)

  total = pop %>% dplyr::transmute(time, pop = "total", id = NA_character_,
                                   z_neg = zminus, z_pos = zplus,
                                   frac_pos = zplus/(zminus + zplus), ccf_neg = 1, ccf_pos = 1)

  ccf = .multirates_observed_ccf(x) %>%
    dplyr::filter(group %in% c("wt", "driver", "dc", "driver_n", "driver_p"),
                  is.na(quantity) | quantity == "driver") %>%
    dplyr::mutate(id = ifelse(group == "wt", "wt", id),
                  epistate = ifelse(epistate == "+", "ccf_pos", "ccf_neg")) %>%
    dplyr::select(time, pop = group, id, epistate, ccf) %>%
    tidyr::pivot_wider(names_from = epistate, values_from = ccf)

  if(!"ccf_pos" %in% colnames(ccf)) ccf$ccf_pos = NA_real_
  if(!"ccf_neg" %in% colnames(ccf)) ccf$ccf_neg = NA_real_

  ccf = ccf %>%
    dplyr::mutate(ccf_pos = ifelse(is.na(ccf_pos), 0, ccf_pos),
                  ccf_neg = ifelse(is.na(ccf_neg), 0, ccf_neg)) %>%
    dplyr::left_join(pop, by = "time") %>%
    dplyr::mutate(z_neg = ccf_neg * zminus, z_pos = ccf_pos * zplus,
                  frac_pos = z_pos / (z_neg + z_pos)) %>%
    dplyr::select(time, pop, id, z_neg, z_pos, frac_pos, ccf_neg, ccf_pos)

  dplyr::bind_rows(total, ccf)

}
