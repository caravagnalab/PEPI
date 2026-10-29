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

#' R version of the Stan approx_Z_clade(), vectorised over posterior draws.
#'
#' Clades start from a single "+" cell (boundary condition), so there is no
#' chi term.
#'
#' @return A list with the "-" (`neg`) and "+" (`pos`) clade sizes.
#' @keywords internal

.approx_Z_clade = function(lambda_n, lambda_p, omega_n, omega_p, t_1, t_2, W1){

  dt = t_2 - t_1

  eps = 1e-10

  delta = sqrt((lambda_p - lambda_n)^2 + eps)

  x1 = 0.5 * (lambda_p + lambda_n) + sqrt(0.25 * (lambda_p - lambda_n)^2 + omega_p * omega_n)

  exp_p = exp(x1 * dt)

  list(neg = (omega_n / delta) * W1 * exp_p,
       pos = W1 * exp_p)

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
#' then generated as in the likelihood of the Stan model: for each population,
#' its total count is drawn from `lognormal(log(Z- + Z+), sigma_count)` and its
#' split in "+" cells from `beta_proportion(Z+ / (Z- + Z+), kappa)`. The
#' private populations driver_n and driver_p live in one compartment and only
#' get the lognormal count, as do wt/dc clades whose observed "-" count is 0.
#' CCFs, within the "-" and "+" samples, have no likelihood of their own and
#' are evaluated on the latent sizes: the count noise is shared by all the
#' populations of a sample and cancels in the CCF;
#' wt clades and dc clades are subpopulations of wt and dc and do not enter
#' the sample sizes. The total population size (`z_tot` of `pop = "total"`)
#' follows the likelihood of the intermediate total counts,
#' `lognormal(log(ztot_lat), sigma_count)`, with `ztot_lat` the latent size
#' of all populations except the clades.
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
#'   `z_neg`, `z_pos` (counts), `z_tot` (total size, only for `"total"`), `frac_pos` (fraction of "+" cells within the
#'   population), `ccf_neg`, `ccf_pos` (CCF in the "-" and "+" samples).
#' @examples
#' \donttest{
#' source(system.file("examples", "fit_example.R", package = "PEPI"))
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

  # populations as list(pop, id, compartment(t), top, size(t) -> list(neg, pos)):
  # `top` = FALSE for the clades, subpopulations that do not enter the totals
  both = function(t) "both"
  clades = list()

  clades[[1]] = list(pop = "wt", id = "wt", compartment = both, top = TRUE, size = function(t){
    .approx_Z(lambda_n, lambda_n * (1 + s_epi), p("omega_n_wt"), p("omega_p_wt"),
              p("tmrca"), t, p("bc_wt"), p("chi_wt"))
  })

  for(j in seq_len(data$N_driver)){
    local({
      j = j
      s = p(paste0("s_driver[", j, "]"))
      clades[[length(clades) + 1]] <<- list(pop = "driver", id = ids("driver", data$N_driver)[j],
                                            compartment = both, top = TRUE, size = function(t){
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
                                            compartment = both, top = TRUE, size = function(t){
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
                                            compartment = function(t) "neg", top = TRUE,
                                            size = function(t){
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
                                            compartment = function(t) "pos", top = TRUE,
                                            size = function(t){
        list(neg = rep(0, ndraws),
             pos = exp(lambda_n * (1 + p(paste0("s_driver_p[", j, "]"))) * (1 + s_epi) *
                         (t - p(paste0("t_driver_p[", j, "]"))) + log(p(paste0("bc_driver_p[", j, "]")))))
      })
    })
  }

  # Clades: as in the likelihood, only the "+" count is observed when the
  # observed "-" count is 0. At the sampling times the observed count at that
  # time is used; at other times, whether it is ever positive.
  sampling_times = data$times[data$sampling_index]

  clade_compartment = function(neg_counts){
    function(t){
      k = which(abs(sampling_times - t) < 1e-10)
      positive = if(length(k) > 0) neg_counts[k[1]] > 0 else any(neg_counts > 0)
      if(positive) "both" else "pos"
    }
  }

  for(j in seq_len(data$N_clades_wt)){
    local({
      j = j
      clades[[length(clades) + 1]] <<- list(pop = "clade_wt", id = ids("clade_wt", data$N_clades_wt)[j],
                                            compartment = clade_compartment(data$ccf_clade_wt[j, , 1] * data$zminus),
                                            top = FALSE, size = function(t){
        .approx_Z_clade(lambda_n, lambda_n * (1 + s_epi), p("omega_n_wt"), p("omega_p_wt"),
                        p(paste0("t_clade_wt[", j, "]")), t, p(paste0("bc_clade_wt[", j, "]")))
      })
    })
  }

  for(j in seq_len(data$N_dc)){
    local({
      j = j
      s = p(paste0("s_dc[", j, "]"))
      clades[[length(clades) + 1]] <<- list(pop = "dc_clade", id = ids("dc", data$N_dc)[j],
                                            compartment = clade_compartment(data$ccf_dc_clade[j, , 1] * data$zminus),
                                            top = FALSE, size = function(t){
        .approx_Z_clade(lambda_n * (1 + s), lambda_n * (1 + s_epi) * (1 + s),
                        p(paste0("omega_n_dc[", j, "]")), p(paste0("omega_p_dc[", j, "]")),
                        p(paste0("t_dc[", j, ",2]")), t, p(paste0("bc_clade_dc[", j, "]")))
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

    Z = lapply(clades, function(cl) cl$size(t))
    z = lapply(seq_along(clades), function(k) observe(Z[[k]], clades[[k]]$compartment(t)))

    top = vapply(clades, function(cl) cl$top, logical(1))

    # sizes of the "-" and "+" samples (named so that they are not masked
    # by the z_neg/z_pos columns inside tibble())
    tot_neg = Reduce(`+`, lapply(z[top], function(v) v$neg))
    tot_pos = Reduce(`+`, lapply(z[top], function(v) v$pos))

    # latent sizes of the "-" and "+" samples, for the CCFs
    tot_neg_lat = Reduce(`+`, lapply(Z[top], function(v) v$neg))
    tot_pos_lat = Reduce(`+`, lapply(Z[top], function(v) v$pos))

    # total population size, as the intermediate counts of the likelihood
    ztot_lat = Reduce(`+`, lapply(Z[top], function(v) v$neg + v$pos))
    z_tot = if(noise) exp(stats::rnorm(ndraws, log(ztot_lat), sigma_count)) else ztot_lat

    res = list(tibble::tibble(samp = seq_len(ndraws), time = t, pop = "total", id = NA_character_,
                              z_neg = tot_neg, z_pos = tot_pos, z_tot = z_tot,
                              frac_pos = tot_pos / (tot_neg + tot_pos),
                              ccf_neg = 1, ccf_pos = 1))

    for(k in seq_along(clades)){
      res[[k + 1]] = tibble::tibble(samp = seq_len(ndraws), time = t, pop = clades[[k]]$pop, id = clades[[k]]$id,
                                    z_neg = z[[k]]$neg, z_pos = z[[k]]$pos,
                                    frac_pos = z[[k]]$pos / (z[[k]]$neg + z[[k]]$pos),
                                    # ccf_neg = z[[k]]$neg / tot_neg,
                                    # ccf_pos = z[[k]]$pos / tot_pos)
                                    ccf_neg = Z[[k]]$neg / tot_neg_lat,
                                    ccf_pos = Z[[k]]$pos / tot_pos_lat)
    }

    dplyr::bind_rows(res)

  }) %>% dplyr::bind_rows()

  out

}


#' Observed data in the format of \code{get_predicted_trajectories()}.
#'
#' Counts of each population in the "-"/"+" samples are its CCF times the
#' size of the sample, as in the Stan model's transformed data; total
#' counts at the intermediate times are added to `pop = "total"`.
#'
#' @return A tibble with columns `time`, `pop`, `id`, `z_neg`, `z_pos`,
#'   `z_tot`, `frac_pos`, `ccf_neg`, `ccf_pos`.
#' @keywords internal

.observed_trajectories = function(x){

  tabs = x$multirates

  pop = tabs$population_sizes %>% dplyr::filter(type == "sampling") %>%
    dplyr::select(time, zminus, zplus)

  total = pop %>% dplyr::transmute(time, pop = "total", id = NA_character_,
                                   z_neg = zminus, z_pos = zplus, z_tot = zminus + zplus,
                                   frac_pos = zplus/(zminus + zplus), ccf_neg = 1, ccf_pos = 1)

  # total counts at the intermediate times
  intermediate = tabs$population_sizes %>% dplyr::filter(type == "intermediate")
  if(nrow(intermediate) > 0){
    total = dplyr::bind_rows(total, intermediate %>% dplyr::transmute(time, pop = "total", id = NA_character_, z_tot = ztot))
  }

  ccf = .multirates_observed_ccf(x) %>%
    dplyr::mutate(group = dplyr::case_when(group == "dc" & quantity == "clade" ~ "dc_clade", TRUE ~ group),
                  id = ifelse(group == "wt", "wt", id),
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
