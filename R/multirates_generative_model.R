#' Synthetic generative model for the unified multirates model
#'
#' Simulates a deterministic mean-field trajectory (the same 2x2 growth/switch
#' ODE the Stan model itself uses as its likelihood mean) for a wt-background
#' lineage plus any number of driver, driver-only, driver+clade, and
#' wt-background-clade sub-lineages, each founded by a single cell at its own
#' introduction time. CCF/fraction/mutation-count observations are then
#' generated with noise matched to the Stan model's own likelihoods
#' (beta_proportion for CCF/fractions, lognormal for counts, Poisson with the
#' same rate formulas as the model's `generated quantities` block for mutation
#' counts), so that data simulated here is recoverable by `fit_multirates()`.

#' Integrate the 2-epistate growth/switch ODE analytically via matrix exponential.
#'
#' @param lambda_n Growth rate in the "-" epistate.
#' @param s Fitness of "+" relative to "-" (lambda_p = lambda_n*(1+s)).
#' @param omega_n,omega_p Switch rates - to + and + to -.
#' @param t1,t2 Start/end time.
#' @param Z0 Length-2 initial state c(z_minus, z_plus) at t1.
#' @return Length-2 state c(z_minus, z_plus) at t2.
#' @keywords internal

simulate_Z = function(lambda_n, s, omega_n, omega_p, t1, t2, Z0){

  if(t2 <= t1) return(Z0)

  lambda_p = lambda_n * (1 + s)

  A = matrix(c(lambda_n, omega_p, omega_n, lambda_p), nrow = 2)

  as.vector(expm::expm(A * (t2 - t1)) %*% Z0)

}

#' @keywords internal
add_ccf_noise = function(p_true, kappa){

  p_true = min(max(p_true, 1e-6), 1 - 1e-6)

  rbeta(1, p_true * kappa, (1 - p_true) * kappa)

}

#' @keywords internal
add_count_noise = function(z_true, sigma_count){

  exp(rnorm(1, log(max(z_true, 1e-9)), sigma_count))

}


#' Simulate synthetic input tables for `init_multirates()`.
#'
#' @param sampling_times Numeric vector of sampling times.
#' @param tmrca Time of most recent common ancestor.
#' @param t_min Time origin (lower bound of the tmrca prior / start of the wt lineage).
#' @param lambda_n Growth rate of "-" cells.
#' @param s_epi Fitness of "+" vs "-" wt cells.
#' @param omega_n_wt,omega_p_wt wt switch rates.
#' @param drivers List of driver specs, each a list with `id`, `t_driver`, `s_driver`,
#'   `omega_n_driver`, `omega_p_driver`, `ms_driver`, `sigma_driver`,
#'   `alpha_n_driver`, `beta_n_driver`, `alpha_p_driver`, `beta_p_driver`.
#' @param driver_n List of driver-only "-"-confined specs: `id`, `t_driver_n`,
#'   `s_driver_n`, `ms_driver_n`, `sigma_driver_n`.
#' @param driver_p List of driver-only "+"-confined specs: `id`, `t_driver_p`,
#'   `s_driver_p`, `ms_driver_p`, `sigma_driver_p`.
#' @param dc List of driver+clade combo specs: `id`, `delta_t1` (tmrca to driver
#'   acquisition), `delta_t2` (driver acquisition to clade/switch event), `s_dc`,
#'   `omega_n_dc`, `omega_p_dc`, `ms_dc`, `sigma_dc`, `alpha_n_dc`, `beta_n_dc`,
#'   `alpha_p_dc`, `beta_p_dc`.
#' @param clades_wt List of wt-background clade specs: `id`, `t_clade_wt`.
#' @param mu Mutation rate per division per bp per allele.
#' @param l Length of the genome.
#' @param kappa Concentration parameter for the beta_proportion CCF/fraction noise.
#' @param sigma_count Lognormal sd for population-size noise.
#' @param seed Optional RNG seed.
#' @return A list with `tables` (ready for `init_multirates()`) and `truth`
#'   (the ground-truth parameter values used to simulate the data).
#' @examples
#' simulate_multirates_tree(sampling_times = c(2,5,8), tmrca = 1, t_min = 0,
#'   lambda_n = 1, s_epi = 0.2, omega_n_wt = 1e-3, omega_p_wt = 1e-4,
#'   drivers = list(list(id = "KRAS", t_driver = 3, s_driver = 0.3,
#'     omega_n_driver = 1e-3, omega_p_driver = 1e-4, ms_driver = -0.5,
#'     sigma_driver = 0.5, alpha_n_driver = 1, beta_n_driver = 10,
#'     alpha_p_driver = 1, beta_p_driver = 10)),
#'   seed = 1908)
#' @export

simulate_multirates_tree = function(sampling_times, tmrca, t_min, lambda_n, s_epi,
    omega_n_wt, omega_p_wt,
    drivers = list(), driver_n = list(), driver_p = list(), dc = list(), clades_wt = list(),
    mu = 1e-7, l = 2.7e9, kappa = 200, sigma_count = 0.1, seed = NULL){

  if(!is.null(seed)) set.seed(seed)

  sampling_times = sort(unique(sampling_times))

  Z_wt_at = function(t) simulate_Z(lambda_n, s_epi, omega_n_wt, omega_p_wt, t_min, t, c(1,0))

  Z_driver_at = function(spec, t){
    if(t < spec$t_driver) return(c(0,0))
    simulate_Z(lambda_n, spec$s_driver, spec$omega_n_driver, spec$omega_p_driver, spec$t_driver, t, c(1,0))
  }

  Z_driver_n_at = function(spec, t){
    if(t < spec$t_driver_n) return(0)
    exp(lambda_n * (1 + spec$s_driver_n) * (t - spec$t_driver_n))
  }

  Z_driver_p_at = function(spec, t){
    if(t < spec$t_driver_p) return(0)
    exp(lambda_n * (1 + spec$s_driver_p) * (t - spec$t_driver_p))
  }

  Z_dc_driver_at = function(spec, t){
    t_dc1 = tmrca + spec$delta_t1
    if(t < t_dc1) return(c(0,0))
    simulate_Z(lambda_n, spec$s_dc, spec$omega_n_dc, spec$omega_p_dc, t_dc1, t, c(1,0))
  }

  Z_dc_clade_at = function(spec, t){
    t_dc1 = tmrca + spec$delta_t1
    t_dc2 = t_dc1 + spec$delta_t2
    if(t < t_dc2) return(c(0,0))
    simulate_Z(lambda_n, spec$s_dc, spec$omega_n_dc, spec$omega_p_dc, t_dc2, t, c(1,0))
  }

  Z_clade_wt_at = function(spec, t){
    if(t < spec$t_clade_wt) return(c(0,0))
    simulate_Z(lambda_n, s_epi, omega_n_wt, omega_p_wt, spec$t_clade_wt, t, c(1,0))
  }

  # ---- total population per sampling time (sum of all additive lineages) --

  totals = lapply(sampling_times, function(t){

    zwt = Z_wt_at(t)
    neg = zwt[1]; pos = zwt[2]

    for(spec in drivers){ z = Z_driver_at(spec, t); neg = neg + z[1]; pos = pos + z[2] }
    for(spec in driver_n){ neg = neg + Z_driver_n_at(spec, t) }
    for(spec in driver_p){ pos = pos + Z_driver_p_at(spec, t) }
    for(spec in dc){ z = Z_dc_driver_at(spec, t); neg = neg + z[1]; pos = pos + z[2] }

    list(neg = neg, pos = pos, zwt = zwt)
  })
  names(totals) = as.character(sampling_times)

  # ---- population_sizes + wt_ccf -------------------------------------------

  population_sizes = tibble::tibble(
    time = sampling_times,
    type = "sampling",
    zminus = sapply(totals, function(x) add_count_noise(x$neg, sigma_count)),
    zplus = sapply(totals, function(x) add_count_noise(x$pos, sigma_count))
  )

  wt_ccf = dplyr::bind_rows(lapply(sampling_times, function(t){
    tot = totals[[as.character(t)]]
    tibble::tibble(time = t, epistate = c("-","+"),
                   ccf = c(add_ccf_noise(tot$zwt[1]/tot$neg, kappa),
                           add_ccf_noise(tot$zwt[2]/tot$pos, kappa)))
  }))

  # ---- driver ---------------------------------------------------------------

  driver_muts_tbl = NULL; driver_ccf_tbl = NULL
  if(length(drivers) > 0){

    driver_muts_tbl = dplyr::bind_rows(lapply(drivers, function(spec){
      m_driver = rpois(1, 4*mu*l*lambda_n*(spec$t_driver - tmrca))
      tibble::tibble(driver_id = spec$id, m_driver = m_driver,
                     ms_driver = spec$ms_driver, sigma_driver = spec$sigma_driver,
                     alpha_n_driver = spec$alpha_n_driver, beta_n_driver = spec$beta_n_driver,
                     alpha_p_driver = spec$alpha_p_driver, beta_p_driver = spec$beta_p_driver)
    }))

    driver_ccf_tbl = dplyr::bind_rows(lapply(drivers, function(spec){
      dplyr::bind_rows(lapply(sampling_times, function(t){
        tot = totals[[as.character(t)]]
        z = Z_driver_at(spec, t)
        tibble::tibble(driver_id = spec$id, time = t, epistate = c("-","+"),
                       ccf = c(add_ccf_noise(z[1]/tot$neg, kappa), add_ccf_noise(z[2]/tot$pos, kappa)))
      }))
    }))
  }

  # ---- driver_n ---------------------------------------------------------------

  driver_n_tbl = NULL; driver_n_ccf_tbl = NULL
  if(length(driver_n) > 0){

    driver_n_tbl = dplyr::bind_rows(lapply(driver_n, function(spec){
      m = rpois(1, 4*mu*l*lambda_n*(spec$t_driver_n - tmrca))
      tibble::tibble(driver_n_id = spec$id, m_driver_n = m,
                     ms_driver_n = spec$ms_driver_n, sigma_driver_n = spec$sigma_driver_n)
    }))

    driver_n_ccf_tbl = dplyr::bind_rows(lapply(driver_n, function(spec){
      dplyr::bind_rows(lapply(sampling_times, function(t){
        tot = totals[[as.character(t)]]
        zn = Z_driver_n_at(spec, t)
        tibble::tibble(driver_n_id = spec$id, time = t, ccf = add_ccf_noise(zn/tot$neg, kappa))
      }))
    }))
  }

  # ---- driver_p ---------------------------------------------------------------

  driver_p_tbl = NULL; driver_p_ccf_tbl = NULL
  if(length(driver_p) > 0){

    driver_p_tbl = dplyr::bind_rows(lapply(driver_p, function(spec){
      tibble::tibble(driver_p_id = spec$id, ms_driver_p = spec$ms_driver_p, sigma_driver_p = spec$sigma_driver_p)
    }))

    driver_p_ccf_tbl = dplyr::bind_rows(lapply(driver_p, function(spec){
      dplyr::bind_rows(lapply(sampling_times, function(t){
        tot = totals[[as.character(t)]]
        zp = Z_driver_p_at(spec, t)
        tibble::tibble(driver_p_id = spec$id, time = t, ccf = add_ccf_noise(zp/tot$pos, kappa))
      }))
    }))
  }

  # ---- dc ---------------------------------------------------------------

  dc_tbl = NULL; dc_ccf_tbl = NULL
  if(length(dc) > 0){

    dc_tbl = dplyr::bind_rows(lapply(dc, function(spec){
      t_dc1 = tmrca + spec$delta_t1
      t_dc2 = t_dc1 + spec$delta_t2
      m1 = rpois(1, 2*mu*l*lambda_n*2*(t_dc1 - tmrca))
      m2 = rpois(1, 2*mu*l*lambda_n*(1 + spec$s_dc)*2*(t_dc2 - t_dc1))
      tibble::tibble(dc_id = spec$id, m_dc_driver = m1, m_dc_clade = m2,
                     ms_dc = spec$ms_dc, sigma_dc = spec$sigma_dc,
                     alpha_n_dc = spec$alpha_n_dc, beta_n_dc = spec$beta_n_dc,
                     alpha_p_dc = spec$alpha_p_dc, beta_p_dc = spec$beta_p_dc)
    }))

    dc_ccf_tbl = dplyr::bind_rows(lapply(dc, function(spec){
      dplyr::bind_rows(lapply(sampling_times, function(t){
        tot = totals[[as.character(t)]]
        zd = Z_dc_driver_at(spec, t)
        zc = Z_dc_clade_at(spec, t)
        dplyr::bind_rows(
          tibble::tibble(dc_id = spec$id, time = t, epistate = c("-","+"), quantity = "driver",
                         ccf = c(add_ccf_noise(zd[1]/tot$neg, kappa), add_ccf_noise(zd[2]/tot$pos, kappa))),
          tibble::tibble(dc_id = spec$id, time = t, epistate = c("-","+"), quantity = "clade",
                         ccf = c(add_ccf_noise(zc[1]/tot$neg, kappa), add_ccf_noise(zc[2]/tot$pos, kappa)))
        )
      }))
    }))
  }

  # ---- clades_wt ---------------------------------------------------------------

  clades_wt_tbl = NULL; clades_wt_ccf_tbl = NULL
  if(length(clades_wt) > 0){

    clades_wt_tbl = dplyr::bind_rows(lapply(clades_wt, function(spec){
      m_clade_wt = rpois(1, 2*mu*l*lambda_n*2*(spec$t_clade_wt - tmrca))
      tibble::tibble(clade_id = spec$id, m_clade_wt = m_clade_wt)
    }))

    clades_wt_ccf_tbl = dplyr::bind_rows(lapply(clades_wt, function(spec){
      dplyr::bind_rows(lapply(sampling_times, function(t){
        tot = totals[[as.character(t)]]
        z = Z_clade_wt_at(spec, t)
        tibble::tibble(clade_id = spec$id, time = t, epistate = c("-","+"),
                       ccf = c(add_ccf_noise(z[1]/tot$neg, kappa), add_ccf_noise(z[2]/tot$pos, kappa)))
      }))
    }))
  }

  m_trunk = rpois(1, 2*mu*l*lambda_n*2*(tmrca - t_min))

  tables = list(
    population_sizes = population_sizes,
    wt_ccf = wt_ccf,
    driver_muts = driver_muts_tbl, driver_ccf = driver_ccf_tbl,
    driver_n_muts = driver_n_tbl, driver_n_ccf = driver_n_ccf_tbl,
    driver_p_muts = driver_p_tbl, driver_p_ccf = driver_p_ccf_tbl,
    dc_muts = dc_tbl, dc_ccf = dc_ccf_tbl,
    clades_wt = clades_wt_tbl, clades_wt_ccf = clades_wt_ccf_tbl
  )

  truth = list(
    lambda_n = lambda_n, s_epi = s_epi, omega_n_wt = omega_n_wt, omega_p_wt = omega_p_wt,
    tmrca = tmrca, t_min = t_min, m_trunk = m_trunk,
    s_driver = stats::setNames(sapply(drivers, function(s) s$s_driver), sapply(drivers, function(s) s$id)),
    s_driver_n = stats::setNames(sapply(driver_n, function(s) s$s_driver_n), sapply(driver_n, function(s) s$id)),
    s_driver_p = stats::setNames(sapply(driver_p, function(s) s$s_driver_p), sapply(driver_p, function(s) s$id)),
    s_dc = stats::setNames(sapply(dc, function(s) s$s_dc), sapply(dc, function(s) s$id))
  )

  list(tables = tables, truth = truth)

}
