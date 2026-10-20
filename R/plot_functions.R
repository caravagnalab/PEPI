#' plotting functions

#' Plot prior and posterior distributions of the model parameters.
#'
#' One histogram per parameter, prior (transparent) and posterior (solid).
#' Growth parameters (`lambda_n`, `s_*`) are shown below 2, times (`tmrca`,
#' `t_*`) between -60 and the largest sampled `tmrca`, and switch rates
#' (`omega_*`) on a log10 scale. Prior draws of growth parameters beyond the
#' 99% posterior quantile + 1.5 are dropped.
#'
#' @param x PEPI_Multirates object, after \code{get_posterior_multirates()}.
#' @param params Vector of Stan parameter names, with indices for vector
#'   parameters (e.g. `"s_driver[1]"`). If NULL, `lambda_n`, `s_epi`,
#'   `omega_p_wt`, `tmrca` and every `s_driver`, `t_driver`, `omega_p_driver`,
#'   `omega_n_driver` element.
#' @param truth Optional named vector of ground-truth values, named as `params`;
#'   shown as dashed vertical lines.
#' @param ncol,nrow Layout of the panels.
#' @return A ggplot object.
#' @examples
#' \donttest{
#' source(system.file("examples", "fit_example.R", package = "PEPI"))
#' x = get_posterior_multirates(x)
#' plot_prior_posterior(x, truth = c(lambda_n = 1, s_epi = 0.5, "s_driver[1]" = 0.35))
#' }
#' @export

plot_prior_posterior = function(x, params = NULL, truth = NULL, ncol = 3, nrow = 3){

  if(is.null(x$posterior$multirates)){
    stop("run get_posterior_multirates() first")
  }

  draws = x$posterior$multirates %>%
    dplyr::filter(family == "param") %>%
    dplyr::mutate(variable = sub("_prior", "", variable))

  if(is.null(params)){
    all = unique(draws$variable)
    params = c("lambda_n", "s_epi", grep("^s_driver\\[", all, value = TRUE),
               "omega_p_wt", "tmrca", grep("^t_driver\\[", all, value = TRUE),
               grep("^omega_p_driver\\[", all, value = TRUE),
               grep("^omega_n_driver\\[", all, value = TRUE))
  }

  missing = setdiff(params, draws$variable)
  if(length(missing) > 0){
    stop("parameters not found: ", paste(missing, collapse = ", "))
  }

  draws = draws %>% dplyr::filter(variable %in% params)

  # drop prior draws of growth parameters far beyond the posterior
  draws = draws %>%
    dplyr::group_by(variable) %>%
    dplyr::mutate(lim_sup = stats::quantile(value[type == "posterior"], 0.99, na.rm = TRUE) + 1.5,
                  is_growth = variable == "lambda_n" | grepl("^s", variable)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is_growth | type == "posterior" | value <= lim_sup)

  t_max = draws %>% dplyr::filter(variable == "tmrca") %>% dplyr::pull(value)
  t_max = if(length(t_max) > 0) round(max(t_max), digits = 3) else Inf

  base_plot = function(d, p){
    ggplot2::ggplot(d) +
      ggplot2::geom_histogram(ggplot2::aes(x = value, alpha = type), fill = "steelblue",
                              position = "identity", bins = 100) +
      ggplot2::scale_alpha_manual(values = c("posterior" = 1, "prior" = 0.5)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none", plot.title = ggplot2::element_text(size = 9)) +
      ggplot2::labs(x = "", title = p)
  }

  growth = params[grepl("^s_|^lambda_", params)]
  times = params[grepl("^t", params)]
  omega = params[grepl("^omega_", params)]
  other = setdiff(params, c(growth, times, omega))

  plots = c(
    lapply(growth, function(p){
      base_plot(draws %>% dplyr::filter(variable == p, value < 2), p) +
        ggplot2::scale_x_continuous(breaks = round(seq(-1, 1.5, 0.3), digits = 2))
    }),
    lapply(times, function(p){
      base_plot(draws %>% dplyr::filter(variable == p, value < t_max, value > -60), p)
    }),
    lapply(omega, function(p){
      base_plot(draws %>% dplyr::filter(variable == p, value > 1e-10), p) + ggplot2::scale_x_log10()
    }),
    lapply(other, function(p){
      base_plot(draws %>% dplyr::filter(variable == p), p)
    })
  )
  names(plots) = c(growth, times, omega, other)

  if(!is.null(truth)){
    for(p in intersect(names(truth), names(plots))){
      plots[[p]] = plots[[p]] + ggplot2::geom_vline(xintercept = truth[[p]], linetype = "dashed")
    }
  }

  ggpubr::ggarrange(plotlist = plots, nrow = nrow, ncol = ncol)

}


# ---- posterior predictive trajectories ---------------------------------------

#' @keywords internal
.ribbon_plot = function(pred, obs, var, probs, ylab, facet = FALSE){

  group_vars = if(facet) c("facet", "time") else "time"

  pred_summary = pred %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::summarise(
      low  = stats::quantile(.data[[var]], probs[1], na.rm = TRUE),
      mid  = stats::quantile(.data[[var]], 0.50, na.rm = TRUE),
      high = stats::quantile(.data[[var]], probs[2], na.rm = TRUE),
      .groups = "drop"
    )

  p = ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = pred_summary, ggplot2::aes(x = time, ymin = low, ymax = high),
                         fill = "steelblue", alpha = 0.2) +
    ggplot2::geom_line(data = pred_summary, ggplot2::aes(x = time, y = mid),
                       color = "steelblue", linetype = "dashed", linewidth = 1) +
    ggplot2::geom_point(data = obs, ggplot2::aes(x = time, y = .data[[var]]),
                        color = "indianred", size = 1.1) +
    ggplot2::labs(y = ylab, x = "Time") +
    ggplot2::theme_minimal(base_size = 14)

  if(facet) p = p + ggplot2::facet_wrap(~facet)

  p

}

# one facet per population; ids are shown only when a type has several
.pop_facet = function(d){
  d %>% dplyr::group_by(pop) %>%
    dplyr::mutate(facet = ifelse(dplyr::n_distinct(id) > 1, paste0(pop, ": ", id), pop)) %>%
    dplyr::ungroup()
}


#' Plot the posterior predictive fraction of "+" cells in the whole population.
#'
#' Median (dashed line) and 5%-95% interval (band) of the fraction of "+"
#' cells over time, from \code{get_predicted_trajectories()}, against the
#' observed fractions at the sampling times (points).
#'
#' @param x PEPI_Multirates object fit with \code{fit_multirates()}.
#' @param pred Output of \code{get_predicted_trajectories()}; computed if NULL.
#' @return A ggplot object.
#' @examples
#' \donttest{
#' source(system.file("examples", "fit_example.R", package = "PEPI"))
#' plot_fractions(x)
#' }
#' @export

plot_fractions = function(x, pred = NULL){

  if(is.null(pred)) pred = get_predicted_trajectories(x)

  .ribbon_plot(pred %>% dplyr::filter(pop == "total"),
               .observed_trajectories(x) %>% dplyr::filter(pop == "total"),
               var = "frac_pos", probs = c(0.05, 0.95), ylab = "Fraction Z+ in the whole population")

}


#' Plot the posterior predictive total population size.
#'
#' Median (dashed line) and 5%-95% interval (band) of the total number of
#' cells over time (log scale), from \code{get_predicted_trajectories()},
#' against the observed sizes at the sampling and intermediate times (points).
#'
#' @param x PEPI_Multirates object fit with \code{fit_multirates()}.
#' @param pred Output of \code{get_predicted_trajectories()}; computed if NULL.
#' @return A ggplot object.
#' @examples
#' \donttest{
#' source(system.file("examples", "fit_example.R", package = "PEPI"))
#' plot_counts(x)
#' }
#' @export

plot_counts = function(x, pred = NULL){

  if(is.null(pred)) pred = get_predicted_trajectories(x)

  pred = pred %>% dplyr::filter(pop == "total")
  obs = .observed_trajectories(x) %>% dplyr::filter(pop == "total")

  .ribbon_plot(pred, obs, var = "z_tot", probs = c(0.05, 0.95), ylab = "Ztot") +
    ggplot2::scale_y_log10()

}


#' Plot the posterior predictive "+"/"-" split of each population.
#'
#' Median (dashed line) and 5%-95% interval (band) of the fraction of "+"
#' cells within every population (wild type, drivers, ...) over time, from
#' \code{get_predicted_trajectories()}, against the observed split at the
#' sampling times (points).
#'
#' @param x PEPI_Multirates object fit with \code{fit_multirates()}.
#' @param pred Output of \code{get_predicted_trajectories()}; computed if NULL.
#' @param pops Populations to show; default all that are sorted in "-" and
#'   "+" (wt, drivers, dc, wt clades, dc clades).
#' @return A ggplot object.
#' @examples
#' \donttest{
#' source(system.file("examples", "fit_example.R", package = "PEPI"))
#' plot_split(x)
#' }
#' @export

plot_split = function(x, pred = NULL, pops = NULL){

  if(is.null(pred)) pred = get_predicted_trajectories(x)

  if(is.null(pops)) pops = c("wt", "driver", "dc", "clade_wt", "dc_clade")

  .ribbon_plot(pred %>% dplyr::filter(pop %in% pops) %>% .pop_facet(),
               .observed_trajectories(x) %>% dplyr::filter(pop %in% pops) %>% .pop_facet(),
               var = "frac_pos", probs = c(0.05, 0.95), ylab = "Z+ / (Z- + Z+) in population",
               facet = TRUE)

}


#' Plot the posterior predictive CCF of each population in the sorted samples.
#'
#' Median (dashed line) and 2%-98% interval (band) of the CCF of every
#' population (wild type, drivers, ...) in the "-" or "+" sample over time,
#' from \code{get_predicted_trajectories()}, against the observed CCF at the
#' sampling times (points).
#'
#' @param x PEPI_Multirates object fit with \code{fit_multirates()}.
#' @param epistate "-" or "+".
#' @param pred Output of \code{get_predicted_trajectories()}; computed if NULL.
#' @param pops Populations to show; default all except `"total"`.
#' @return A ggplot object.
#' @examples
#' \donttest{
#' source(system.file("examples", "fit_example.R", package = "PEPI"))
#' plot_ccf(x, epistate = "-")
#' }
#' @export

plot_ccf = function(x, epistate = c("-", "+"), pred = NULL, pops = NULL){

  epistate = match.arg(epistate)
  var = ifelse(epistate == "-", "ccf_neg", "ccf_pos")

  if(is.null(pred)) pred = get_predicted_trajectories(x)

  # driver-only clades live in one sample only
  if(is.null(pops)){
    pops = c("wt", "driver", "dc", "clade_wt", "dc_clade", ifelse(epistate == "-", "driver_n", "driver_p"))
  }

  .ribbon_plot(pred %>% dplyr::filter(pop %in% pops) %>% .pop_facet(),
               .observed_trajectories(x) %>% dplyr::filter(pop %in% pops) %>% .pop_facet(),
               var = var, probs = c(0.02, 0.98),
               ylab = paste0("Population fraction in Z", ifelse(epistate == "-", "\u2212", "+"), " sample"),
               facet = TRUE)

}


#' Plot all posterior predictive checks.
#'
#' Arranges \code{plot_fractions()}, \code{plot_counts()},
#' \code{plot_split()} and \code{plot_ccf()} for both samples, all computed
#' from the same posterior predictive trajectories.
#'
#' @param x PEPI_Multirates object fit with \code{fit_multirates()}.
#' @param pred Output of \code{get_predicted_trajectories()}; computed if NULL.
#' @return A ggplot object.
#' @examples
#' \donttest{
#' source(system.file("examples", "fit_example.R", package = "PEPI"))
#' plot_posterior_predictive(x)
#' }
#' @export

plot_posterior_predictive = function(x, pred = NULL){

  if(is.null(pred)) pred = get_predicted_trajectories(x)

  ggpubr::ggarrange(
    ggpubr::ggarrange(plot_fractions(x, pred), plot_counts(x, pred), ncol = 2),
    plot_split(x, pred),
    plot_ccf(x, "-", pred),
    plot_ccf(x, "+", pred),
    nrow = 4, ncol = 1
  )

}
