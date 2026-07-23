#' plotting functions

#' Plot multivariate VAF distributions with cluster associated to tree nodes.
#'
#' A multivariate plot is generated from a labelled dataset
#'
#' @param spectrum VAF spectrum with cluster labels
#' @return A multivariate plot
#' @examples
#' \dontrun{
#' plot_multivariate(spectrum)
#' }
#' @export

plot_multivariate = function(spectrum){


if(!"node" %in% colnames(spectrum)){

  stop("no cluster labels")

}

max_level = spectrum %>% pull(node) %>% nchar() %>% max() - 1
cls = get_colors(max_level)

ggplot(spectrum %>% mutate(vaf_x = Nx/DPx, vaf_y = Ny/DPy)) + geom_point(aes(x = vaf_x, y = vaf_y,color = node)) +
  get_pepi_theme() + scale_colour_manual(values = cls) +
  labs(title = "Multivariate Spectrum", x = "VAF -", y = "VAF +")

}

#' Plot marginal VAF distributions with cluster associated to tree nodes.
#'
#' Two marginal histograms are generated from a labelled dataset.
#'
#' @param spectrum VAF spectrum with cluster labels
#' @return Two marginal histograms.
#' @examples
#' \dontrun{
#' plot_marginal(spectrum)
#' }
#' @export

plot_marginal = function(spectrum){

  if(!"node" %in% colnames(spectrum)){

    stop("no cluster labels")

  }

  max_level = spectrum %>% pull(node) %>% nchar() %>% max() - 1
  cls = get_colors(max_level)

px =   ggplot(spectrum %>% mutate(vaf_x = Nx/DPx) %>% filter(vaf_x > 0)) +
  geom_histogram(aes(x = vaf_x, fill = node), bins = 50) +
    get_pepi_theme() + scale_fill_manual(values = cls) + labs(title = "Marginal -", x = "VAF")

py =   ggplot(spectrum %>% mutate(vaf_y = Ny/DPy) %>% filter(vaf_y > 0)) +
  geom_histogram(aes(x = vaf_y, fill= node), bins = 50) +
  get_pepi_theme() + scale_fill_manual(values = cls) + labs(title = "Marginal +", x = "VAF")

 ggarrange(plotlist = list(px,py),ncol = 2, nrow = 1)

 }

#' Plot predicted vs observed counts for the unified multirates model.
#'
#' Faceted by group/id and epistate, using \code{get_predicted_counts()}.
#'
#' @param x PEPI_Multirates object, after \code{get_posterior_multirates()}.
#' @return Plot of predicted vs observed counts.
#' @examples
#' \dontrun{
#' plot_counts_multirates(x)
#' }
#' @export

plot_counts_multirates = function(x){

  d = get_predicted_counts(x)

  d$facet = ifelse(is.na(d$id), d$group, paste0(d$group, ": ", d$id))

  ggplot(d) +
    geom_point(aes(x = sampling_time, y = z_observed), color = "black") +
    geom_pointrange(aes(x = sampling_time, y = mean, ymin = lower, ymax = upper), color = "dodgerblue") +
    facet_grid(facet ~ epistate, scales = "free_y") +
    get_pepi_theme() +
    labs(x = "time", y = "count", title = "Predicted (blue) vs observed (black) counts")

}

#' Plot predicted vs observed fractions for the unified multirates model.
#'
#' Faceted by group/id and epistate, using \code{get_predicted_fractions()}.
#'
#' @param x PEPI_Multirates object, after \code{get_posterior_multirates()}.
#' @return Plot of predicted vs observed fractions.
#' @examples
#' \dontrun{
#' plot_fractions(x)
#' }
#' @export

plot_fractions = function(x){

  d = get_predicted_fractions(x)

  d$facet = ifelse(is.na(d$id), d$group, paste0(d$group, ": ", d$id))

  ggplot(d) +
    geom_point(aes(x = sampling_time, y = frac_observed), color = "black") +
    geom_pointrange(aes(x = sampling_time, y = mean, ymin = lower, ymax = upper), color = "dodgerblue") +
    facet_grid(facet ~ epistate, scales = "free_y") +
    get_pepi_theme() +
    labs(x = "time", y = "fraction", title = "Predicted (blue) vs observed (black) fractions")

}

#' Plot predicted vs observed CCF for the unified multirates model.
#'
#' Faceted by group/id and epistate, using \code{get_predicted_ccf()}.
#'
#' @param x PEPI_Multirates object, after \code{get_posterior_multirates()}.
#' @return Plot of predicted vs observed CCF.
#' @examples
#' \dontrun{
#' plot_ccf(x)
#' }
#' @export

plot_ccf = function(x){

  d = get_predicted_ccf(x)

  d$facet = ifelse(is.na(d$id), d$group, paste0(d$group, ": ", d$id))

  ggplot(d) +
    geom_point(aes(x = sampling_time, y = ccf_observed), color = "black") +
    geom_pointrange(aes(x = sampling_time, y = mean, ymin = lower, ymax = upper), color = "dodgerblue") +
    facet_grid(facet ~ epistate, scales = "free_y") +
    get_pepi_theme() +
    labs(x = "time", y = "CCF", title = "Predicted (blue) vs observed (black) CCF")

}

#' Plot predicted vs observed mutation counts for the unified multirates model.
#'
#' Posterior-predictive distributions (violin) vs observed counts (point),
#' using \code{get_posterior_multirates()}'s "m" family.
#'
#' @param x PEPI_Multirates object, after \code{get_posterior_multirates()}.
#' @return Plot of predicted vs observed mutation counts.
#' @examples
#' \dontrun{
#' plot_mutations(x)
#' }
#' @export

plot_mutations = function(x){

  if(is.null(x$posterior$multirates)){
    stop("run get_posterior_multirates() first")
  }

  draws = x$posterior$multirates %>% dplyr::filter(family == "m", type == "posterior")

  draws = dplyr::left_join(draws, .multirates_observed_mutations(x), by = c("group","id","event"))

  draws$label = ifelse(is.na(draws$id), draws$group,
                       ifelse(is.na(draws$event), paste0(draws$group, ": ", draws$id),
                              paste0(draws$group, ": ", draws$id, " (", draws$event, ")")))

  ggplot(draws) +
    geom_violin(aes(x = label, y = value), fill = "steelblue", alpha = 0.5) +
    geom_point(aes(x = label, y = m_observed), color = "black", size = 2) +
    coord_flip() +
    get_pepi_theme() +
    labs(x = NULL, y = "mutation count", title = "Predicted (violin) vs observed (point) mutation counts")

}

#' Plot posterior and prior distributions of the unified multirates model's parameters.
#'
#' Posterior and prior draws histograms are plotted for any required parameter,
#' using the raw ("param" family) rows of \code{get_posterior_multirates()}'s output.
#'
#' @param x PEPI_Multirates object, after \code{get_posterior_multirates()}.
#' @param params A vector of canonical parameter names (e.g. "lambda_n", "s_driver").
#' @param groups A vector of groups to include (wt/driver/driver_n/driver_p/dc/clade_wt/global).
#' @return A plot with posterior and prior distributions.
#' @examples
#' \dontrun{
#' plot_inference(x,params = c("lambda_n","s_epi"))
#' }
#' @export

plot_inference = function(x, params = NULL, groups = NULL){

  if(is.null(x$posterior$multirates)){
    stop("run get_posterior_multirates() first")
  }

  d = x$posterior$multirates %>% dplyr::filter(family == "param")

  d$canonical = sub("_prior$", "", d$base)

  if(!is.null(params)){
    d = d %>% dplyr::filter(canonical %in% params)
  }
  if(!is.null(groups)){
    d = d %>% dplyr::filter(group %in% groups)
  }

  if(nrow(d) == 0){
    stop("required parameters are not present")
  }

  d$facet = ifelse(is.na(d$id), d$canonical, paste0(d$canonical, "[", d$id, "]"))

  nr = round(length(unique(d$facet))/4 + 1)
  nc = min(length(unique(d$facet)) + 1, 4)

  ggplot(d) +
    geom_histogram(aes(x = value, alpha = type, fill = group), bins = 40) +
    facet_wrap(~facet, scales = "free", nrow = nr, ncol = nc) +
    scale_alpha_manual(values = c("prior" = 0.4, "posterior" = 1)) +
    get_pepi_theme() + theme(legend.position = "bottom")

}
