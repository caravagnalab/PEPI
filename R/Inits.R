#' Create a PEPI object for the unified multirates model.
#'
#' A PEPI object holding the tables required to build the data list for
#' \code{inst/multirates_positive_s.stan} is created. See
#' \code{check_input_multirates()} for the required table schema.
#'
#' @param tables Named list of tibbles: \code{population_sizes}, \code{wt_ccf}
#'   (required), and any subset of \code{driver_muts}/\code{driver_ccf},
#'   \code{driver_n_muts}/\code{driver_n_ccf}, \code{driver_p_muts}/\code{driver_p_ccf},
#'   \code{dc_muts}/\code{dc_ccf}, \code{clades_wt}/\code{clades_wt_ccf} (optional).
#' @param m_trunk Number of truncal mutations.
#' @return PEPI object of class "PEPI_Multirates"
#' @examples
#' sim = simulate_multirates_tree(
#'   sampling_times = c(3, 6, 9), tmrca = 1, t_min = 0,
#'   lambda_n = 1, s_epi = 0.15, omega_n_wt = 5e-3, omega_p_wt = 2e-3,
#'   clades_wt = list(list(id = "c1", t_clade_wt = 1.5)),
#'   mu = 1e-7, l = 2.7e9, kappa = 20, sigma_count = 0.1, seed = 42
#' )
#' x = init_multirates(sim$tables, m_trunk = sim$truth$m_trunk)
#' @export

init_multirates = function(tables, m_trunk = 0L){

  check_input_multirates(tables)

  pepi = list(multirates = tables, m_trunk = m_trunk)

  class(pepi) = "PEPI_Multirates"

  return(pepi)

}
