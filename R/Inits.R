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
#' \dontrun{
#' init_multirates(tables, m_trunk = 120)
#' }
#' @export

init_multirates = function(tables, m_trunk = 0L){

  check_input_multirates(tables)

  pepi = list(multirates = tables, m_trunk = m_trunk)

  class(pepi) = "PEPI_Multirates"

  return(pepi)

}
