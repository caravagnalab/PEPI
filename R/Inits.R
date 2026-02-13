#' Initialize a PEPI object
#'
#' This function initializes a PEPI object containing clade-level VAF statistics,
#' population counts across time points, and genomic constants used by the model.
#'
#' @param clade_statistics Tibble containing clade-level mutation and VAF information.
#' @param counts Tibble containing population counts over time.
#' @param genomic_constants Tibble with columns `variable` and `value`,
#'   including at least `genome_length` and `mu`.
#'
#' @return An object of class \code{"PEPI"}.
#' @export
#'
#' @examples
#' library(tibble)
#'
#' clade_statistics <- tibble(
#'   cluster_name = "S",
#'   n_muts = 400,
#'   vaf_1_n = 0,
#'   vaf_1_p = 0.25,
#'   vaf_2_n = 0,
#'   vaf_2_p = 0.28,
#'   is_clade = TRUE,
#'   names = "first clade",
#'   driver_type = "",
#'   has_driver = FALSE,
#'   rewind = FALSE
#' )
#'
#' counts <- tibble(
#'   count_n = c(1e5, 1.2e5),
#'   count_p = c(7e4, 9e4),
#'   time = c(5, 6)
#' )
#'
#' genomic_constants <- tibble(
#'   variable = c("genome_length", "mu"),
#'   value = c(3e9, 1e-7)
#' )
#'
#' pepi <- init(clade_statistics, counts, genomic_constants)

init <- function(clade_statistics, counts, genomic_constants) {
  
  check_input_vaf(clade_statistics)
  check_input_count(counts)
  check_genomic_constants(genomic_constants)
  
  pepi <- list(
    clade_statistics = clade_statistics,
    counts = counts,
    genomic_constants = genomic_constants
  )
  
  class(pepi) <- "PEPI"
  return(pepi)
}


