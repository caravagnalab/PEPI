#' Check clade_statistics input for arbitrary time points
#'
#' @param clade_statistics A tibble or data.frame with clade-level statistics.
#'   Must include columns `cluster_name`, `n_muts`, `is_clade`, `names`, `driver_type`, `has_driver`, `rewind`
#'   and VAF columns named `vaf_k_n` and `vaf_k_p` for k = 1, ..., N.
#'
#' @return Invisibly returns TRUE if checks pass; stops with error otherwise.
#' 
#' @export
check_input_vaf <- function(clade_statistics) {
  required_base <- c("cluster_name", "n_muts", "is_clade", "names",
                     "driver_type", "has_driver", "rewind")
  missing_base <- setdiff(required_base, names(clade_statistics))
  if(length(missing_base) > 0) stop(paste("Missing columns in clade_statistics:", paste(missing_base, collapse=", ")))
  
  # Automatically check all vaf_k_n and vaf_k_p columns
  vaf_cols <- grep("^vaf_\\d+_[np]$", names(clade_statistics), value = TRUE)
  if(length(vaf_cols) == 0) stop("No VAF columns detected (expected names like vaf_1_n, vaf_1_p, ...)")
  
  # check numeric and [0,1]
  for(col in vaf_cols) {
    if(!is.numeric(clade_statistics[[col]])) stop(paste(col,"must be numeric"))
    if(any(clade_statistics[[col]] < 0 | clade_statistics[[col]] > 1)) stop(paste(col,"values out of [0,1] range"))
  }
  
  if(!is.numeric(clade_statistics$n_muts)) stop("n_muts must be numeric")
  if(!is.logical(clade_statistics$is_clade)) stop("is_clade must be logical")
  
  invisible(TRUE)
}

#' Check counts input
#'
#' @param counts A tibble or data.frame with columns `count_n`, `count_p`, `time`
#'
#' @return Invisibly returns TRUE if checks pass; stops with error otherwise.
#'
#' @export
check_input_count <- function(counts) {
  required_cols <- c("count_n", "count_p", "time")
  missing_cols <- setdiff(required_cols, names(counts))
  if(length(missing_cols) > 0) stop(paste("Missing columns in counts:", paste(missing_cols, collapse=", ")))
  
  for(col in required_cols) {
    if(!is.numeric(counts[[col]])) stop(paste(col,"must be numeric"))
  }
  
  invisible(TRUE)
}


# check_genomic_constants.R
#' Check Genomic Constants in a PEPI object
#'
#' This function validates that the genomic constants tibble contains the required
#' variables `genome_length` and `mu` with positive numeric values.
#'
#' @param genomic_constants A tibble with columns `variable` and `value`.
#'
#' @return TRUE if valid; otherwise stops with an error.
#' @export
#'
#' @examples
#' library(tibble)
#' # Valid example
#' genomic_constants <- tibble(
#'   variable = c("genome_length", "mu"),
#'   value = c(3e9, 1e-7)
#' )
#' check_genomic_constants(genomic_constants)
#'
#' # Invalid example: negative mutation rate
#' \dontrun{
#' genomic_constants_bad <- tibble(
#'   variable = c("genome_length", "mu"),
#'   value = c(3e9, -1e-7)
#' )
#' check_genomic_constants(genomic_constants_bad)
#' }
check_genomic_constants <- function(genomic_constants) {
  
  # Check class
  if (!tibble::is_tibble(genomic_constants)) {
    stop("Genomic_constants must be a tibble")
  }
  
  # Check columns
  required_cols <- c("variable", "value")
  missing_cols <- setdiff(required_cols, colnames(genomic_constants))
  if (length(missing_cols) > 0) {
    stop("Genomic_constants tibble must have columns: ", paste(required_cols, collapse = ", "))
  }
  
  # Check required variables
  required_vars <- c("genome_length", "mu")
  missing_vars <- setdiff(required_vars, genomic_constants$variable)
  if (length(missing_vars) > 0) {
    stop("Genomic_constants tibble is missing variables: ", paste(missing_vars, collapse = ", "))
  }
  
  # Check that values are positive
  for (var in required_vars) {
    value <- genomic_constants$value[genomic_constants$variable == var]
    if (!is.numeric(value) || value <= 0) {
      stop(sprintf("Value of %s must be numeric and positive", var))
    }
  }
  
  return(TRUE)
}

