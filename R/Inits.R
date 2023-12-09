#' Create a PEPI object of VAF type.
#'
#' A PEPI object vaf type is created.
#'
#' @param data Dataset to analyze. It must be a dataframe with number of variants and depth for any mutation and sample
#' @return PEPI object
#' @examples
#' init_vaf(data)
#' @export


init_vaf = function(data){
  
  check_input_vaf(data)
  
  pepi = list(VAF = data)
  
  class(pepi) = "PEPI_VAF"
  
  return(pepi)
  
}

#' Create a PEPI object of Counts type.
#'
#' A PEPI object of vaf type or counts type is created.
#'
#' @param data Dataframe containing with number of cell counts of - and + cells at different time points
#' @return PEPI object
#' @examples
#' init_counts(data)
#' @export


init_counts = function(data){
  
  check_input_counts(data)
  
  pepi = list(counts = data)
  
  class(pepi) = "PEPI_Counts"
  
  return(pepi)
  
}


