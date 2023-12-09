#' Check if VAF input data has the correct form.

#' Check if the input data has the correct fields and stop if someone is missing.
#'
#' @param data Dataset to analyze. It must be a dataframe with number of variants and depth 
#'                                  for any mutation and sample.
#' @return PEPI object
#' @examples
#' check_input_vaf(data)
#' @export

check_input_vaf = function(data){
  
    if(! "Nx" %in% colnames(data)){
      
      stop("Missing Nx")
    }
    if(! "Ny" %in% colnames(data)){
      
      stop("Missing Ny")
    }
    if(! "DPx" %in% colnames(data)){
      
      stop("Missing DPx")
    }
    if(! "DPy" %in% colnames(data)){
      
      stop("Missing DPy")
    }
    if(! is.integer(data$Nx) & sum(data$Nx > 0)/nrow(data)){
      
      stop("Nx must be non negative integer")
    }
    if(! is.integer(data$Ny) & sum(data$Ny > 0)/nrow(data)){
      
      stop("Ny must be non negative integer")
    }
    if(! is.integer(data$DPx) & sum(data$DPx > 0)/nrow(data) ){
      
      stop("DPx must be non negative integer")
    }
    if(! is.integer(data$DPy) & sum(data$DPy > 0)/nrow(data)){
      
      stop("DPy must be non negative integer")
    }
    if(! sum(data$Nx <= data$DPx)/nrow(data)){
      
      stop(" Nx <= DPx not true")
    }
    if(! sum(data$Ny <= data$DPy)/nrow(data)){
      
      stop(" Ny <= DPy not true")
    }
    
}



#' Check if counts input data has the correct form.

#' Check if the input data has the correct fields and stop if someone is missing.
#'
#' @param data Dataframe containing with number of cell counts of - and + cells 
#'             at different time points in case of Counts application. 
#' @return PEPI object
#' @examples
#' check_input_counts(data)
#' @export
 
check_input_counts = function(data){
    
    if(! "epistate" %in% colnames(data)){
      
      stop("Missing epistate")
    }
    if(! "counts" %in% colnames(data)){
      
      stop("Missing counts")
    }
    if(! "time" %in% colnames(data)){
      
      stop("Missing time")
    }
    if(! "+" %in% data$epistate & "-" %in% data$epistate){
      
      stop("Epistates should be + and -")
    }
    if(! sum(data$counts >= 0)/nrow(data)){
      
      stop("Counts must be positive")
    }
    
}
  


