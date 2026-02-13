#' Build Stan data list for PEPI model
#'
#' @param pepi A PEPI object.
#'
#' @return Named list for Stan input.
#' @export

get_stan_data_pepi <- function(pepi) {
  
  clade_statistics <- pepi$clade_statistics
  counts <- pepi$counts
  genomic_constants <- pepi$genomic_constants
  
  n_times <- get_n_times(pepi)
  
  ## ---------------------------
  ## Counts of each type
  ## ---------------------------
  N_driver_n <- clade_statistics %>%
    dplyr::filter(has_driver, driver_type == "driver_n") %>%
    dplyr::pull(names) %>% unique() %>% length()
  
  N_dc <- clade_statistics %>%
    dplyr::filter(has_driver, driver_type == "dc") %>%
    dplyr::pull(names) %>% unique() %>% length()
  
  N_cd <- clade_statistics %>%
    dplyr::filter(has_driver, driver_type == "cd") %>%
    dplyr::pull(names) %>% unique() %>% length()
  
  N_clades <- clade_statistics %>%
    dplyr::filter(is_clade, !has_driver) %>%
    dplyr::pull(names) %>% unique() %>% length()
  
  ## ---------------------------
  ## Helper to extract VAF arrays
  ## ---------------------------
  extract_ccf <- function(df) {
    if (nrow(df) == 0) return(array(0, dim = c(0, n_times, 2)))
    arr <- array(0, dim = c(nrow(df), n_times, 2))
    for (j in seq_len(n_times)) {
      col_n <- paste0("vaf_", j, "_n")
      col_p <- paste0("vaf_", j, "_p")
      if (!(col_n %in% colnames(df)) || !(col_p %in% colnames(df))) {
        stop("VAF columns missing: ", col_n, " or ", col_p)
      }
      arr[, j, 1] <- 2 * df[[col_n]]
      arr[, j, 2] <- 2 * df[[col_p]]
    }
    arr
  }
  
  ## ---------------------------
  ## Clades
  ## ---------------------------
  if (N_clades > 0) {
    clades <- clade_statistics %>% dplyr::filter(is_clade, !has_driver)
    ccf_clade <- extract_ccf(clades)
    m_clade <- clades$n_muts
  } else {
    ccf_clade <- vector(mode="numeric", length=0)
    m_clade <- vector(mode="numeric", length=0)
  }
  
  ## ---------------------------
  ## Driver negative
  ## ---------------------------
  if (N_driver_n > 0) {
    drivers_n <- clade_statistics %>% dplyr::filter(has_driver, driver_type == "driver_n")
    ccf_driver_n <- extract_ccf(drivers_n)
    m_driver_n <- drivers_n$n_muts
  } else {
    ccf_driver_n <- vector(mode="numeric", length=0)
    m_driver_n <- vector(mode="numeric", length=0)
  }
  
  ## ---------------------------
  ## DC and CD
  ## ---------------------------
  build_pair <- function(type) {
    df <- clade_statistics %>% dplyr::filter(has_driver, driver_type == type)
    if (nrow(df) == 0) {
      list(
        n_muts = vector(mode="numeric", length=0),
        ccf_driver = vector(mode="numeric", length=0),
        ccf_clade = vector(mode="numeric", length=0)
      )
    } else {
      ccf_d <- extract_ccf(df)
      ccf_c <- extract_ccf(df)
      # CD: clade minus driver
      if (type == "cd") {
        ccf_c <- pmax(0, ccf_c - ccf_d)
      }
      list(
        n_muts = matrix(df$n_muts, ncol = 2, byrow = TRUE),
        ccf_driver = ccf_d,
        ccf_clade = ccf_c
      )
    }
  }
  
  dc <- build_pair("dc")
  cd <- build_pair("cd")
  
  ## ---------------------------
  ## Wild type (subtract only existing drivers)
  ## ---------------------------
  ccf_wt <- array(1, dim = c(n_times, 2))
  
  subtract_drivers <- function(ccf_array) {
    if (length(ccf_array) > 0) {
      colSums_ccf <- apply(ccf_array, c(2, 3), sum)
      ccf_wt <<- pmax(0, ccf_wt - colSums_ccf)
    }
  }
  
  subtract_drivers(ccf_driver_n)
  subtract_drivers(dc$ccf_driver)
  subtract_drivers(cd$ccf_driver)
  
  ## ---------------------------
  ## Counts and genomic constants
  ## ---------------------------
  zminus <- counts$count_n
  zplus  <- counts$count_p
  
  mu <- genomic_constants %>%
    dplyr::filter(variable == "mu") %>% dplyr::pull(value)
  
  l <- genomic_constants %>%
    dplyr::filter(variable == "genome_length") %>% dplyr::pull(value)
  
  ## ---------------------------
  ## Return Stan-ready list
  ## ---------------------------
  list(
    n_times = n_times,
    N_clades = N_clades,
    N_driver_n = N_driver_n,
    N_dc = N_dc,
    N_cd = N_cd,
    
    m_clade = m_clade,
    ccf_clade = ccf_clade,
    
    m_driver_n = m_driver_n,
    ccf_driver_n = ccf_driver_n,
    
    m_dc = dc$n_muts,
    ccf_dc_driver = dc$ccf_driver,
    ccf_dc_clade = dc$ccf_clade,
    
    m_cd = cd$n_muts,
    ccf_cd_driver = cd$ccf_driver,
    ccf_cd_clade = cd$ccf_clade,
    
    ccf_wt = ccf_wt,
    zminus = zminus,
    zplus = zplus,
    times = counts$time,
    
    mu = mu,
    l = l
  )
}




#' Infer number of time points from a PEPI object
#'
#' @param pepi A PEPI object.
#'
#' @return Integer number of time points.
#' @export
get_n_times <- function(pepi) {
  
  clade_statistics <- pepi$clade_statistics
  counts <- pepi$counts
  
  vaf_cols <- grep("^vaf_[0-9]+_[np]$", colnames(clade_statistics), value = TRUE)
  
  if (length(vaf_cols) == 0) {
    stop("No VAF columns found in clade_statistics.")
  }
  
  n_times <- length(vaf_cols) / 2
  
  if (!n_times %% 1 == 0) {
    stop("Uneven number of VAF columns (expected vaf_k_n / vaf_k_p pairs).")
  }
  
  if (nrow(counts) != n_times) {
    stop(
      "Mismatch between inferred number of time points (", n_times,
      ") and number of rows in counts (", nrow(counts), ")."
    )
  }
  
  as.integer(n_times)
}

#' Extract and save posterior in a PEPI object
#'
#' @param pepi A PEPI object containing a Stan fit (`fit`).
#' @return The updated PEPI object with posterior saved in pepi$posterior.
#' @export
get_posterior <- function(pepi) {
  if (is.null(pepi$fit)) stop("Stan fit not found in the PEPI object.")
  
  sm <- pepi$fit
  n_chains <- sm$num_chains()
  n_draws <- nrow(as_tibble(sm$draws()))
  
  posterior <- lapply(1:n_chains, function(i) {
    ps <- as_tibble(sm$draws()) %>% dplyr::select(starts_with(paste0(i)))
    colnames(ps) <- gsub(paste0(i, "\\."), "", colnames(ps))
    ps <- ps %>% mutate(iter = 1:n_draws)
    
    ps %>%
      reshape2::melt() %>%
      filter(!grepl("prior", variable), variable != "lp__") %>%
      mutate(chain = i, type = "posterior")
  }) %>% bind_rows()
  
  pepi$posterior <- posterior
  return(pepi)
}

#' Extract and save prior in a PEPI object
#'
#' @param pepi A PEPI object containing a Stan fit (`fit`).
#' @return The updated PEPI object with prior saved in pepi$prior.
#' @export
get_prior <- function(pepi) {
  if (is.null(pepi$fit)) stop("Stan fit not found in the PEPI object.")
  
  sm <- pepi$fit
  n_chains <- sm$num_chains()
  n_draws <- nrow(as_tibble(sm$draws()))
  
  prior <- lapply(1:n_chains, function(i) {
    ps <- as_tibble(sm$draws()) %>% dplyr::select(starts_with(paste0(i)))
    colnames(ps) <- gsub(paste0(i, "\\."), "", colnames(ps))
    ps <- ps %>% mutate(iter = 1:n_draws)
    
    ps %>%
      reshape2::melt() %>%
      filter(grepl("prior", variable)) %>%
      mutate(variable = gsub("_prior", "", variable),
             chain = i,
             type = "prior")
  }) %>% bind_rows()
  
  pepi$prior <- prior
  return(pepi)
}

