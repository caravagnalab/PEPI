#' Build Stan data list for PEPI model
#'
#' @param pepi A PEPI object.
#'
#' @return Named list for Stan input.
#' @export

get_stan_data_pepi <- function(pepi) {
  
  clade_statistics <- pepi$clade_statistics
  counts <- pepi$Counts
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
    if(nrow(df) == 0) return(array(0, dim = c(0, n_times, 2)))
    arr <- array(0, dim = c(nrow(df), n_times, 2))
    for (j in seq_len(n_times)) {
      arr[, j, 1] <- 2 * df[[paste0("vaf_", j, "_n")]]
      arr[, j, 2] <- 2 * df[[paste0("vaf_", j, "_p")]]
    }
    arr
  }
  
  ## ---------------------------
  ## Clades
  ## ---------------------------
  clades <- clade_statistics %>% filter(is_clade, !has_driver)
  ccf_clade <- extract_ccf(clades)
  m_clade <- clades$m
  
  ## ---------------------------
  ## Driver negative
  ## ---------------------------
  drivers_n <- clade_statistics %>% filter(has_driver, driver_type == "driver_n")
  ccf_driver_n <- extract_ccf(drivers_n)
  m_driver_n <- drivers_n$m
  
  ## ---------------------------
  ## DC and CD
  ## ---------------------------
  build_pair <- function(type) {
    df <- clade_statistics %>% filter(has_driver, driver_type == type)
    list(
      m = matrix(df$m, ncol = 2, byrow = TRUE),
      ccf_driver = extract_ccf(df),
      ccf_clade = extract_ccf(df)
    )
  }
  
  dc <- build_pair("dc")
  cd <- build_pair("cd")
  
  ## ---------------------------
  ## Wild type
  ## ---------------------------
  ccf_wt <- array(1, dim = c(n_times, 2))
  
  ## ---------------------------
  ## Counts and genomic constants
  ## ---------------------------
  zminus <- counts$count_n
  zplus  <- counts$count_p
  
  mu <- genomic_constants %>%
    dplyr::filter(variable == "mutation_rate") %>% pull(value)
  
  l <- genomic_constants %>%
    dplyr::filter(variable == "genome_length") %>% pull(value)
  
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
    
    m_dc = dc$m,
    ccf_dc_driver = dc$ccf_driver,
    ccf_dc_clade = dc$ccf_clade,
    
    m_cd = cd$m,
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
  counts <- pepi$Counts
  
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

