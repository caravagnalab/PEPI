#' Check if multirates input data has the correct form.
#'
#' Validates a named list of tibbles against the schema required to build
#' the data list for \code{inst/multirates_positive_s.stan}: required tables
#' and columns are present, every id referenced in a `_ccf` table exists in
#' its metadata table (and vice versa), every `time` referenced in a `_ccf`/
#' `wt_ccf` table exists in `population_sizes` with the right `type`, and
#' every id x time (x epistate) combination is present exactly once.
#'
#' @param tables Named list of tibbles, see \code{init_multirates()}.
#' @return Invisibly TRUE if valid, otherwise stops with a descriptive message.
#' @examples
#' \dontrun{
#' check_input_multirates(tables)
#' }
#' @export

check_input_multirates = function(tables){

  if(!is.list(tables)){
    stop("tables must be a named list of tibbles")
  }

  # ---- required tables ----------------------------------------------------

  if(is.null(tables$population_sizes)){
    stop("Missing required table: population_sizes")
  }
  if(is.null(tables$wt_ccf)){
    stop("Missing required table: wt_ccf")
  }

  req_cols = function(df, cols, name){
    missing = setdiff(cols, colnames(df))
    if(length(missing) > 0){
      stop(paste0("Table '", name, "' is missing columns: ", paste(missing, collapse = ", ")))
    }
  }

  req_cols(tables$population_sizes, c("time","type"), "population_sizes")

  if(!all(tables$population_sizes$type %in% c("sampling","intermediate"))){
    stop("population_sizes$type must be one of 'sampling' or 'intermediate'")
  }
  if(any(duplicated(tables$population_sizes$time))){
    stop("population_sizes has duplicated time values")
  }

  sampling_rows = tables$population_sizes %>% dplyr::filter(type == "sampling")

  if(nrow(sampling_rows) == 0){
    stop("population_sizes must contain at least one row with type == 'sampling'")
  }

  req_cols(sampling_rows, c("zminus","zplus"), "population_sizes (sampling rows)")

  if(any(is.na(sampling_rows$zminus)) || any(is.na(sampling_rows$zplus))){
    stop("population_sizes: zminus/zplus must be provided for every 'sampling' row")
  }

  intermediate_rows = tables$population_sizes %>% dplyr::filter(type == "intermediate")
  if(nrow(intermediate_rows) > 0){
    req_cols(intermediate_rows, "ztot", "population_sizes (intermediate rows)")
    if(any(is.na(intermediate_rows$ztot))){
      stop("population_sizes: ztot must be provided for every 'intermediate' row")
    }
  }

  tmax_row = tables$population_sizes %>% dplyr::filter(time == max(time))
  if(tmax_row$type != "sampling"){
    stop("the row with the largest time in population_sizes must have type == 'sampling'")
  }

  req_cols(tables$wt_ccf, c("time","epistate","ccf"), "wt_ccf")

  sampling_times = sampling_rows$time

  check_id_time_epistate = function(muts, ccf, id_col, name, has_epistate = TRUE, extra_muts_cols = c()){

    if(is.null(muts) & is.null(ccf)) return(invisible(TRUE))

    if(is.null(muts) | is.null(ccf)){
      stop(paste0(name, ": both metadata and ccf tables must be supplied together"))
    }

    req_cols(muts, c(id_col, extra_muts_cols), paste0(name, "_muts"))
    ccf_cols = c(id_col,"time","ccf", if(has_epistate) "epistate")
    req_cols(ccf, ccf_cols, paste0(name, "_ccf"))

    ids = muts[[id_col]]

    if(any(duplicated(ids))){
      stop(paste0(name, "_muts: duplicated ", id_col))
    }

    unknown_ids = setdiff(unique(ccf[[id_col]]), ids)
    if(length(unknown_ids) > 0){
      stop(paste0(name, "_ccf references unknown ", id_col, ": ", paste(unknown_ids, collapse = ", ")))
    }

    unknown_times = setdiff(unique(ccf$time), sampling_times)
    if(length(unknown_times) > 0){
      stop(paste0(name, "_ccf references time values not present as 'sampling' rows in population_sizes: ",
                  paste(unknown_times, collapse = ", ")))
    }

    expected = expand.grid(
      id_val = ids,
      time = sampling_times,
      epistate = if(has_epistate) c("-","+") else NA_character_,
      stringsAsFactors = FALSE
    )
    colnames(expected)[1] = id_col
    if(!has_epistate) expected$epistate = NULL
    expected = tibble::as_tibble(expected)

    missing_rows = dplyr::anti_join(expected, ccf, by = colnames(expected))
    if(nrow(missing_rows) > 0){
      stop(paste0(name, "_ccf is missing ", nrow(missing_rows), " required row(s), e.g.:\n",
                  paste(utils::capture.output(print(utils::head(missing_rows))), collapse = "\n")))
    }

    key = do.call(paste, c(as.list(ccf[colnames(expected)]), sep = ""))
    if(any(duplicated(key))){
      stop(paste0(name, "_ccf has duplicated id/time", if(has_epistate) "/epistate", " combinations"))
    }

    invisible(TRUE)
  }

  check_id_time_epistate(tables$driver_muts, tables$driver_ccf, "driver_id", "driver",
                          extra_muts_cols = c("m_driver","ms_driver","sigma_driver",
                                              "alpha_n_driver","beta_n_driver","alpha_p_driver","beta_p_driver"))

  check_id_time_epistate(tables$driver_n_muts, tables$driver_n_ccf, "driver_n_id", "driver_n", has_epistate = FALSE,
                          extra_muts_cols = c("m_driver_n","ms_driver_n","sigma_driver_n"))

  check_id_time_epistate(tables$driver_p_muts, tables$driver_p_ccf, "driver_p_id", "driver_p", has_epistate = FALSE,
                          extra_muts_cols = c("ms_driver_p","sigma_driver_p"))

  check_id_time_epistate(tables$clades_wt, tables$clades_wt_ccf, "clade_id", "clades_wt",
                          extra_muts_cols = c("m_clade_wt"))

  if(!is.null(tables$dc_muts) | !is.null(tables$dc_ccf)){

    if(is.null(tables$dc_muts) | is.null(tables$dc_ccf)){
      stop("dc: both dc_muts and dc_ccf must be supplied together")
    }

    req_cols(tables$dc_muts, c("dc_id","m_dc_driver","m_dc_clade","ms_dc","sigma_dc",
                               "alpha_n_dc","beta_n_dc","alpha_p_dc","beta_p_dc"), "dc_muts")
    req_cols(tables$dc_ccf, c("dc_id","time","epistate","quantity","ccf"), "dc_ccf")

    if(!all(tables$dc_ccf$quantity %in% c("driver","clade"))){
      stop("dc_ccf$quantity must be one of 'driver' or 'clade'")
    }

    ids = tables$dc_muts$dc_id
    if(any(duplicated(ids))) stop("dc_muts: duplicated dc_id")

    unknown_ids = setdiff(unique(tables$dc_ccf$dc_id), ids)
    if(length(unknown_ids) > 0){
      stop(paste0("dc_ccf references unknown dc_id: ", paste(unknown_ids, collapse = ", ")))
    }

    expected = tidyr::expand_grid(dc_id = ids, time = sampling_times, epistate = c("-","+"), quantity = c("driver","clade"))
    missing_rows = dplyr::anti_join(expected, tables$dc_ccf, by = colnames(expected))
    if(nrow(missing_rows) > 0){
      stop(paste0("dc_ccf is missing ", nrow(missing_rows), " required row(s), e.g.:\n",
                  paste(utils::capture.output(print(utils::head(missing_rows))), collapse = "\n")))
    }
  }

  invisible(TRUE)
}
