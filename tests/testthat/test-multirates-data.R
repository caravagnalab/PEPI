tiny_tables = function(){

  pop = tibble::tibble(
    time = c(2,5,9),
    type = "sampling",
    zminus = c(100,150,300),
    zplus = c(50,120,250)
  )

  wt_ccf = tidyr::expand_grid(time = pop$time, epistate = c("-","+")) %>%
    dplyr::mutate(ccf = ifelse(epistate == "-", 0.7, 0.6))

  driver_muts = tibble::tibble(driver_id = "KRAS", m_driver = 40, ms_driver = -0.5, sigma_driver = 0.5,
                               alpha_n_driver = 1, beta_n_driver = 10, alpha_p_driver = 1, beta_p_driver = 10)
  driver_ccf = tidyr::expand_grid(driver_id = "KRAS", time = pop$time, epistate = c("-","+")) %>%
    dplyr::mutate(ccf = 0.2)

  clades_wt = tibble::tibble(clade_id = "c1", m_clade_wt = 20)
  clades_wt_ccf = tidyr::expand_grid(clade_id = "c1", time = pop$time, epistate = c("-","+")) %>%
    dplyr::mutate(ccf = 0.1)

  list(population_sizes = pop, wt_ccf = wt_ccf,
       driver_muts = driver_muts, driver_ccf = driver_ccf,
       clades_wt = clades_wt, clades_wt_ccf = clades_wt_ccf)

}

test_that("check_input_multirates accepts a well-formed tables list", {

  expect_true(check_input_multirates(tiny_tables()))

})

test_that("check_input_multirates rejects a missing ccf row", {

  bad = tiny_tables()
  bad$driver_ccf = bad$driver_ccf[-1,]

  expect_error(check_input_multirates(bad), "missing")

})

test_that("build_stan_data_multirates builds the correct time grid", {

  x = init_multirates(tiny_tables(), m_trunk = 100)

  built = build_stan_data_multirates(x, mu = 1e-7, l = 2.7e9, t_min = 0,
            ms_epi = 0, sigma_epi = 0.5, alpha_lambda = 1, beta_lambda = 1,
            alpha_n_wt = 1, beta_n_wt = 10, alpha_p_wt = 1, beta_p_wt = 10,
            min_kappa = 10, max_kappa = 1000, min_sigma_count = 0.01, max_sigma_count = 1)

  data = built$data

  expect_equal(data$n_sampling_times, 3L)
  expect_equal(data$times[1], max(c(2,5,9)))
  expect_equal(data$sampling_index[data$n_sampling_times], 1L)
  expect_equal(data$zminus, c(100,150,300))
  expect_equal(data$zplus, c(50,120,250))

})

test_that("build_stan_data_multirates produces correctly shaped arrays", {

  x = init_multirates(tiny_tables(), m_trunk = 100)

  built = build_stan_data_multirates(x, mu = 1e-7, l = 2.7e9, t_min = 0,
            ms_epi = 0, sigma_epi = 0.5, alpha_lambda = 1, beta_lambda = 1,
            alpha_n_wt = 1, beta_n_wt = 10, alpha_p_wt = 1, beta_p_wt = 10,
            min_kappa = 10, max_kappa = 1000, min_sigma_count = 0.01, max_sigma_count = 1)

  data = built$data

  expect_equal(data$N_driver, 1L)
  expect_equal(data$N_clades_wt, 1L)
  expect_equal(data$N_dc, 0L)
  expect_equal(dim(data$ccf_driver), c(1L,3L,2L))
  expect_equal(dim(data$ccf_clade_wt), c(1L,3L,2L))
  expect_equal(dim(data$ccf_wt), c(3L,2L))
  expect_equal(dim(data$m_dc), c(0L,2L))
  expect_equal(built$index_maps$driver_ids, "KRAS")

})

test_that("pivot_ccf round-trips known values", {

  long = tidyr::expand_grid(id = c("a","b"), time = c(1,2), epistate = c("-","+")) %>%
    dplyr::mutate(ccf = as.numeric(paste0(match(id, c("a","b")), match(time, c(1,2)), match(epistate, c("-","+")))))

  arr = pivot_ccf(long, "id", c("a","b"), c(1,2), epistate = TRUE)

  expect_equal(arr[1,1,1], 111)
  expect_equal(arr[2,2,2], 222)
  expect_equal(dim(arr), c(2L,2L,2L))

})

test_that("build_stan_data_multirates handles empty optional groups", {

  tabs = tiny_tables()
  tabs$driver_muts = NULL
  tabs$driver_ccf = NULL
  tabs$clades_wt = NULL
  tabs$clades_wt_ccf = NULL

  x = init_multirates(tabs, m_trunk = 100)

  built = build_stan_data_multirates(x, mu = 1e-7, l = 2.7e9, t_min = 0,
            ms_epi = 0, sigma_epi = 0.5, alpha_lambda = 1, beta_lambda = 1,
            alpha_n_wt = 1, beta_n_wt = 10, alpha_p_wt = 1, beta_p_wt = 10,
            min_kappa = 10, max_kappa = 1000, min_sigma_count = 0.01, max_sigma_count = 1)

  expect_equal(built$data$N_driver, 0L)
  expect_equal(built$data$N_clades_wt, 0L)
  expect_equal(dim(built$data$ccf_driver), c(0L,3L,2L))

})
