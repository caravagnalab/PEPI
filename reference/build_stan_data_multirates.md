# Build the Stan data list for the unified multirates model.

Turns the tables stored on a \`PEPI_Multirates\` object (see
[`init_multirates()`](https://caravagnalab.github.io/PEPI/reference/init_multirates.md))
into the exact data list required by `inst/multirates_positive_s.stan`,
plus an \`index_maps\` list that lets getters translate Stan array
indices back into ids/times.

## Usage

``` r
build_stan_data_multirates(
  x,
  mu,
  l,
  t_min,
  ms_epi,
  sigma_epi,
  alpha_lambda,
  beta_lambda,
  alpha_n_wt,
  beta_n_wt,
  alpha_p_wt,
  beta_p_wt,
  min_kappa,
  max_kappa,
  min_sigma_count,
  max_sigma_count,
  include_poisson = 1L,
  include_trunk = 1L
)
```

## Arguments

- x:

  PEPI_Multirates object.

- mu:

  Mutation rate per division per bp per allele.

- l:

  Length of the genome.

- t_min:

  Lower bound for the tmrca prior.

- ms_epi, sigma_epi:

  Lognormal prior hyperparameters for s_epi.

- alpha_lambda, beta_lambda:

  Gamma prior hyperparameters for lambda_n.

- alpha_n_wt, beta_n_wt, alpha_p_wt, beta_p_wt:

  Gamma prior hyperparameters for the wt switch rates
  omega_n_wt/omega_p_wt.

- min_kappa, max_kappa:

  Uniform prior bounds for kappa.

- min_sigma_count, max_sigma_count:

  Uniform prior bounds for sigma_count.

- include_poisson, include_trunk:

  Stan model switches (0/1).

## Value

A list with elements \`data\` (the Stan data list) and \`index_maps\`.

## Examples

``` r
sim = simulate_multirates_tree(
  sampling_times = c(3, 6, 9), tmrca = 1, t_min = 0,
  lambda_n = 1, s_epi = 0.15, omega_n_wt = 5e-3, omega_p_wt = 2e-3,
  clades_wt = list(list(id = "c1", t_clade_wt = 1.5)),
  mu = 1e-7, l = 2.7e9, kappa = 20, sigma_count = 0.1, seed = 42
)
x = init_multirates(sim$tables, m_trunk = sim$truth$m_trunk)
build_stan_data_multirates(x, mu = 1e-7, l = 2.7e9, t_min = 0,
  ms_epi = 0, sigma_epi = 0.5, alpha_lambda = 1, beta_lambda = 1,
  alpha_n_wt = 1, beta_n_wt = 10, alpha_p_wt = 1, beta_p_wt = 10,
  min_kappa = 5, max_kappa = 200, min_sigma_count = 0.01, max_sigma_count = 1)
#> $data
#> $data$N_clades_wt
#> [1] 1
#> 
#> $data$n_times
#> [1] 3
#> 
#> $data$n_sampling_times
#> [1] 3
#> 
#> $data$n_intermediate_times
#> [1] 0
#> 
#> $data$sampling_index
#> [1] 2 3 1
#> 
#> $data$intermediate_index
#> integer(0)
#> 
#> $data$N_driver
#> [1] 0
#> 
#> $data$N_dc
#> [1] 0
#> 
#> $data$N_driver_n
#> [1] 0
#> 
#> $data$N_driver_p
#> [1] 0
#> 
#> $data$include_poisson
#> [1] 1
#> 
#> $data$m_trunk
#> [1] 1035
#> 
#> $data$include_trunk
#> [1] 1
#> 
#> $data$m_driver
#> integer(0)
#> 
#> $data$ccf_driver
#> , , 1
#> 
#>      [,1] [,2] [,3]
#> 
#> , , 2
#> 
#>      [,1] [,2] [,3]
#> 
#> 
#> $data$ms_driver
#> numeric(0)
#> 
#> $data$sigma_driver
#> numeric(0)
#> 
#> $data$m_driver_n
#> integer(0)
#> 
#> $data$ccf_driver_n
#>      [,1] [,2] [,3]
#> 
#> $data$ms_driver_n
#> numeric(0)
#> 
#> $data$sigma_driver_n
#> numeric(0)
#> 
#> $data$ccf_driver_p
#>      [,1] [,2] [,3]
#> 
#> $data$ms_driver_p
#> numeric(0)
#> 
#> $data$sigma_driver_p
#> numeric(0)
#> 
#> $data$m_dc
#>      [,1] [,2]
#> 
#> $data$ccf_dc_driver
#> , , 1
#> 
#>      [,1] [,2] [,3]
#> 
#> , , 2
#> 
#>      [,1] [,2] [,3]
#> 
#> 
#> $data$ccf_dc_clade
#> , , 1
#> 
#>      [,1] [,2] [,3]
#> 
#> , , 2
#> 
#>      [,1] [,2] [,3]
#> 
#> 
#> $data$ms_dc
#> numeric(0)
#> 
#> $data$sigma_dc
#> numeric(0)
#> 
#> $data$times
#> [1] 9 3 6
#> 
#> $data$zminus
#>          3          6          9 
#>   23.03808  381.37337 8408.36760 
#> 
#> $data$zplus
#>           3           6           9 
#>   0.1621439   8.1756684 305.5005431 
#> 
#> $data$ztot
#> numeric(0)
#> 
#> $data$m_clade_wt
#> [1] 500
#> 
#> $data$ccf_clade_wt
#> , , 1
#> 
#>           [,1]      [,2]     [,3]
#> [1,] 0.2176044 0.1906469 0.266163
#> 
#> , , 2
#> 
#>           [,1]      [,2]      [,3]
#> [1,] 0.1610224 0.0726404 0.1787246
#> 
#> 
#> $data$ccf_wt
#>      [,1] [,2]
#> [1,]    1    1
#> [2,]    1    1
#> [3,]    1    1
#> 
#> $data$t_min
#> [1] 0
#> 
#> $data$ms_epi
#> [1] 0
#> 
#> $data$sigma_epi
#> [1] 0.5
#> 
#> $data$alpha_lambda
#> [1] 1
#> 
#> $data$beta_lambda
#> [1] 1
#> 
#> $data$alpha_n_wt
#> [1] 1
#> 
#> $data$beta_n_wt
#> [1] 10
#> 
#> $data$alpha_p_wt
#> [1] 1
#> 
#> $data$beta_p_wt
#> [1] 10
#> 
#> $data$alpha_n_driver
#> numeric(0)
#> 
#> $data$beta_n_driver
#> numeric(0)
#> 
#> $data$alpha_p_driver
#> numeric(0)
#> 
#> $data$beta_p_driver
#> numeric(0)
#> 
#> $data$alpha_n_dc
#> numeric(0)
#> 
#> $data$beta_n_dc
#> numeric(0)
#> 
#> $data$alpha_p_dc
#> numeric(0)
#> 
#> $data$beta_p_dc
#> numeric(0)
#> 
#> $data$mu
#> [1] 1e-07
#> 
#> $data$l
#> [1] 2.7e+09
#> 
#> $data$min_kappa
#> [1] 5
#> 
#> $data$max_kappa
#> [1] 200
#> 
#> $data$max_sigma_count
#> [1] 1
#> 
#> $data$min_sigma_count
#> [1] 0.01
#> 
#> 
#> $index_maps
#> $index_maps$sampling_times
#> [1] 3 6 9
#> 
#> $index_maps$intermediate_times
#> numeric(0)
#> 
#> $index_maps$driver_ids
#> character(0)
#> 
#> $index_maps$driver_n_ids
#> character(0)
#> 
#> $index_maps$driver_p_ids
#> character(0)
#> 
#> $index_maps$dc_ids
#> character(0)
#> 
#> $index_maps$clade_wt_ids
#> [1] "c1"
#> 
#> $index_maps$epistates
#> [1] "-" "+"
#> 
#> 
```
