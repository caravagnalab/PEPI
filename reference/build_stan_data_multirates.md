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
source(system.file("examples", "tables_example.R", package = "PEPI"))
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union
build_stan_data_multirates(x, mu = mu$mean, l = 51e06, t_min = -20,
  ms_epi = -0.8, sigma_epi = 0.9, alpha_lambda = 10, beta_lambda = 10,
  alpha_n_wt = 1, beta_n_wt = 1e4, alpha_p_wt = 0.2, beta_p_wt = 10,
  min_kappa = 3000, max_kappa = 5000, min_sigma_count = 0.05, max_sigma_count = 0.15,
  include_poisson = 1, include_trunk = 0)
#> Warning: replacing previous import ‘cli::num_ansi_colors’ by ‘crayon::num_ansi_colors’ when loading ‘easypar’
#> Warning: replacing previous import ‘S4Arrays::makeNindexFromArrayViewport’ by ‘DelayedArray::makeNindexFromArrayViewport’ when loading ‘SummarizedExperiment’
#> Warning: replacing previous import ‘Biostrings::translate’ by ‘seqinr::translate’ when loading ‘dndscv’
#> ✔ Loading mobster, 'Model-based clustering in cancer'. Support :
#> <https://caravagn.github.io/mobster/>.
#> $data
#> $data$N_clades_wt
#> [1] 0
#> 
#> $data$n_times
#> [1] 2
#> 
#> $data$n_sampling_times
#> [1] 2
#> 
#> $data$n_intermediate_times
#> [1] 0
#> 
#> $data$sampling_index
#> [1] 1 2
#> 
#> $data$intermediate_index
#> integer(0)
#> 
#> $data$N_driver
#> [1] 1
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
#> [1] 0
#> 
#> $data$include_trunk
#> [1] 0
#> 
#> $data$m_driver
#> [1] 424
#> 
#> $data$ccf_driver
#> , , 1
#> 
#>           [,1]      [,2]
#> [1,] 0.6622551 0.7733073
#> 
#> , , 2
#> 
#>           [,1]      [,2]
#> [1,] 0.6973267 0.7752791
#> 
#> 
#> $data$ms_driver
#> [1] -0.8
#> 
#> $data$sigma_driver
#> [1] 0.9
#> 
#> $data$m_driver_n
#> integer(0)
#> 
#> $data$ccf_driver_n
#>      [,1] [,2]
#> 
#> $data$ms_driver_n
#> numeric(0)
#> 
#> $data$sigma_driver_n
#> numeric(0)
#> 
#> $data$ccf_driver_p
#>      [,1] [,2]
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
#>      [,1] [,2]
#> 
#> , , 2
#> 
#>      [,1] [,2]
#> 
#> 
#> $data$ccf_dc_clade
#> , , 1
#> 
#>      [,1] [,2]
#> 
#> , , 2
#> 
#>      [,1] [,2]
#> 
#> 
#> $data$ms_dc
#> numeric(0)
#> 
#> $data$sigma_dc
#> numeric(0)
#> 
#> $data$times
#> [1] 12.44268 13.05335
#> 
#> $data$zminus
#> [1]  397729.6 1030392.6
#> 
#> $data$zplus
#> [1]  613447.4 1863968.6
#> 
#> $data$ztot
#> numeric(0)
#> 
#> $data$m_clade_wt
#> integer(0)
#> 
#> $data$ccf_clade_wt
#> , , 1
#> 
#>      [,1] [,2]
#> 
#> , , 2
#> 
#>      [,1] [,2]
#> 
#> 
#> $data$ccf_wt
#>           [,1]      [,2]
#> [1,] 0.3377449 0.3026733
#> [2,] 0.2266927 0.2247209
#> 
#> $data$t_min
#> [1] -20
#> 
#> $data$ms_epi
#> [1] -0.8
#> 
#> $data$sigma_epi
#> [1] 0.9
#> 
#> $data$alpha_lambda
#> [1] 10
#> 
#> $data$beta_lambda
#> [1] 10
#> 
#> $data$alpha_n_wt
#> [1] 1
#> 
#> $data$beta_n_wt
#> [1] 10000
#> 
#> $data$alpha_p_wt
#> [1] 0.2
#> 
#> $data$beta_p_wt
#> [1] 10
#> 
#> $data$alpha_n_driver
#> [1] 1.5
#> 
#> $data$beta_n_driver
#> [1] 4
#> 
#> $data$alpha_p_driver
#> [1] 0.3
#> 
#> $data$beta_p_driver
#> [1] 6
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
#> [1] 5.509475e-07
#> 
#> $data$l
#> [1] 5.1e+07
#> 
#> $data$min_kappa
#> [1] 3000
#> 
#> $data$max_kappa
#> [1] 5000
#> 
#> $data$max_sigma_count
#> [1] 0.15
#> 
#> $data$min_sigma_count
#> [1] 0.05
#> 
#> 
#> $index_maps
#> $index_maps$sampling_times
#> [1] 12.44268 13.05335
#> 
#> $index_maps$intermediate_times
#> numeric(0)
#> 
#> $index_maps$driver_ids
#> [1] "B"
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
#> character(0)
#> 
#> $index_maps$epistates
#> [1] "-" "+"
#> 
#> 
```
