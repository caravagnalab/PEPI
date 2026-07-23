# Simulate synthetic input tables for \`init_multirates()\`.

Simulate synthetic input tables for \`init_multirates()\`.

## Usage

``` r
simulate_multirates_tree(
  sampling_times,
  tmrca,
  t_min,
  lambda_n,
  s_epi,
  omega_n_wt,
  omega_p_wt,
  drivers = list(),
  driver_n = list(),
  driver_p = list(),
  dc = list(),
  clades_wt = list(),
  mu = 1e-07,
  l = 2.7e+09,
  kappa = 200,
  sigma_count = 0.1,
  seed = NULL
)
```

## Arguments

- sampling_times:

  Numeric vector of sampling times.

- tmrca:

  Time of most recent common ancestor.

- t_min:

  Time origin (lower bound of the tmrca prior / start of the wt
  lineage).

- lambda_n:

  Growth rate of "-" cells.

- s_epi:

  Fitness of "+" vs "-" wt cells.

- omega_n_wt, omega_p_wt:

  wt switch rates.

- drivers:

  List of driver specs, each a list with \`id\`, \`t_driver\`,
  \`s_driver\`, \`omega_n_driver\`, \`omega_p_driver\`, \`ms_driver\`,
  \`sigma_driver\`, \`alpha_n_driver\`, \`beta_n_driver\`,
  \`alpha_p_driver\`, \`beta_p_driver\`.

- driver_n:

  List of driver-only "-"-confined specs: \`id\`, \`t_driver_n\`,
  \`s_driver_n\`, \`ms_driver_n\`, \`sigma_driver_n\`.

- driver_p:

  List of driver-only "+"-confined specs: \`id\`, \`t_driver_p\`,
  \`s_driver_p\`, \`ms_driver_p\`, \`sigma_driver_p\`.

- dc:

  List of driver+clade combo specs: \`id\`, \`delta_t1\` (tmrca to
  driver acquisition), \`delta_t2\` (driver acquisition to clade/switch
  event), \`s_dc\`, \`omega_n_dc\`, \`omega_p_dc\`, \`ms_dc\`,
  \`sigma_dc\`, \`alpha_n_dc\`, \`beta_n_dc\`, \`alpha_p_dc\`,
  \`beta_p_dc\`.

- clades_wt:

  List of wt-background clade specs: \`id\`, \`t_clade_wt\`.

- mu:

  Mutation rate per division per bp per allele.

- l:

  Length of the genome.

- kappa:

  Concentration parameter for the beta_proportion CCF/fraction noise.

- sigma_count:

  Lognormal sd for population-size noise.

- seed:

  Optional RNG seed.

## Value

A list with \`tables\` (ready for \`init_multirates()\`) and \`truth\`
(the ground-truth parameter values used to simulate the data).

## Examples

``` r
simulate_multirates_tree(sampling_times = c(2,5,8), tmrca = 1, t_min = 0,
  lambda_n = 1, s_epi = 0.2, omega_n_wt = 1e-3, omega_p_wt = 1e-4,
  drivers = list(list(id = "KRAS", t_driver = 3, s_driver = 0.3,
    omega_n_driver = 1e-3, omega_p_driver = 1e-4, ms_driver = -0.5,
    sigma_driver = 0.5, alpha_n_driver = 1, beta_n_driver = 10,
    alpha_p_driver = 1, beta_p_driver = 10)),
  seed = 1908)
#> $tables
#> $tables$population_sizes
#> # A tibble: 3 × 4
#>    time type      zminus   zplus
#>   <dbl> <chr>      <dbl>   <dbl>
#> 1     2 sampling    8.36 0.00163
#> 2     5 sampling  165.   0.136  
#> 3     8 sampling 2981.   5.33   
#> 
#> $tables$wt_ccf
#> # A tibble: 6 × 3
#>    time epistate   ccf
#>   <dbl> <chr>    <dbl>
#> 1     2 -        1    
#> 2     2 +        1    
#> 3     5 -        0.916
#> 4     5 +        0.992
#> 5     8 -        0.955
#> 6     8 +        0.961
#> 
#> $tables$driver_muts
#> # A tibble: 1 × 8
#>   driver_id m_driver ms_driver sigma_driver alpha_n_driver beta_n_driver
#>   <chr>        <int>     <dbl>        <dbl>          <dbl>         <dbl>
#> 1 KRAS          2152      -0.5          0.5              1            10
#> # ℹ 2 more variables: alpha_p_driver <dbl>, beta_p_driver <dbl>
#> 
#> $tables$driver_ccf
#> # A tibble: 6 × 4
#>   driver_id  time epistate       ccf
#>   <chr>     <dbl> <chr>        <dbl>
#> 1 KRAS          2 -        8.85e-281
#> 2 KRAS          2 +        1.11e-312
#> 3 KRAS          5 -        5.71e-  2
#> 4 KRAS          5 +        1.68e-  2
#> 5 KRAS          8 -        3.45e-  2
#> 6 KRAS          8 +        2.62e-  2
#> 
#> $tables$driver_n_muts
#> NULL
#> 
#> $tables$driver_n_ccf
#> NULL
#> 
#> $tables$driver_p_muts
#> NULL
#> 
#> $tables$driver_p_ccf
#> NULL
#> 
#> $tables$dc_muts
#> NULL
#> 
#> $tables$dc_ccf
#> NULL
#> 
#> $tables$clades_wt
#> NULL
#> 
#> $tables$clades_wt_ccf
#> NULL
#> 
#> 
#> $truth
#> $truth$lambda_n
#> [1] 1
#> 
#> $truth$s_epi
#> [1] 0.2
#> 
#> $truth$omega_n_wt
#> [1] 0.001
#> 
#> $truth$omega_p_wt
#> [1] 1e-04
#> 
#> $truth$tmrca
#> [1] 1
#> 
#> $truth$t_min
#> [1] 0
#> 
#> $truth$m_trunk
#> [1] 1076
#> 
#> $truth$s_driver
#> KRAS 
#>  0.3 
#> 
#> $truth$s_driver_n
#> named list()
#> 
#> $truth$s_driver_p
#> named list()
#> 
#> $truth$s_dc
#> named list()
#> 
#> 
```
